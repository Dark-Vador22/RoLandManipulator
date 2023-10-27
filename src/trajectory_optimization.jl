"""
Creates a naive control trajectory in three steps:
    1. Linearly interpolating joint space coordinates for configurations between x0 and xf
    2. Assuming a given maximal acceleration until half distance, followed by maximal deceleration 
    3. Using inverse dynamics
    Optionally using forward simulation to obtain new states, otherwise take ideal states
"""
function naive_trajectory(man::Manipulator, op::OptimizationParameters; sim_forward::Bool = false)
    
    X = Vector{Vector{Float64}}(fill(zeros(man.n), op.N))
    τ = Vector{Vector{Float64}}(fill(zeros(man.m), op.N-1))
    X[1] .= op.x0

    # interpolation of joint coordinates 
    Δθ = op.xf[1:man.m] - op.x0[1:man.m]

    # assuming "bang-bang" acceleration vector
    a = 4Δθ/op.tf^2

    # fake_state is used to get the acceleration vector as SegmentedVector, which is not sufficiently documented to be created alone
    fake_state = MechanismState(man.rbd_model, zeros(man.m), a)
    fake_state_ = MechanismState(man.rbd_model, zeros(man.m), -a)
    
    for i = 1:op.N-1
        t = op.times[i]
        if t <= op.tf/2
            pos = 0.5*t^2*a + op.x0[1:man.m]
            vel = t*a
            current_state = MechanismState(man.rbd_model, pos, vel)
            if op.aspo
                τ[i] = man.S_'*inverse_dynamics(current_state, fake_state.v)
            else 
                τ[i] = inverse_dynamics(current_state, fake_state.v)
            end 
            X[i] = [pos; vel]
        else
            pos = op.x0[1:man.m] - a/4*op.tf^2 + a*(t*op.tf - 0.5t^2)
            vel = a*(op.tf-t)
            current_state = MechanismState(man.rbd_model, pos, vel) 
            if op.aspo
                τ[i] = man.S_'*inverse_dynamics(current_state, fake_state_.v)
            else 
                τ[i] = inverse_dynamics(current_state, fake_state_.v)
            end 
            X[i] = [pos; vel]
        end
    end
    X[end] = op.xf

    # forward simulation with torque inputs instead of theoretical states
    if sim_forward
        for k = 1:op.N-1
            X[k+1] = rk4((x, u) -> dynamics(man, op.aspo, x, u), X[k], τ[k], op.times[k+1] - op.times[k])
        end
    end
 
    return X, τ
end

"""
Uses the discrete dynamics (under RK4 with zero-order hold) to obtain the linear transforms 'A' and 'B'
Here, everything happens in joint-space
    x:          state variables that contain position and velocity
    u:          controls that contain the joint torques
"""
function linearize_model!(man::Manipulator, x::AbstractVector, u::Union{AbstractVector, Real}, h::Float64; aspo=false)
    discrete_dynamics(x, u) = rk4((x, u) -> dynamics(man, aspo, x, u), x, u, h)

    ForwardDiff.jacobian!(man.A, x_ -> discrete_dynamics(x_, u), x)
    ForwardDiff.jacobian!(man.B, u_ -> discrete_dynamics(x, u_), u)

end

"""
Returns the additive cost plus final cost of an entire trajectory wrt. the goal state
"""
function trajectory_cost(X, U, op::OptimizationParameters)
    J = 0.0
    
    # stage cost
    for k = 1:(op.N-1)
        J += 0.5*((X[k] - op.xf)'*op.Q*(X[k] - op.xf)) + 0.5*(U[k]'*op.R*U[k])
    end

    # terminal cost
    J += 0.5*((X[end] - op.xf)'*op.Qf*(X[end] - op.xf))

    return J
end

function backward_pass(X, U, man::Manipulator, op::OptimizationParameters; fpf = false, verbose=false)
    # preallocating feedback gains and feed-forward terms
    K = Vector{Matrix{Float64}}(fill(zeros(op.m, op.n), op.N-1))
    d = Vector{Vector{Float64}}(fill(zeros(op.m), op.N-1))
    
    # storing two terms for expected cost change
    ΔVu = 0.
    ΔVuu = 0.

    # optimal terminal cost to go
    p = op.Qf*(X[end] - op.xf)
    P = op.Qf

    for k = (op.N-1):-1:1
        # getting A and B matrices
        linearize_model!(man, X[k], U[k], op.times[k+1] - op.times[k], aspo=op.aspo)
        
        # computing linearized action-value function
        Qx = op.Q*(X[k] - op.xf) + man.A'*p
        Qu = op.R*U[k] + man.B'*p

        # initial value for regularization parameter
        op.β = 0
        
        # increase if previous forward pass fails
        if fpf 
            op.β = op.Φ_β
        end

        # avoid unnecessary matrix operations
        if iszero(op.β)
            Iβ = 0.
        else 
            Iβ = man.B'*op.β*I*man.B
        end 

        Quu = op.R + man.B'*P*man.B .+ Iβ 
        Qxx = op.Q + man.A'*P*man.A
        Qux = man.B'*P*man.A 

        c = 0
        while !isposdef(Quu)
            
            op.β = (op.β+1)*op.Φ_β
            op.β > op.β_max && (op.β = op.β_max)

            Quu += man.B'*op.β*I*man.B 
            if c > op.r_max
                if verbose
                    @warn "Regularization failed at k = $k"
                end 
                break
            end
            c += 1
        end
        # op.β > op.β_max && (op.β = op.β_max)
        
        # if partial hessian is rank deficient, reinitialize backward pass with increased regularization
        if !isinvertible(Quu)
            return K, d, ΔVu, ΔVuu, true
        end

        d[k] = -Quu\Qu
        K[k] = -Quu\Qux

        # expected change in cost-to-go
        ΔVu += d[k]'*Qu
        ΔVuu += 0.5*d[k]'*Quu*d[k]

        # updating p and P
        p = Qx + K[k]'*Quu*d[k] + K[k]'*Qu + Qux'*d[k]
        P = Qxx + K[k]'*Quu*K[k] + K[k]'*Qux + Qux'*K[k]
    end

    return K, d, ΔVu, ΔVuu, false 
end

function forward_pass(X, U, J, K, d, ΔVu, ΔVuu, man::Manipulator, op::OptimizationParameters; bpf = false, verbose=false)
    # skip execution and directly jump to backward pass if previous backward pass was stopped due to rank deficiency
    if bpf 
        return X, U, J, true
    end 
    # capture success of backward pass  
    fpf = false

    # initializing problem
    Xn = deepcopy(X)
    Un = deepcopy(U)
    α = 1.
     
    # line search
    count = 0
    Jn = 0.
    while true
        for k = 1:op.N-1
            Un[k] = U[k] + K[k]*(Xn[k] - X[k]) + α*d[k]
            Xn[k+1] = rk4((x, u) -> dynamics(man, op.aspo, x, u), Xn[k], Un[k], op.times[k+1] - op.times[k])
        end
        Jn = trajectory_cost(Xn, Un, op)
        # skip line search if cost limit surpassed
        if Jn > op.J_max
            fpf = true
            break 
        end 
        # new line search 
        z = -(J - Jn)/(α*ΔVu + α^2*ΔVuu)
        if op.lb < z < op.ub
            break
        elseif count > op.l_max
            if verbose
                @warn "Line search does not converge!"
            end 
            fpf = true 
            break
        end
        α = 0.5*α
        count += 1
    end
    return Xn, Un, Jn, fpf
end 

function iLQR(man::Manipulator, op::OptimizationParameters, inverse_dyn::Bool = false; verbose = false)
    # checking first the conformity of Manipulator and OptimizationParameters
    if op.n != man.n
        throw(error("State dimensions do not match cost matrix dimensions. x ∈ ℝ^{$(op.n) × 1} but Q ∈ ℝ^{$(man.n) × $(man.n)} "))
    end 
    if op.m != man.m
        throw(error("Control vector dimensions do not match cost matrix dimensions. u ∈ ℝ^{$(op.m) × 1} but R ∈ ℝ^{$(man.m) × $(man.m)} "))
    end 
    # checking reachability of configurations
    # checking x-value of zero and goal configuration and throwing warning, if not identical
    if inverse_dyn
        X, U = initialize_trajectory(man, op)
    else
        # initializing with random controls
        X = [op.x0 for _ = 1:op.N]
        U = [randn(op.m) for _ = 1:op.N-1].*1
        # initial rollout to have consistent values
        for k = 1:op.N-1
            X[k+1] = rk4((x, u) -> dynamics(man, op.aspo, x, u), X[k], U[k], op.times[k+1] - op.times[k])
        end
    end

    return iLQR(X, U, man, op, verbose=verbose)
end

function iLQR(X::AbstractArray, U::AbstractArray, man::Manipulator, op::OptimizationParameters; verbose = false)
    # getting feedbak and feed-forward terms from initial trajectory
    K, d, ΔVu, ΔVuu = backward_pass(X, U, man, op)
    J = trajectory_cost(X, U, op)
    c = 1
    bpf = false
    while maximum(abs.(collect(Iterators.flatten(d)))) > op.tol
        X, U, J, fpf = forward_pass(X, U, J, K, d, ΔVu, ΔVuu, man, op, bpf=bpf,verbose=verbose)
        K, d, ΔVu, ΔVuu, bpf = backward_pass(X, U, man, op, fpf=fpf,verbose=verbose);
        c += 1
        println("iteration: $c")
        if c >= op.i_max
            break
        end
    end
    return X, U
end


"""
4th order Runge-Kutta method with state and control vector 'x', 'u' and time-step 'h'
"""
function rk4(f::Any, x::AbstractVector, u::Union{AbstractVector, Real}, h::Float64)
   
    k1 = f(x, u)
    k2 = f(x + 0.5h*k1, u)
    k3 = f(x + 0.5h*k2, u)
    k4 = f(x + h*k3, u)
    
    return x + h/6*(k1 + 2k2 + 2k3 + k4)
end

"""
Accurate friction model considering coulomb and viscous friction 
"""
function friction_model(q̇)
    fs = 1.0        # static coulomb  
    fv = 0.5        # viscous friction 
    τ_f = fs*sign.(q̇) + fv*q̇
    return τ_f
end   

"""
Continuous dynamics of the manipulator from RigidBodyDynamics
"""
function dynamics(man::Manipulator, aspo::Bool, x::AbstractVector{T1}, u::AbstractVector{T2}) where {T1,T2}
    T = promote_type(T1, T2)
    state = man.statecache[T]
    res = man.dyncache[T]

    copyto!(state, x)

    # TODO: decide if artificially making dynamics underactuated should still be an option
    if size(u) == (1,)
        u_ = [0., u[1]]
    else
        if aspo 
            u_ = man.S'*u
        else 
            u_ = u
        end 
    end

    # friction model only partly velocity dependent 
    friction = friction_model.(parent(velocity(state)))
    
    dynamics!(res, state, u_ - friction)
    q̇ = res.q̇
    v̇ = res.v̇
    return [q̇; v̇]
end 

"""
Gives the ranges to acces the different parts of the dynamic models
"""
function get_partition(man::Manipulator)
    n, m = man.n, man.m
    return  1:n, n .+ (1:m), n+m .+ (1:n)
end


