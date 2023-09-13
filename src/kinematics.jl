"""
Forward kinematics expressed in Lie group formalism relative to robots base frame
man:        manipulator geometry and configuration
θ:          5 joint angles
"""
function forward_kinematics(man::Manipulator, θ::AbstractArray)
    
    @assert length(θ) == 5 "Wrong joint angle input"
    
    conf = deepcopy(man.Cl[2])

    # POM formula - first config does not change
    for (i, C) in enumerate(man.Cl[2:end])
        conf = conf*exp_map(man.joint_axis[i]*θ[i])*C
    end

    # returning last frame - end-effector
    return conf
end

"""
Forward kinematics expressed in Lie group formalism relative to robots base frame - alters input
    man:        manipulator geometry and configuration
    θ:          5 joint angles
"""
function forward_kinematics!(man::Manipulator, θ::AbstractArray)
    
    @assert length(θ) == 5 "Wrong joint angle input"

    man.θ = θ
 
    # POM formula - first config does not change
    for (i, C) in enumerate(man.Cl[2:end])
        man.config[i+1] = man.config[i]*exp_map(man.joint_axis[i]*θ[i])*C
    end

    # returning last frame - end-effector
    return deepcopy(man.config[end])

end

"""
Inverse kinematics relative to robots COM
    man:        manipulator geometry and configuration
    C:          configuration in task space relative to COM
    sol:        inverse solution of elbow: 1 or 2
"""
function inverse_kinematics(man::Manipulator, C::SE3; sol::Int = 1)
    
    θ = fill(NaN, 5)
    # checking whether rotated x-vector lies in y-z-plane - kinematic constraint
    if isapprox(C.R.mat[1,1], 0, atol = 1e-7) 
        # extracting geometric information from relative transforms
        ee = abs(man.Cl[end].r[1])
        il = abs(man.Cl[3].r[3])
        ol = abs(man.Cl[4].r[3])

        # x-position can be directly deduced from second frame
        θ[1] = C.r[1] - man.config[1].r[1]

        # x-vector of end-effector is deemed to lie in parallel plane to global y-z-plane
        # γ is angle between end-effector x-axis and global y-axis
        # the sign component ensures to be in correct quadrant
        γ = acos(dot([0, 1, 0], C.R.mat[:,1]))*sign(C.R.mat[3,1])
        # projecting point of 4th joint
        y4 = C.r[2] - cos(γ)*ee
        z4 = C.r[3] - sin(γ)*ee
        # solving traditional 2R problem from shoulder base (Cl[1])
        Δy = y4 - man.Cl[1].r[2]; Δz = z4 - man.Cl[1].r[3]
        ϕ = atan(Δz, Δy)
        d = √(Δy^2 + Δz^2)
        α = uacos((il^2 + d^2 - ol^2)/(2*il*d))
        β = uacos((il^2 + ol^2 - Δy^2 - Δz^2)/(2*il*ol))
        if sol == 1
            θ[2:3] .= [ϕ - α - 3π/2, π - β]
        elseif sol == 2
            θ[2:3] .= [ϕ + α - 3π/2, β - π]
        else
            @assert false "sol must be of 1 or 2"
        end

        # chaining up
        θ[4] = π/2 + γ -θ[2] - θ[3]

        # last rotation in wrist is angle between global x-axis and end-effector z-axis
        # the sign component ensures to be in correct quadrant
        θ[5] = acos(dot([1, 0, 0], C.R.mat[:,3]))*sign(C.R.mat[1,2])

        θ[2:end] .= θ[2:end] .% 2π
    end

    return θ
end

"""
Inverse kinematics relative to robots COM - alters input
    man:        manipulator geometry and configuration
    C:          configuration in task space relative to COM
    sol:        inverse solution of elbow: 1 or 2
"""
function inverse_kinematics!(man::Manipulator, C::SE3; sol::Int = 1)
    forward_kinematics!(man, inverse_kinematics(man, C, sol = sol))
    return nothing
end

"""
Inverse kinematics relative to robots COM
    man:        manipulator geometry and configuration
    C:          configuration in task space relative to COM
    θ:          initial guess for Newton scheme
    ωϵ:         error threshold for rotation
    uϵ:         error threshold for translation
"""
function inverse_kinematics_numerical(man::Manipulator, C::SE3; θ::AbstractArray = zeros(5,), ωϵ::Real = 1e-6, uϵ::Real = 1e-6)
    @assert length(θ) == 5 "Initial values are vector of length 5"

    # creating copy of manipulator for iteration
    mand = deepcopy(man)

    Sd, θd = log_map(inv(forward_kinematics!(mand, θ))*C)

    i = 0
    while norm(Sd.ω) > ωϵ #|| norm(Sd.u) > uϵ # interestingly the translation criterium leads to non-convergence
        θ += pinv(body_jacobian(mand))*Sd.vec
        Sd, θd = log_map(inv(forward_kinematics!(mand, θ))*C)
        i += 1
        if i == 50
            println("No solution found")
            break
        end
    end
    println("IK needed $(i) iterations")

    θ[2:end] .= θ[2:end] .% 2π
    
    return θ
end

"""
Relates end-effector twist with joint velocities, expressed in robot frame
"""
function space_jacobian(man::Manipulator)
    
    J = Matrix{Real}(undef, (6,5))
    J[:,1] .= man.joint_axis[1].vec

    for i = 1:length(man.joint_axis)-1
        J[:, i+1] .= Ad(man.config[i+1])*man.joint_axis[i+1].vec
    end

    return J
end

"""
Gives the manipulator Jacobian w.r.t. the end-effector fixed frame
"""
function body_jacobian(man::Manipulator)
    return Matrix{Real}(Ad(inv(man.config[end]))*space_jacobian(man))
end

"""
Returning NaN when out of bounds for input of acos
==========
    x:       input for function [ ]
"""
function uacos(x::Number)
    if abs(x) > 1.0
        return NaN
    else
        return acos(x)
    end
end

"""
Returns proper location of index in array
"""
function whereis(ar::AbstractArray, a::Any)
    sol = Vector{Int}(undef, 0)
    for (i, val) in enumerate(ar)
        if val == a
            push!(sol, i)
        end
    end
    return sol
end
