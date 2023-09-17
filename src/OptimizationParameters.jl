mutable struct OptimizationParameters
    
    # general parameters
    const x0::Vector                      # initial configuration in state-space
    const xf::Vector                      # final configuration in state-space
    const tf::Union{Float64, Nothing}     # final time
    const times::AbstractArray            # time stamps either as Vector or StepRange etc.
    const N::Int                          # number of steps 
    const n::Int                          # number of states
    const m::Int                          # number of controls
    const aspo::Bool                      # decides over performing the optimization in the actuation space 

    # optimal control parameters
    const Q::Matrix
    const R::Union{Matrix, Real}
    const Qf::Matrix

    const tol::Float64                    # stopping criteria for biggest term in feed-forward vector
    const i_max::Int                      # max. Iterations of forward/backward pass 

    β::Float64                            # regularization parameter
    const Φ_β::Float64                    # regularization scaling paramterer
    const β_max::Float64                  # max. regularization value
    const r_max::Float64                  # max. regularization iterations  

    const J_max::Float64                  # maximum cost allowed in rollout
    const lb::Float64                     # line search lower bound
    const ub::Float64                     # line search upper bound            
    const l_max::Int                      # maximum line search iterations

    """
    Standard constructor with predefined final time and step-size
    """
    function OptimizationParameters(x0::Vector, xf::Vector, h::Float64, tf::Float64; Q::Matrix, R::Union{Matrix, Real}, Qf::Matrix, tol::Float64, aspo=false, i_max = 50, β=0, Φ_β=2, β_max=1e8, r_max=50, J_max=1e8, lb=1e-2, ub=10, l_max=5)
        if size(Qf) != size(Q)
            throw(error("Final and intermediate cost matrix dimensions do not match. Qf ∈ ℝ^{$(size(Qf)[1]) × $(size(Qf)[2])} but Q ∈ ℝ^{$(size(Q)[1]) × $(size(Q)[2])}"))
        end 
        times = 0:h:tf
        N = length(times)
        n = size(Q)[1]
        if R isa Matrix
            m = size(R)[1]
        else
            m = 1
        end
        new(x0, xf, tf, times, N, n, m, aspo, Q, R, Qf, tol, i_max, β, Φ_β, β_max, r_max, J_max, lb, ub, l_max)
    end

    function OptimizationParameters(x0::Vector, xf::Vector, times::AbstractArray; Q::Matrix, R::Union{Matrix, Real}, Qf::Matrix, tol::Float64, aspo=false, i_max = 50, β=0, Φ_β=2, β_max=1e8, r_max=50, J_max=1e8, lb=1e-2, ub=10, l_max=5)
        if size(Qf) != size(Q)
            throw(error("Final and intermediate cost matrix dimensions do not match. Qf ∈ ℝ^{$(size(Qf)[1]) × $(size(Qf)[2])} but Q ∈ ℝ^{$(size(Q)[1]) × $(size(Q)[2])}"))
        end 

        N = length(times)
        
        if !all(isapprox.([times[k+1] - times[k] for k = 1:N-1], times[2] - times[1]))
            @info "Non-equidistant time array used for computation"
        end
        new(x0, xf, times[end], times, N, n, m, aspo, Q, R, Qf, tol, i_max, β, Φ_β, β_max, r_max, J_max, lb, ub, l_max)
    end

end
