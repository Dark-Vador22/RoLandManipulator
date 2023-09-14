"""
Defines geometry and carries the URDF model of the RoLand manipulator. 
Kinematic configurations and stiffness model parameters are stored.

Whereas overall kinematics are considered 5DOF, dynamics of the arm are 4DOF without prismatic joint in the shoulder.
"""
mutable struct Manipulator
    
    # geometric parameters
    const C0::SE3                   # robot center relative to other world frame
    const Cl::Vector{SE3}           # link transformations

    # joint space quantities
    const joint_axis::Vector{se3}   # joint screws in robot frame
    const K_θ::Matrix{<:Real}       # joint space stiffness matrixj 
    config::Vector{SE3}             # configuration of joint frames       
    θ::Vector{<:Real}               # joint space vector    
    
    # actuation space quantities
    const S::Matrix{<:Real}         # transpose of structure matrix: relates actuation and joint space
    const S_::Matrix{<:Real}        # inverse of S
    const τ_limit::Float64          # actuator torque limits
    const K_μ::Matrix{<:Real}       # actuation space stiffness matrix
    μ::Vector{<:Real}               # actuation space vector
    
    # dynamic model - only planar and therefore 4DOF
    rbd_model::Union{Mechanism, Nothing}
    statecache
    dyncache::DynamicsResultCache{Float64}
    n::Int                          # number of states
    m::Int                          # number of controls
    A::Matrix{<:Real}               # linearized dynamics after state
    B::Matrix{<:Real}               # linearized dynamics after controls
    friction::Vector{Float64}       # friction coefficients

    # visualization
    vis::Visualizer
    mvis::MechanismVisualizer

    """
    Standard constructor - geometry of manipulator is defined inside
    """
    function Manipulator(C0::SE3, Cl::Vector{SE3}, joint_axis::Vector{se3}, S::Matrix{<:Real}, τ_limit::Float64, K_μ::Matrix{<:Real}, urdf::String = "default"; n::Int = 8, m::Int = 4, friction::Float64 = 0.)
        
        # checking input
        DOF = length(joint_axis)
        @assert length(Cl) == DOF + 1 "Number of relative transforms and joints inconsistent"
        @assert size(S) == (DOF-1, DOF-1) "Wrong dimension of structure matrix - has to be $(DOF-1)×$(DOF-1)"
        DOF == 5 ? nothing : @warn "Different DOFs than expected - kinematic functions do not work!"

        rbd_model = parse_urdf(urdf)
        
        # THIS DOES NOT WORK AT THE MOMENT UNDFER julia 1.9.2 !!!
        # try 
        #     rbd_model = parse_urdf(urdf)
        # catch
        # else
        #     rbd_model = nothing
        #     @warn "No URDF model was loaded - dynamic functions do not work!"
        # end

        # defining absolute configurations from link transforms, assuming zero pose
        config = Vector{SE3}(undef, length(Cl))
        
        config[1] = Cl[1]
        for (i, C) in enumerate(Cl[2:end])
            config[i+1] = config[i]*C
        end

        # zero pose in both spaces
        θ = zeros(DOF)
        μ = zeros(DOF)

        # mapping stiffness from actuation space to joint space and reducing the first prismatic joint
        S_ = inv(S)
        K_θ = S'*K_μ*S

        # caches for dynamics from RigidBodyDynamics
        statecache = StateCache(rbd_model)
        dyncache = DynamicsResultCache(rbd_model)

        # creating visualizer type and rendering it
        vis = Visualizer()
        render(vis)
        # create mechanism visualizer object
        mvis = MechanismVisualizer(rbd_model, URDFVisuals(urdf), vis) 

        new(C0, Cl, joint_axis, K_θ, config, θ, S, S_, τ_limit, K_μ, μ, rbd_model, statecache, dyncache, n, m, fill(NaN, (n, n)), fill(NaN, (n, m)), friction*ones(m), vis, mvis);
    
    end
end
