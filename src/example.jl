drake_traj = CSV.read("./move_up_traj.csv", DataFrame); 
# parsing values to matrices 
pos = [Vector{Float64}(undef, 4) for _ in 1:nrow(drake_traj)];
tor = [Vector{Float64}(undef, 4) for _ in 1:nrow(drake_traj)];
vel = [Vector{Float64}(undef, 4) for _ in 1:nrow(drake_traj)];
for i = 1:nrow(drake_traj)
    pos[i] .=  drake_traj.joint1_pos[i], drake_traj.joint2_pos[i], drake_traj.joint3_pos[i], drake_traj.joint4_pos[i]
    tor[i] .=  drake_traj.joint1_torque[i], drake_traj.joint2_torque[i], drake_traj.joint3_torque[i], drake_traj.joint4_torque[i]
    vel[i] .=  drake_traj.joint1_vel[i], drake_traj.joint2_vel[i], drake_traj.joint3_vel[i], drake_traj.joint4_vel[i]
end 
   
# building old and new manipulator structure (needed for all robots)
old_man = build_old_manipulator(); 

######## RoLand manipulator
# these work !!
# op = OptimizationParameters([pos[1]; zeros(4)], [pos[end]; zeros(4)], 0.001, 3., diagm([10*ones(4); 0.01*ones(4)]), 1.0diagm(ones(4)), diagm([300000*ones(4);ones(4)]), 1e-8); 
# op = OptimizationParameters([pos[1]; zeros(4)], [pos[end]; zeros(4)], 0.001, 2., diagm([10*ones(4); 0.01*ones(4)]), 1.0diagm(ones(4)), diagm([300000*ones(4);ones(4)]), 1e-8); 
# op = OptimizationParameters([pos[1]; zeros(4)], [pos[end]; zeros(4)], 0.001, 0.5, diagm([10*ones(4); 0.01*ones(4)]), 1.0diagm(ones(4)), diagm([350000*ones(4);ones(4)]), 1e-8); 
# the perfect trajectory
# op = OptimizationParameters([pos[1]; zeros(4)], [pos[end]; zeros(4)], 0.001, 0.8, diagm([10*ones(4); 0.01*ones(4)]), 1.0diagm(ones(4)), diagm([550000*ones(4);ones(4)]), 1e-8); 
# op = OptimizationParameters([pos[1]; zeros(4)], [pos[end]; zeros(4)], 0.001, 0.8, diagm([10*ones(4); 0.01*ones(4)]), 1.0diagm([ones(2)*10;10;100]), diagm([550000*ones(4);ones(4)]), 1e-8); 
op = OptimizationParameters([pos[1]; zeros(4)], [pos[end]; zeros(4)], 0.001, 0.8, diagm([10*ones(4); 0.01*ones(4)]), 1.0diagm([ones(2);1;1000]), diagm([550000*ones(4);ones(4)]), 1e-8, r_max=60, β_max=1e8,i_max=50); 

X, U = iLQR(old_man, op)     
 
piece = op.βtraj[1:end,1:801]  
heatmap(1:size(piece,1),1:size(piece,2), piece,right_margin=15mm) 


data = rand(4,4)
heatmap(1:size(data,1),1:size(data,2), data)

plot_states(op.times, X) 
plot_torques(op.times, U) 

xs = [X[i][1:4] for i=1:op.N]; 
animation = MeshCat.Animation(old_man.mvis, op.times, xs);
setanimation!(old_man.mvis, animation);  

diff = xs[end]-pos[end] .|> abs .|> rad2deg 

######## acrobot
urdf = abspath(joinpath(@__DIR__, "..", "urdf/double_pendulum.urdf"))
man = Manipulator(old_man.C0, old_man.Cl, old_man.joint_axis, old_man.S, old_man.τ_limit, old_man.K_μ, urdf, n = 4, m = 1, friction = 1.);
x0 = [0., 0., 0., 0.];  
xf = [π, π, 0., 0.]; 
# op = OptimizationParameters(x0, xf, 0.001, 5., diagm([0.005*ones(2); 0.001*ones(2)]), 0.002ones(1,1), diagm(400*ones(4)), 1e-3, 0.); # the working parameters
op = OptimizationParameters(x0, xf, 0.01, 5., diagm([0.5*ones(2); 0.1*ones(2)]), 0.02ones(1,1), diagm(4000*ones(4)), 1e-3);
X, U = iLQR(man, op)          
plot_states(op.times, X)
plot_torques(op.times, U)
xs = [X[i][1:2] for i=1:op.N];
animation = MeshCat.Animation(man.mvis, op.times, xs);
setanimation!(man.mvis, animation); 

######## triple pendulum
urdf = abspath(joinpath(@__DIR__, "..", "urdf/triple_pendulum.urdf"))
man = Manipulator(old_man.C0, old_man.Cl, old_man.joint_axis, old_man.S, old_man.τ_limit, old_man.K_μ, urdf, n = 6, m = 3, friction = 0.1);
x0 = [0., 0., 0., 0., 0., 0.];
xf = [π, 0., 0., 0., 0., 0.];
op = OptimizationParameters(x0, xf, 0.01, 5., diagm([3*ones(3); 0.01*ones(3)]), 0.01diagm(ones(3)), diagm([300000*ones(3);ones(3)]), 1e-8); # these work perfect
# op = OptimizationParameters(x0, xf, 0.01, 5., diagm([1*ones(3); 0.01*ones(3)]), 0.01diagm(ones(3)), diagm(300000*ones(6)), 1e-8); # these work quite good
X, U = iLQR(man, op)        

op.βtraj
plot(op.βtraj[7,:])


plot_states(op.times, X) 
plot_torques(op.times, U) 

xs = [X[i][1:3] for i=1:op.N];
animation = MeshCat.Animation(man.mvis, op.times, xs); 
setanimation!(man.mvis, animation); 


# THAT WORKS!!!!!!! For Acrobot
# x0 = [0., 0., 0., 0.];
# xf = [[π, 0]; [0., 0.]];
# op = OptimizationParameters(x0, xf, 0.001, 5., diagm([0.005*ones(2); 0.001*ones(2)]), 0.002ones(1,1), diagm(400*ones(4)), 1e-3, 0.);
  
# mapping into actuation space
ζ = Vector{Vector{Float64}}(fill(zeros(old_man.m), op.N-1))
for i = 1:op.N-1
    ζ[i] = old_man.S[2:end,2:end]'*U[i]
end
plot_torques(op.times, U, ζ)
plot_states(op.times, X) 

# making videos... 
# MeshCat.convert_frames_to_video("meshcat.tar")  


# function control!(τ::AbstractVector, t, state::MechanismState)
#     k = Int(round(t/op.tf*(op.N-1)))
#     if k < 1
#         k = 1
#     end
#     friction = parent(velocity(state)) .* old_man.friction
#     τ .= U[k] - friction
# end

# pos0 = deepcopy(pos[500])
# state0 = MechanismState(old_man.rbd_model, pos0, zeros(4))
# t, q, v = simulate(state0, op.tf, control!);
# # t, q, v = simulate(state0, op.tf, Δt = 1e-2);
# animation = MeshCat.Animation(old_man.mvis, t, q);
# setanimation!(old_man.mvis, animation) 

# plot_states(t, parent.(q))

# conny = []
# for i = 1:10000
#     state = MechanismState(old_man.rbd_model, 2π*rand(4), zeros(4))
#     push!(conny,1/cond(mass_matrix(state)))
# end
# scatter(conny, alpha = 0.3)
