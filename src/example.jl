# initial and final configuration
x0 = [0.17; -1.22; -0.7; 0; zeros(4)]
xf = [3.12; 1.22; 1.57; 0; zeros(4)]

# building old and new manipulator structure (needed for all robots)
old_man = build_old_manipulator(); 

######## RoLand manipulator
# these work !!
op = OptimizationParameters(x0, xf, 0.001, 1.5, diagm([10*ones(4); 0.01*ones(4)]), 15.0diagm(ones(4)), diagm([30000*ones(4);ones(4)]), 1e-8); 
# op = OptimizationParameters([pos[1]; zeros(4)], [pos[end]; zeros(4)], 0.001, 2., diagm([10*ones(4); 0.01*ones(4)]), 1.0diagm(ones(4)), diagm([300000*ones(4);ones(4)]), 1e-8); 
# op = OptimizationParameters([pos[1]; zeros(4)], [pos[end]; zeros(4)], 0.001, 0.5, diagm([10*ones(4); 0.01*ones(4)]), 1.0diagm(ones(4)), diagm([350000*ones(4);ones(4)]), 1e-8); 
# the perfect trajectory
# op = OptimizationParameters([pos[1]; zeros(4)], [pos[end]; zeros(4)], 0.001, 0.8, diagm([10*ones(4); 0.01*ones(4)]), 1.0diagm(ones(4)), diagm([550000*ones(4);ones(4)]), 1e-8); 
# op = OptimizationParameters([pos[1]; zeros(4)], [pos[end]; zeros(4)], 0.001, 0.8, diagm([10*ones(4); 0.01*ones(4)]), 1.0diagm([ones(2)*10;10;100]), diagm([550000*ones(4);ones(4)]), 1e-8); 
op = OptimizationParameters(x0, xf, 0.01, 1.8, diagm([10*ones(4); 0.01*ones(4)]), 1.0diagm([ones(2);1;1000]), diagm([550000*ones(4);ones(4)]), 1e-8, r_max=60, β_max=1e8,i_max=50); 


X, U = naive_trajectory(old_man, op);

xs = [X[i][1:4] for i=1:op.N]; 
animation = MeshCat.Animation(old_man.mvis, op.times, xs);
setanimation!(old_man.mvis, animation); 

X, U = iLQR(old_man, op);

plot_states(op.times, X) 
plot_torques(op.times, U) 


diff = xs[end]-pos[end] .|> abs .|> rad2deg 

######## acrobot
urdf = abspath(joinpath(@__DIR__, "..", "urdf/double_pendulum.urdf"))
man = Manipulator(old_man.C0, old_man.Cl, old_man.joint_axis, old_man.S, old_man.τ_limit, old_man.K_μ, urdf, n = 4, m = 1, friction = 1.);
x0 = [0., 0., 0., 0.];  
xf = [π, π, 0., 0.]; 
# op = OptimizationParameters(x0, xf, 0.001, 5., diagm([0.005*ones(2); 0.001*ones(2)]), 0.002ones(1,1), diagm(400*ones(4)), 1e-3, 0.); # the working parameters
op = OptimizationParameters(x0, xf, 0.01, 5., diagm([0.5*ones(2); 0.1*ones(2)]), 0.02ones(1,1), diagm(4000*ones(4)), 1e-3);

X, U = iLQR(man, op);
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
op = OptimizationParameters(x0, xf, 0.01, 5., diagm([3*ones(3); 0.01*ones(3)]), 0.01diagm(ones(3)), diagm([300000*ones(3);ones(3)]), 1e-8); # these work perfectly
# op = OptimizationParameters(x0, xf, 0.01, 5., diagm([1*ones(3); 0.01*ones(3)]), 0.01diagm(ones(3)), diagm(300000*ones(6)), 1e-8); # these work quite good
X, U = iLQR(man, op);

plot_states(op.times, X) 
plot_torques(op.times, U) 

xs = [X[i][1:3] for i=1:op.N];
animation = MeshCat.Animation(man.mvis, op.times, xs); 
setanimation!(man.mvis, animation); 
