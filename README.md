# RoLandManipulator.jl
**This work is currently under review:**
C. Stoeffler, J. Janzen, A. del Rio and H. Peters, Design Analysis of a Novel Belt-Driven Manipulator for Fast Movements, IEEE International Conferenceon Robotics and Automation, 2024, Submitted September 2023

* Be aware that the kinematic functions are written for $`\mathbf{\theta} \in \mathbb{R}^5`$, whereas the trajectory optimization works for $`\mathbf{\theta} \in \mathbb{R}^4`$

## Installing julia and installing package

* Visit https://julialang.org and download the appropriate version
* Use of `const` inside `mutable struct` requires at least julia version 1.8 here!
* A symlink can be created in `.bashrc` with `export PATH="$PATH:/path/to/<Julia directory>/bin"`
* Download the repo
* Start julia and open the package manager with `]` and run `dev --local RELATIVE_PATH_TO/RoLandManipulator`

## Using the package
The package can now be loaded in a Julia session and a variable of type `Manipulator` can be created by calling a function from a config file.
```jl
julia> using RoLandManipulator

julia> manipulator = build_old_manipulator()
```

This creates a MeshCat server for the visualization of the manipulator in its initial configuration
```jl
[ Info: Listening on: 127.0.0.1:8704, thread id: 1
┌ Info: MeshCat server started. You can open the visualizer by visiting the following URL in your browser:
└ http://127.0.0.1:8700

julia> # MORE DOC HERE!
```
### Simple kinematics
Forward and inverse kinematics can be computed, where the elbow joint possesses two inverse solutions. Note that the functions `forward_kinematics!` and `inverse_kinematics!` alter the variable of type `Manipulator`:
```jl
θ = 1 .- 2rand(5,)
C = forward_kinematics!(manipulator, θ)
update_visualization!(manipulator)
θn = inverse_kinematics(manipulator, C, sol = 2)
```

### Trajectory generation
The optimization is performed in the _joint space_, which implies the actuation torques are obtained through a mapping. Optionally the optimization can be carried out directly in the _actuation space_. 
Both options are explained below. 
```jl
# initial and final configuration
x0 = [0.17; -1.22; -0.7; 0; zeros(4)]
xf = [3.12; 1.22; 1.57; 0; zeros(4)]

# Setting weighting matrices for the optimization in the joint space. Hyperparameters are set to comman values but can be accessed through keyword arguments (see docs). 
# The optimzation with these weights has been tested and leads to a reasonable trajectory where the goal state is reached.  
op_joint = OptimizationParameters(x0, xf, 0.0001, 0.8, Q = diagm([10*ones(4);0.01*ones(4)]), R = diagm([1;1;2.5;1000]), Qf = diagm([20*1e6*ones(4);1e6*ones(4)]), 1e-8); 

X_joint, U_joint = iLQR(manipulator, op_joint);

# mapping torques into actuation space
ζ = Vector{Vector{Float64}}(fill(zeros(manipulator.m), op_joint.N-1))
for i = 1:op_joint.N-1
    ζ[i] = manipulator.S_'*U_joint[i]
end 

stp = 100   # specify a step length to smoothen curves
plot_states(op_joint.times[begin:stp:end], X_joint[begin:stp:end]) 
plot_torques(op_joint.times[begin:stp:end], U_joint[begin:stp:end], ζ[begin:stp:end]) 

# animating the trajectory (may take some time in the first run)
animate_manipulator!(manipulator, op_joint.times, X_joint)
```
Creating a trajectory in actuation space, the keyword **aspo** is set true in the `OptimizationParameters` variable. To compare resulting joint and actuation torques again the torques are mapped between the two spaces.  
```jl
# Setting weigths and hyperparameters
op_act = OptimizationParameters(x0, xf, 0.0001, 0.8, Q = diagm([ones(4);0.1ones(4)]), R = 10diagm([1;1;1;1]), Qf = diagm(1e7ones(8)), 1e-8, aspo=true); 

X_act, U_act = iLQR(manipulator, op_act);

# mapping into joint space 
τ = Vector{Vector{Float64}}(fill(zeros(op_act.m), op_act.N-1))
for i = 1:op_act.N-1
    τ[i] = manipulator.S[2:end,2:end]'\U_act[i]
end 
# plot state and torque trajectory
plot_states(op_act.times[begin:stp:end], X_act[begin:stp:end]) 
plot_torques(op_act.times[begin:stp:end], τ[begin:stp:end], U_act[begin:stp:end]) 

animate_manipulator!(manipulator, op_act.times, X_act)
```
Alternatively, a simple trajectory can be created as well
```jl
# making time steps bigger for that case
op = OptimizationParameters(x0, xf, 0.01, 0.8, Q = diagm([ones(4);0.1ones(4)]), R = 10diagm([1;1;1;1]), Qf = diagm(1e7ones(8)), 1e-8, aspo=true);

X, U = naive_trajectory(manipulator, op);

plot_states(op.times, X) 
plot_torques(op.times, U)

animate_manipulator!(manipulator, op.times, X)
```

### Eigenmodes
A plot of the frequency modes for 500 random poses is created with
```jl
julia> plt = show_frequencies(manipulator, 500)
```

## Joint space kinematics
![test](./images/arm_ik.png?raw=true "Zero pose and multiple solutions")
![test](./images/relative_transforms.png?raw=true "Relative rigid-body transforms")

## Belt kinematics of the old design

![test](./images/belt_drives.png?raw=true "Belt routing in robot arm")

Due to belts and bevel gear, there is a linear transform between joint coordinates and motor coordinates, involving the gear ratios. The underlying schematics looks as follows:
![test](./images/arm_schematics_old.png?raw=true "Schematics of the arm for the last 4 DOFs in planar depiction")

From this, one can obtain the connected graph that represents the system of linear equations

![test](./images/graph_old.png?raw=true "Unerlying graph of the manipulator")

The transfer between motor and joint coordinates can be represented in the structure matrix $\mathbf{S}^T$
```math
\mathbf{\mu} = \mathbf{S}\mathbf{\theta}
```
where $`\mathbf{\theta}`$ denotes joint coordinates and $`\mathbf{\mu}`$ motor coordinates. The matrix $`\mathbf{S}`$ is obtained from the set of circuit- and coaxial-equations and writes
```math
\mathbf{S}^T = 
        \begin{matrix}
                1/g_1 & 1 & 1 & 1 \\
                0 & 1/g_2 & 1/g_3 & 1/g_4 \\
                0 & 0 & 1/g_5 & 1/g_6 \\
                0 & 0 & 1/g_5 & -1/g_6
        \end{matrix}
```
with the following gear ratios $`g_i`$ related to schemactics and graph:
```math
\mathbf{g} = \begin{matrix}
                g_1 \\
                g_2 \\
                g_3 \\
                g_4 \\
                g_5 \\
                g_6 \\
        \end{matrix} = 
        \begin{matrix}
                N_1^5N_5^6 \\
                N_2^7N_7^8 \\
                N_3^9 \\
                N_4^{10} \\
                N_3^9N_9^{11} \\
                N_4^{10}N_{10}^{12}
        \end{matrix}
```
where the gear ratio computes from the wheel belt radi s.t. $N_j^i = r_j/r_i$ and $1/N_j^i = N_i^j$. 

## Belt kinematics of the new design

![test](./images/arm_schematics_new.png?raw=true "Schematics of the arm for the last 4 DOFs in planar depiction")

From this, one can obtain the connected graph that represents the system of linear equations

![test](./images/graph_new.png?raw=true "Unerlying graph of the manipulator")

