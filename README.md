# RoLandArm.jl

* Be aware that the kinematic functions are written for $`\mathbf{\theta} \in \mathbb{R}^5`$, whereas the trajectory optimization works for $`\mathbf{\theta} \in \mathbb{R}^4`$

## Installing julia and installing package

* Visit https://julialang.org and download the appropriate version
* Use of `const` inside `mutable struct` requires julia version 1.8 here!
* A symlink can be created in `.bashrc` with `export PATH="$PATH:/path/to/<Julia directory>/bin"`
* Download the repo
* Start julia and open the package manager with `]` and run `dev --local RELATIVE_PATH_TO/RoLandArm`

## Using the package
The package can now be loaded in a Julia session and an object of type `Manipulator` can be created
```jl
julia> using RoLandManipulator

julia> # MORE DOC HERE!
```

Depiction of the manipulator in its initial configuration can be done
```jl
julia> # MORE DOC HERE!
```

Forward and inverse kinematics can be computed, where the elbow joint possesses two inverse solutions. Note that the function `forward_kinematics!` alters the Manipulator object and returns the end-effector frame, whereas `inverse_kinematics_analytical` simply returns the joint angles given from a end-effector pose.
```jl
θ = 1 .- 2rand(5,)
C = forward_kinematics!(man, θ)
show_manipulator(man)
θn = inverse_kinematics_analytical(man, C, sol = 2)
```
## Joint space kinematics
![test](./images/arm_ik.png?raw=true "Zero pose and multiple solutions")
![test](./images/relative_transforms.png?raw=true "Relative rigid-body transforms")

## Belt kinematics of the old design

![test](./images/belt_drives.png?raw=true "Belt routing in robot arm")

Due to belts and bevel gear, there is a linear transform between joint coordinates and motor coordinates, involving the gear ratios - see also `src/manipulator_kinematics.jl -> motor2joint()` The underlying schematics looks as follows:
![test](./images/arm_schematics.png?raw=true "Schematics of the arm for the last 4 DOFs in planar depiction")

From this, one can obtain the connected graph that represents the system of linear equations

![test](./images/graph.png?raw=true "Unerlying graph of the manipulator")

The transfer between motor and joint coordinates can be represented in the structure matrix $\mathbf{A}$
```math
\mathbf{\theta} = \mathbf{A}\mathbf{\mu}
```
where $`\mathbf{\theta}`$ denotes joint coordinates and $`\mathbf{\mu}`$ motor coordinates. However, the first joint $`\theta_1`$ relates to the prismatic joint and no coupling to other motors occurs. The matrix $`\mathbf{A}`$ is obtained from the set of circuit- and coaxial-equations and writes
```math
\mathbf{A} = 
\begin{bmatrix}
g_1 & 0 & 0 & 0 & 0 \\
0 & g_2 & 0 & 0 & 0 \\
0 & -g_2g_3 & g_3 & 0 & 0 \\
0 & 0.5(g_4 + g_7) & -0.5(g_5 + g_8) & 0.5g_6 & 0.5g_9 \\
0 & 0.5(g_4 - g_7) & -0.5(g_5 - g_8) & 0.5g_6 & -0.5g_9 \\
\end{bmatrix}
```
with the following gear ratios $`g_i`$ related to schemactics and graph:
```math
\begin{matrix}
g_2 = & N_1^5N_5^6 \\
g_3 = & N_2^7N_7^8 \\
g_4 = & (N_5^6N_1^5N_7^8N_2^7 - N_5^6N_1^5N_3^9)N_9^{11} \\
g_5 = & N_7^8N_2^7N_9^{11} \\
g_6 = & N_3^9N_9^{11} \\
g_7 = & (N_5^6N_1^5N_7^8N_2^7 - N_5^6N_1^5N_4^{10})N_{10}^{12} \\
g_8 = & N_7^8N_2^7N_{10}^{12} \\
g_9 = & N_4^{10}N_{10}^{12}
\end{matrix}
```
where the gear ratio computes from the wheel belt radi s.t. $N_j^i = r_j/r_i$ and $1/N_j^i = N_i^j$. One can see, that if the transmissions in the entire power line of the bevel gear are equal, such that $N_3^9 = N_4^{10}$ and $N_9^{11} = N_{10}^{12}$, it turns out that $g_4 = g_7$, $g_5 = g_8$ and $g_6 = g_9$. The structure matrix then reduces to
```math
\mathbf{A} = 
\begin{bmatrix}
g_1 & 0 & 0 & 0 & 0 \\
0 & g_2 & 0 & 0 & 0 \\
0 & -g_2g_3 & g_3 & 0 & 0 \\
0 & g_4 & -g_5  & 0.5g_6 & 0.5g_6 \\
0 & 0 & 0 & 0.5g_6 & -0.5g_6 \\
\end{bmatrix}
```

Looking at the time derivative of $\mathbf{\theta}$, one gets
```math
\dot{\mathbf{\theta}} = \mathbf{A}\dot{\mathbf{\mu}} + \underbrace{\dot{\mathbf{A}}\mathbf{\mu}}_{=0}
```
The power from virtual torque $`\mathbf{\tau}_v`$ in joint space must be equal to the power in the motors:
```math
\mathbf{\tau}_v^T\dot{\mathbf{\theta}} = \mathbf{\tau}^T\dot{\mathbf{\mu}}
```
```math
\mathbf{\tau}_v^T\mathbf{A}\dot{\mathbf{\mu}} = \mathbf{\tau}^T\dot{\mathbf{\mu}}
```
what can only hold true (no zeros in $`\dot{\mathbf{\mu}}`$ assumed) if
```math
\mathbf{\tau} = \mathbf{A}^T\mathbf{\tau}_v
```
and shows that between motors and joints, there is a linear constant map, given by the motor gear ratios and belt transmissions. The available torque on joint level then writes
```math
\mathbf{\tau}_v = \mathbf{A}^{-T}\mathbf{\tau}
```
And in detail
```math
\mathbf{\tau}_v = 
\begin{bmatrix}
1/g_1 & 0 & 0 & 0 & 0 \\
0 & 1/g_2 & 1 & (g_2g_5 - g_4)/(g_2g_6) & (g_2g_5 - g_4)/(g_2g_6) \\
0 & 0 & 1/g_3 & g_5/(g_3g_6) & g_5/(g_3g_6) \\
0 & 0 & 0 & 1/g_6 & 1/g_6 \\
0 & 0 & 0 & 1/g_6 & -1/g_6 \\
\end{bmatrix}
\mathbf{\tau}
```
what shows that the first revolute joint in the base $`\theta_2`$ obtains the summed torque of all subsequent motors.

## Belt kinematics of the new design

![test](./images/belt_drives.png?raw=true "Belt routing in robot arm")

Due to belts and bevel gear, there is a linear transform between joint coordinates and motor coordinates, involving the gear ratios - see also `src/manipulator_kinematics.jl -> motor2joint()` The underlying schematics looks as follows:

![test](./images/arm_schematics_new.png?raw=true "Schematics of the arm for the last 4 DOFs in planar depiction")

From this, one can obtain the connected graph that represents the system of linear equations

![test](./images/graph_new.png?raw=true "Unerlying graph of the manipulator")

For the new system, the structure matrix $\mathbf{A}$ writes
```math
\mathbf{A} = 
\begin{bmatrix}
g_1 & 0 & 0 & 0 & 0 \\
0 & g_2 & 0 & 0 & 0 \\
0 & -g_2N_{15}^8 & g_3N_{15}^8 & 0 & 0 \\
0 & 0.5(g_4 + g_7) & -0.5(g_5 + g_8) & 0.5g_6 & 0.5g_9 \\
0 & 0.5(g_4 - g_7) & -0.5(g_5 - g_8) & 0.5g_6 & -0.5g_9 \\
\end{bmatrix}
```
where $N_{15}^8 = 1$ was picked to have a parallel moving second link. The remaining gear ratios then write
```math
\begin{matrix}
g_2 = & N_1^5N_5^6 \\
g_3 = & N_2^7N_7^{15} \\
g_4 = & (N_5^6N_1^5N_{15}^8 - N_5^6N_1^5N_3^9)N_9^{11} \\
g_5 = & N_2^7N_7^{15}N_{15}^8N_9^{11} \\
g_6 = & N_3^9N_9^{11} \\
g_7 = & (N_5^6N_1^5N_{15}^8 - N_5^6N_1^5N_4^{10})N_{10}^{12} \\
g_8 = & N_2^7N_7^{15}N_{15}^8N_{10}^{12} \\
g_9 = & N_4^{10}N_{10}^{12}
\end{matrix}
```
Again, by symmetric considrations, such as $N_3^9 = N_4^{10}$ and $N_9^{11} = N_{10}^{12}$, it turns out that $g_4 = g_7$, $g_5 = g_8$ and $g_6 = g_9$. The structure matrix again reduces to
```math
\mathbf{A} = 
\begin{bmatrix}
g_1 & 0 & 0 & 0 & 0 \\
0 & g_2 & 0 & 0 & 0 \\
0 & -g_2 & g_3 & 0 & 0 \\
0 & g_4 & -g_5  & 0.5g_6 & 0.5g_6 \\
0 & 0 & 0 & 0.5g_6 & -0.5g_6 \\
\end{bmatrix}
```

The torque mapping between motor and joint space then writes
```math
\mathbf{\tau}_v = 
\begin{bmatrix}
1/g_1 & 0 & 0 & 0 & 0 \\
0 & 1/g_2 & 1/g_3 & (g_2g_5 - g_3g_4)/(g_2g_3g_6) & (g_2g_5 - g_3g_4)/(g_2g_3g_6) \\
0 & 0 & 1/g_3 & g_5/(g_3g_6) & g_5/(g_3g_6) \\
0 & 0 & 0 & 1/g_6 & 1/g_6 \\
0 & 0 & 0 & 1/g_6 & -1/g_6 \\
\end{bmatrix}
\mathbf{\tau}
```
what shows that the first revolute joint in the base $`\theta_2`$ obtains the summed torque of all subsequent motors.
