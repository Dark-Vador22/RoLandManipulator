module RoLandManipulator

using LinearAlgebra
export diagm

using ForwardDiff
using StaticArrays

# data stuff
using LightXML
using JLD2, FileIO
using ImageIO
using Measures 
using CSV
using DataFrames

# visualization
using Plots
using MeshCat
using MeshCatMechanisms
using RigidBodyDynamics

# we are going to overload these ones
import Base: +, -, *, /
import LinearAlgebra: det, inv

# package internal helper functions
include("./utils.jl")

# basic types and functions for representation by Lie groups
include("./geometric_types.jl")
export so3, SO3, se3, SE3
include("./geometric_functions.jl")
export exp_map, dexp, dexp_, log_map, Ad, ad

# core type that stores all data
include("./Manipulator.jl")
export Manipulator

# hyperparameters for OptimizationParameters
include("./OptimizationParameters.jl")
export OptimizationParameters

# basic kinematic functions
include("./kinematics.jl")
export inverse_kinematics!, inverse_kinematics, forward_kinematics!, forward_kinematics

include("./stiffness.jl")

include("./trajectory_optimization.jl")
export dynamics, naive_trajectory, iLQR

include("./visualization.jl")
export plot_frame, plot_frame!, update_visualization!, animate_manipulator!, plot_torques, plot_states, show_frequencies

# creates two manipulator types from stored data
include("./config.jl")
export build_manipulator
   
end # module
