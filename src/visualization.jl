"""
Depiction of robot URDF by MeshCat visualizer
"""
function update_visualization!(man::Manipulator)
    set_configuration!(man.mvis, man.θ[2:end])
end

"""
Shows torques from collected control array 
"""
function plot_torques(times::AbstractVector, τ::AbstractVector)
    plt = plot(xlabel = "time [s]", ylabel = "torque [Nm]")
    for i = 1:length(τ[1])
        plot!(plt, times[1:end-1], [τ[j][i] for j = 1:length(τ)], lw = 2, ls = :dash, alpha = 0.5, label = "joint $i")
    end
    return plt
end

function plot_torques(times::AbstractVector, τ::AbstractVector, ζ::AbstractVector)
    plt = plot_torques(times, τ)
    for i = 1:length(ζ[1])
        plot!(plt, times[1:end-1], [ζ[j][i] for j = 1:length(ζ)], lw = 2, label = "actuator $i")
    end
    return plt
end
        

"""
Shows position and velocity of collected state vector
"""
function plot_states(times::AbstractVector, xs::AbstractVector; size = (600,400))
    plt = plot(xlabel = "time [s]", ylabel = "x, v [rad, rad/sec]", size = size)
    pal = palette(:lightrainbow)
    n = Int(length(xs[1])/2)
    for i = 1:n
        plot!(plt, times, [xs[j][i] for j = 1:length(xs)], lw = 2, c = pal[i], label = "pos $i")
        plot!(plt, times, [xs[j][n+i] for j = 1:length(xs)], lw = 2, ls = :dash, c = pal[i], label = "vel $i")
    end
    return plt
end

"""
Simply plots a frame defined by SE3
"""
function plot_frame!(plt::Plots.Plot, fr::SE3; scale::Real = 0.1, lw::Int = 2, style::Symbol = :solid)
    for i = 1:3
        if i == 1
            color = :red
        elseif i == 2
            color = :green
        elseif i == 3
            color = :blue
        end
       
        plot!(plt, [fr.r[1], fr.r[1] + fr.R.mat[1,i]*scale], 
                   [fr.r[2], fr.r[2] + fr.R.mat[2,i]*scale], 
                   [fr.r[3], fr.r[3] + fr.R.mat[3,i]*scale], c = color, lw = lw, ls = style)
    end
     # lets all axes appear equally long
     enforce_isometry!(plt)
    return plt
end

function plot_frame(fr::SE3; scale::Real = 0.1, lw::Int = 2, style::Symbol = :solid)
    plt = plot3d(xlabel="x", ylabel="y", zlabel="z", camera = (30,30), legend = false, aspec_ratio = :equal)
    plot_frame!(plt, fr, scale = scale, lw = lw, style = style)
     # lets all axes appear equally long
     enforce_isometry!(plt)
    return plt
end

"""
Enforces equal scaling of the axis on 3D-plots
"""
function enforce_isometry!(plt::Plots.Plot)
    x_limits = xlims(plt)
    y_limits = ylims(plt)
    z_limits = zlims(plt)

    x_range = abs(x_limits[2] - x_limits[1])
    x_middle = sum(x_limits)/2
    y_range = abs(y_limits[2] - y_limits[1])
    y_middle = sum(y_limits)/2
    z_range = abs(z_limits[2] - z_limits[1])
    z_middle = sum(z_limits)/2

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max(x_range, y_range, z_range)

    X = [x_middle - plot_radius, x_middle + plot_radius]
    Y = [y_middle - plot_radius, y_middle + plot_radius]
    Z = [z_middle - plot_radius, z_middle + plot_radius]

    plot!(plt, X, Y, Z, linealpha = 0)

    return plt
end

"""
Saves screenshots of specified trajectory, adjust properties beforehand in the browser
"""
function vis_frame_saver(man::Manipulator, X::AbstractArray, indices)
    if maximum(indices) > length(X)
        throw(error("indices are out of trajectory bounds"))
    end 
    for i in indices
        set_configuration!(man.mvis, X[i][1:4])
        MeshCat.save_image(man.vis)
    end 
end

"""
Scales transparency by chaning alpha channel 
"""
scaletransparency(image, scale) = RGBA.(RGB.(image), scale.*alpha.(image))

"""
Create overlay plot with alpha gradient 
"""
function frame_overlay(path, n_files; flipdir=false, alpha_start=0, alpha_end=1)
    alpha = range(alpha_start,alpha_end,length=n_files)
    if flipdir
        alpha = reverse(alpha)
    end     
    img = FileIO.load(path*"meshcat.png")    
    img=scaletransparency(img, alpha[1])
    plt =  plot(img, size=(1200,800),xticks=false,yticks=false,axis=([], false))   
    for i in eachindex(alpha[2:end])
        img=FileIO.load(path*"meshcat($i).png")
        img=scaletransparency(img, alpha[i])
        plot!(img)        
    end 
    return plt
end 




