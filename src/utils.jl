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

"""
Returns moving average of array vs. Mean value is calucated using n points  
"""
moving_average(vs,n) = [sum(@view vs[i:(i+n-1)])/n for i in 1:(length(vs)-(n-1))]

"""
see if matrix is getting rank deficient 
"""
isinvertible(A::Matrix{Float64}) = !isapprox(det(BigFloat.(A)), 0, atol = 1e-10)

