"""
Lie algebra of rotation matrices
"""
struct so3
    mat::AbstractArray{T, 2} where T <: Real
    ω::AbstractArray{T, 1} where T <: Real

    function so3(ωx::Real, ωy::Real, ωz::Real)
        ω = [ωx, ωy, ωz]
        mat = [0 -ωz ωy; ωz 0 -ωx; -ωy ωx 0]
        return new(mat, ω)
    end

    function so3(ω::AbstractArray{T, 1} where T <: Real)
        @assert length(ω) == 3
        mat = [0 -ω[3] ω[2]; ω[3] 0 -ω[1]; -ω[2] ω[1] 0]
        return new(mat, ω)
    end

    function so3(ama::AbstractArray{T, 2} where T <: Real)
        @assert size(ama) == (3, 3)
        ω = [ama[3,2], ama[1,3], ama[2,1]]
        mat = ama
        return new(mat, ω)
    end
end
################################################################################

"""
Lie group of rotation matrices
"""
struct SO3
    mat::AbstractArray{T, 2} where T <: Real

    function SO3(ω::AbstractArray{T, 1} where T <: Real)
        @assert length(ω) == 3
        return exp_map(so3(ω))
    end

    function SO3(ω_t::so3)
        return exp_map(ω_t)
    end

    function SO3(mat::AbstractArray{T, 2} where T <: Real)
        @assert size(mat) == (3, 3)
        @assert abs(1 - det(mat)) <= 1e-5
        return new(mat)
    end
end
################################################################################

"""
Lie algebra of spatial transformations - among the matrix representation, there
is also a vector representation for the isomorphism to Rᵏ and seperate values of
u and ω
"""
struct se3
    mat::AbstractArray{T, 2} where T <: Real
    vec::AbstractArray{T, 1} where T <: Real
    ω::AbstractArray{T, 1} where T <: Real
    u::AbstractArray{T, 1} where T <: Real

    function se3(v::AbstractArray{T, 1} where T <: Real)
        @assert length(v) == 6
        ω = v[1:3]
        u = v[4:6]
        mat = [0 -ω[3] ω[2] u[1]; ω[3] 0 -ω[1] u[2]; -ω[2] ω[1] 0 u[3]; 0 0 0 0]
        vec = v
        return new(mat, vec, ω, u)
    end

    function se3(mat::AbstractArray{T, 2} where T <: Real)
        @assert size(mat) == (4, 4)
        ω = [mat[3,2], mat[1,3], mat[2,1]]
        u = mat[1:3, 4]
        vec = vcat(ω, u)
        return new(mat, vec, ω, u)
    end
end
################################################################################

"""
Lie group of se(3) - special Euclidean group
"""
struct SE3
    mat::AbstractArray{T, 2} where T <: Real
    r::AbstractArray{T, 1} where T <: Real
    R::SO3

    function SE3(v::AbstractArray{T, 1} where T <: Real)
        @assert length(v) == 6
        return exp_map(se3(v))
    end

    function SE3(v_t::se3)
        return exp_map(v_t)
    end

    function SE3(mat::AbstractArray{T, 2} where T <: Real)
        @assert size(mat) == (4, 4)
        r = mat[1:3, 4]
        R = SO3(mat[1:3, 1:3])
        return new(mat, r, R)
    end

end
################################################################################
