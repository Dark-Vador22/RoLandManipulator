"""
Overload for addition - on Lie algebra
"""
function +(a::so3, b::so3)
    return so3(a.mat + b.mat)
end

function +(a::se3, b::se3)
    return se3(a.mat + b.mat)
end
################################################################################

"""
Overload for subtraction - on Lie algebra
"""
function -(a::so3, b::so3)
    return so3(a.mat - b.mat)
end

function -(a::se3, b::se3)
    return se3(a.mat - b.mat)
end
################################################################################

"""
Overload for scalar and matrix multiplication on Lie algebra and Lie group resp.
"""
function *(a::Number, b::so3)
    return so3(a * b.mat)
end

function *(a::so3, b::Number)
    return so3(a.mat * b)
end

function *(a::Number, b::se3)
    return se3(a * b.mat)
end

function *(a::se3, b::Number)
    return se3(a.mat * b)
end

function *(a::so3, b::so3)
    return so3(a.mat * b.mat)
end

function *(a::se3, b::se3)
    return se3(a.mat * b.mat)
end

function *(a::Number, b::SO3)
    return SO3(a * b.mat)
end

function *(a::SO3, b::Number)
    return SO3(a.mat * b)
end

function *(a::Number, b::SE3)
    return SE3(a * b.mat)
end

function *(a::SE3, b::Number)
    return SE3(a.mat * b)
end

function *(a::SO3, b::SO3)
    return SO3(a.mat * b.mat)
end

function *(a::SO3, b::Array{T}) where T <: Any
    return a.mat * b
end

function *(a::Array{T}, b::SO3) where T <: Any
    return a * b.mat
end

function *(a::SE3, b::SE3)
    return SE3(a.mat * b.mat)
end
################################################################################

"""
Overload for division
"""
function /(a::so3, b::Number)
    return so3(a.mat/b)
end

function /(a::se3, b::Number)
    return se3(a.mat/b)
end

function /(a::SO3, b::Number)
    return SO3(a.mat/b)
end

function /(a::SE3, b::Number)
    return SE3(a.mat/b)
end
################################################################################

"""
Overload for determinant
"""
function det(M::T) where T <: Union{SO3, SE3}
    return det(M.mat)
end
################################################################################

"""
Overload for inverse - analytic formulation
"""
function inv(R::SO3)
    return SO3(transpose(R.mat))
end

function inv(C::SE3)
    mat = [transpose(C.R.mat) -transpose(C.R.mat)*C.r; [0 0 0 1]]
    return SE3(mat)
end
################################################################################

"""
Lie bracket operator for two elements of the Lie algebra
"""
function lie_bracket(a_t::so3, b_t::so3)
    mat = a_t.mat*b_t.mat - b_t.mat*a_t.mat
    return so3(mat)
end

function lie_bracket(a_t::se3, b_t::se3)
    mat = ad(a_t)*b_t.mat
    return se3(mat)
end
################################################################################

"""
Exponential map, analytically defined for so(3) and se(3)
"""
function exp_map(ω_t::so3)
    if iszero(ω_t.ω)
        mat = [1 0 0; 0 1 0; 0 0 1]
    else
        ω = norm(ω_t.ω)
        mat = I + sin(ω)/ω*ω_t.mat + (1-cos(ω))/ω^2*ω_t.mat^2
    end
    return SO3(mat)
end

function exp_map(v_t::se3)
    ω_t = so3(v_t.ω)
    u = v_t.u
    mat = [exp_map(ω_t).mat dexp(ω_t)*u; zeros(1,3) 1]
    return SE3(mat)
end
################################################################################

"""
differential of the exponential map for spatial coordinates
"""
function dexp(ω_t::so3)
    if iszero(ω_t.ω)
        return [1 0 0; 0 1 0; 0 0 1]
    else
        ω = norm(ω_t.ω)
        return I + (1 - cos(ω))/ω^2*ω_t.mat + (ω - sin(ω))/ω^3*ω_t.mat^2
    end
end

function dexp(v_t::se3)
    return [dexp(so3(v_t.ω))  zeros(3,3); mixed_dexp(v_t)  dexp(so3(v_t.ω))]
end
################################################################################

"""
differential of the exponential map for body attached coordinates
"""
function dexp_(x_t::T) where T <: Union{so3, se3}
    return dexp(x_t*(-1))
end
################################################################################

"""
mixed matrix of dexp(se(3))
GEOMETRIC METHODS AND FORMULATIONS IN COMPUTATIONAL MULTIBODY SYSTEM DYNAMICS -
Andreas Müller and Zdravko Terze, 2016
"""
function mixed_dexp(v_t::se3)
    ω = norm(v_t.ω)
    if iszero(v_t.u)
        return zeros(3,3)
    elseif iszero(ω)
        return so3(-v_t.u).mat
    else
        ω_t = so3(v_t.ω)
        u_t = so3(v_t.u)
        a = 2/ω*sin(ω/2)*cos(ω/2)
        b = 4sin(ω/2)^2/ω^2
        h = dot(v_t.ω, v_t.u)/ω^2
        return b/2*u_t.mat + (1-a)/ω^2*lie_cracket(u_t, ω_t).mat + h*(a-b)/ω*ω_t.mat + h/ω^2*(b/2 - 3(1-a)/ω)*ω_t.mat^2
    end
end
################################################################################

"""
Closed form solutions for the logarithmic map
"""
function log_map(R::SO3)
    if near_zero(R.mat - I)
        θ = 0
        mat = zeros(3,)
    elseif near_zero(tr(R.mat) + 1)
        θ = π
        if ! near_zero(1 + R.mat[3,3])
            mat = 1/√(2(1 + R.mat[3,3]))*[R.mat[1,3], R.mat[2,3], 1 + R.mat[3,3]]
        elseif ! near_zero(1 + R.mat[2,2])
            mat = 1/√(2(1 + R.mat[2,2]))*[R.mat[1,2], 1 + R.mat[2,2], R.mat[3,2]]
        elseif ! near_zero(1 + R.mat[1,1])
            mat = 1/√(2(1 + R.mat[1,1]))*[1 + R.mat[1,1], R.mat[2,1], R.mat[3,1]]
        end
    else
        acosinput = (tr(R.mat) - 1)/2
        if acosinput > 1
            acosinput = 1
        elseif acosinput < -1
            acosinput = -1
        end
        θ = acos(acosinput)
        mat = θ/(2*sin(θ))*(R.mat - transpose(R.mat))
    end
    return so3(mat), θ
end

function log_map(C::SE3)
    if near_zero(C.R.mat - I) && ! near_zero(C.r)
        θ = norm(C.r)
        u = C.r/θ
        mat = [zeros(3,3) u; zeros(1,3) 0]
    else
        ω_t, θ = log_map(C.R)
        G_ = transpose(inv(dexp_(ω_t)))
        # G_ = 1/θ*I - 0.5ω_t.mat + (1/θ - 0.5cot(θ/2))*ω_t.mat^2
        u = G_*C.r
        mat = [ω_t.mat u; zeros(1,3) 0]
    end
    return se3(mat), θ
end
################################################################################

"""
Vectorial map of the special Euclidean group - linearized form of the
logarithmic map
"""
function vect(C::SE3)
    ϕ_t = (C.R.mat - transpose(C.R.mat))/2
    return vcat([ϕ_t[3,2], ϕ_t[1,3], ϕ_t[2,1]], C.r)
end

function vecto(C::SE3)
    ϕ_t, θ = log_map(C)
    return ϕ_t.vec
end
################################################################################

"""
Adjoint representation for SE(3) in linear space - generally, one has:
        AdR(̃a) = RãRᵀ, AdR(a) = Ra
        AdH(̃a) = HãH⁻¹, AdH(a) = Xa
where this function returns X
"""
function Ad(C::SE3)
    mat = [C.R.mat zeros(3,3); so3(C.r).mat*C.R.mat C.R.mat]
    return mat
end
################################################################################

"""
Adjoint for Lie algebra (hat-operator)
"""
function ad(v_t::se3)
    mat = [so3(v_t.ω).mat zeros(3,3); so3(v_t.u).mat so3(v_t.ω).mat]
    return mat
end
################################################################################

################################################################################
############################ HELPER FUNCTIONS ##################################

"""
Function that computes the product of two Lie algebras by summing
"""
function lie_cracket(a_t::so3, b_t::so3)
    return a_t*b_t + b_t*a_t
end
################################################################################

"""
Determines values close to zero
"""
function near_zero(val::Real)
    return abs(val) < 1e-8
end

function near_zero(mat::AbstractArray)
    return all(near_zero.(mat))
end
