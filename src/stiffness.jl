"""
Computes the actuation space stiffness for one actuator based on gear wheel radi and belt stiffness
    r0:     radius of actuated wheel
    ks:     vector of belt stiffness values from root to end
    ρ:      vector of transmission ratios r_out/r_in from root to end
"""
function actuation_space_stiffness(r0::Real, ks::Vector{<:Real}, ρ::Vector{<:Real})
    kj_ = 0
    
    for i=1:length(ks)
        prod = 1
        if i >= 2
            for l=1:i-1
                prod = prod*ρ[i-1]^2
            end
        end
        kj_ += 1/ks[i]*prod
    end

    return r0^2/kj_
end