# Copyright (c) 2024 Quan-feng WU <wuquanfeng@ihep.ac.cn>
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

# Metric tensor
MetricTensor() = diagm([1, -1, -1, -1])
MT = MetricTensor

# Källén function
Källén_λ_function(x, y, z) = (x - y - z)^2 - 4 * y * z
λ = Källén_λ_function

# scalar product
scalar_product(p, q) = transpose(p) * MT() * q
scalar_product(p) = scalar_product(p, p)
SP = scalar_product
mass(p) = try
    (sqrt ∘ SP)(p)
catch
    @error "Invalid time-like four-momentum: $p"
end

# Lorentz boost
function get_Lorentz_boost_matrix(p)
    m = mass(p)
    γ = first(p) / m
    pp = sum(p[begin+1:end].^2)

    if iszero(pp)
        @assert iszero(γ - 1)
        return I(4)
    end

    Λ = Matrix{Real}(undef, 4, 4)
    Λ[1, 1] = γ
    for ii ∈ 2:4
        Λ[1, ii] = Λ[ii, 1] = p[ii] * sqrt((γ^2 - 1) / pp)
        for jj ∈ ii:4
            Λ[ii, jj] = Λ[jj, ii] = (ii == jj ? 1 : 0) + (γ - 1) * p[ii] * p[jj] / pp
        end
    end
    return Λ
end

Lorentz_boost(p, p_ref) = get_Lorentz_boost_matrix(p_ref) * p
