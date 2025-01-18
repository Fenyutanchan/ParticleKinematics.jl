# Copyright (c) 2024 Quan-feng WU <wuquanfeng@ihep.ac.cn>
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

# function generate_single_four_momentum(msqr; debug::Bool=false, log10_psqr_max_ratio=10, log10_psqr_min_ratio=-20)
#     @assert msqr ≥ 0
#     @assert log10_psqr_min_ratio ≤ log10_psqr_max_ratio

#     ρ = 1

#     q̄, t̄, φbar = rand(3)

#     # q = tan(π * q̄ / 2)
#     # ρ *= π / (1 + cos(π * q̄))
#     q = exp10((q̄ * (log10_psqr_max_ratio - log10_psqr_min_ratio) + log10_psqr_min_ratio) / 2)
#     ρ *= (log10_psqr_max_ratio - log10_psqr_min_ratio) * log(10) * q / 2

#     cosθ = 2 * t̄ - 1
#     sinθ = sqrt(1 - cosθ^2)
#     # ρ *= 2

#     φ = 2 * π * φbar
#     # ρ *= 2 * π

#     E = sqrt(msqr + q^2)
#     # ρ *= q^2 / ((2 * π)^3 * 2 * E)
#     ρ *= q^2 / (4 * π^2 * E)

#     mom = [E, q * sinθ * cos(φ), q * sinθ * sin(φ), q * cosθ]
#     if debug && !(scalar_product(mom) ≈ msqr)
#         @warn """
#         Input: $msqr
#         Output: $ρ, $mom
#         Check: scalar_product(mom) = $(scalar_product(mom))
#         """
#     end

#     return ρ, mom
# end

function generate_single_four_momentum(m, t1, t2, t3; log10_p_min_ratio=-10, log10_p_max_ratio=5)
    @assert m ≥ 0
    @assert log10_p_min_ratio ≤ log10_p_max_ratio

    ρ, p = generate_value_from_log10_range(log10_p_min_ratio, log10_p_max_ratio, t1)
    E = sqrt(p^2 + m^2)

    ρ *= p^2 / (4 * π^2 * E)

    return ρ, [E, (p * __generate_n̂(t2, t3))...]
end
generate_single_four_momentum(m; log10_p_min_ratio=-10, log10_p_max_ratio=5) =
    generate_single_four_momentum(m, rand(), rand(), rand();
        log10_p_min_ratio=log10_p_min_ratio, log10_p_max_ratio=log10_p_max_ratio
    )
