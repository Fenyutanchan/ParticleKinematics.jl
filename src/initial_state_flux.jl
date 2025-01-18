# Copyright (c) 2024 Quan-feng WU <wuquanfeng@ihep.ac.cn>
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

function initial_state_flux(p_list::Vector{<:Vector})
    @assert all(==(4) ∘ length, p_list) "Four-momenta are required, but got $p_list."

    n = length(p_list)
    
    n == 0 && return 0
    n == 1 && return 2 * (mass ∘ first)(p_list)
    if n == 2
        s = (scalar_product ∘ sum)(p_list)
        ma_sqr, mb_sqr = map(scalar_product, p_list)

        return 2 * (sqrt ∘ λ)(s, ma_sqr, mb_sqr)
    end

    @error "We do not implement the case of n ≥ 2."
end
