# Copyright (c) 2024 Quan-feng WU <wuquanfeng@ihep.ac.cn>
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

function __generate_n̂(t1, t2)
    cosθ = 2 * t1 - 1
    sinθ = sqrt(1 - cosθ^2)
    sinϕ, cosϕ = sincospi(2 * t2)

    return [sinθ * cosϕ, sinθ * sinϕ, cosθ]
end
__generate_n̂() = __generate_n̂(rand(), rand())

@doc raw"""
    generate_value_from_semi_infinite_range(v_end, t=rand(); p=1, direction=:+)

Generate ``v`` in the range of ``[v_\mathrm{end}, \infty)`` (for `direction=:+`) or ``(-\infty, v_\mathrm{end}]`` (for `direction=:-`) according to
``
    v = v_\mathrm{end} \pm \left( \frac{t}{1 - t} \right)^p
``
where ``t \in [0, 1)`` and $p > 0$.
"""
function generate_value_from_semi_infinite_range(v_end, t=rand(); p=1, direction=:+)
    @assert 0 ≤ t < 1 "The parameter `t` should be in the range of [0, 1)."
    @assert p > 0 "The parameter `p` should be positive."

    if direction ∉ [:+, :-]
        @warn """
        The direction should be either `:+` or `:-`, but got `$(direction)`.
        Default to `:+`.
        """
        direction = :+
    end

    v = v_end + (direction == :+ ? 1 : -1) * (t / (1 - t))^p
    ρ = p * t^(p - 1) / (1 - t)^(p + 1)

    return ρ, v
end

function generate_value_from_log10_range(log10_v_min, log10_v_max, t=rand())
    @assert log10_v_min <= log10_v_max "The minimum value should be less than or equal to the maximum value."

    log_10 = (log ∘ convert)(typeof(t), 10)

    v = exp10(t * (log10_v_max - log10_v_min) + log10_v_min)
    ρ = (log10_v_max - log10_v_min) * log_10 * v

    return ρ, v
end
