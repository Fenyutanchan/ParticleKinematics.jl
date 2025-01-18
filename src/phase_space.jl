# Copyright (c) 2024 Quan-feng WU <wuquanfeng@ihep.ac.cn>
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

function phase_space_point_Monte_Carlo(P, mass_list)
    @assert length(P) == 4 "Four-momentum are required, but got $P."

    n = length(mass_list)

    # set_random_seed_from_time()

    iszero(n) && return 0, []
    n == 1 && return 1, [P]

    ρn = 1

    M_list = []
    push!(M_list, mass(P))
    for ii ∈ n-1:-1:2
        μi = sum(mass_list[1:ii])
        integration_range = first(M_list) - mass_list[ii+1] - μi
        pushfirst!(M_list, rand() * integration_range + μi)
        ρn *= integration_range
    end
    pushfirst!(M_list, first(mass_list))

    ρ_list = [(sqrt ∘ λ)(M_list[ii]^2, M_list[ii-1]^2, mass_list[ii]^2) / (2 * M_list[ii]) for ii ∈ 2:n]
    ρn *= prod(ρ_list) / (2^n * M_list[end] * (2 * pi)^(3 * n - 4))

    pi_list = []
    ki_list = []
    push!(ki_list, P)
    θ_list = acos.(rand(n) * 2 .- 1)
    ϕ_list = rand(n) * 2 * pi
    ρn *= (4 * pi)^(n - 1)

    for ii ∈ n:-1:2
        Ei = (M_list[ii]^2 + mass_list[ii]^2 - M_list[ii-1]^2) / (2 * M_list[ii])
        pp = sqrt(Ei^2 - mass_list[ii]^2)
        Ek = (M_list[ii]^2 + M_list[ii-1]^2 - mass_list[ii]^2) / (2 * M_list[ii])

        rest_pi = [
            Ei,
            pp * sin(θ_list[ii]) * cos(ϕ_list[ii]),
            pp * sin(θ_list[ii]) * sin(ϕ_list[ii]),
            pp * cos(θ_list[ii])
        ]
        pushfirst!(pi_list, Lorentz_boost(rest_pi, first(ki_list)))

        rest_ki = [Ek, -rest_pi[2:4]...]
        pushfirst!(ki_list, Lorentz_boost(rest_ki, first(ki_list)))
    end
    pushfirst!(pi_list, first(ki_list))

    return ρn, pi_list
end

"""
    phase_space_point(P, mass_list, t_list)

Calculate the phase space point according to the four-momentum `P`, the mass list `mass_list` and the parameter list `t_list`.
Notice that the number of parameters should be ``3 n - 4``, where ``n`` is the number of particles in the final state.
For all `t` in `t_list`, `0 < t < 1` is required.
"""
function phase_space_point(P, mass_list, t_list; check_t_list=true)
    @assert length(P) == 4 "Four-momentum are required, but got $P."
    n = length(mass_list)

    iszero(n) && return 0, []
    n == 1 && return 1, [P]

    n_parameter = length(t_list)
    t_list_index = 1
    @assert n_parameter == 3 * n - 4 "The number of parameters should be $(3 * n - 4), but got $n_parameter."
    check_t_list && @assert all(0 .<= t_list .<= 1) "All `t` in `t_list` should be `0 <= t <= 1`, but got $t_list."
    # t_list_copy = copy(t_list)

    ρ = 1

    M_list = []
    push!(M_list, mass(P))
    for ii ∈ n-1:-1:2
        μi = sum(mass_list[1:ii])
        integration_range = first(M_list) - mass_list[ii+1] - μi
        pushfirst!(M_list, t_list[t_list_index] * integration_range + μi)
        t_list_index += 1
        ρ *= integration_range
    end
    pushfirst!(M_list, first(mass_list))

    ρ_list = [(sqrt ∘ λ)(M_list[ii]^2, M_list[ii-1]^2, mass_list[ii]^2) / (2 * M_list[ii]) for ii ∈ 2:n]
    ρ *= prod(ρ_list) / (2^n * M_list[end] * (2 * π)^(3 * n - 4))

    pi_list = []
    ki_list = []
    push!(ki_list, P)
    ρ *= (4 * pi)^(n - 1)
    for ii ∈ n:-1:2
        Ei = (M_list[ii]^2 + mass_list[ii]^2 - M_list[ii-1]^2) / (2 * M_list[ii])
        pp = sqrt(Ei^2 - mass_list[ii]^2)
        Ek = (M_list[ii]^2 + M_list[ii-1]^2 - mass_list[ii]^2) / (2 * M_list[ii])

        rest_pi = [Ei, (pp * __generate_n̂(t_list[t_list_index], t_list[t_list_index+1]))...]
        t_list_index += 2
        pushfirst!(pi_list, Lorentz_boost(rest_pi, first(ki_list)))

        rest_ki = [Ek, -rest_pi[2:4]...]
        pushfirst!(ki_list, Lorentz_boost(rest_ki, first(ki_list)))
    end
    pushfirst!(pi_list, first(ki_list))

    @assert t_list_index == n_parameter + 1

    return ρ, pi_list
end
function phase_space_point(P, mass_list)
    n = length(mass_list)
    t_list = rand(3 * n - 4)

    return phase_space_point(P, mass_list, t_list)
end
