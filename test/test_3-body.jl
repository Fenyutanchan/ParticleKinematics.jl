# Copyright (c) 2024 Quan-feng WU <wuquanfeng@ihep.ac.cn>
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

function MC_3_body(num_points::Int)
    MeV = 1
    mu = 931.49410242 * MeV
    M = 8.0246073 * mu
    Mp = 8.00530510 * mu

    P = [M, 0, 0, 0]
    mass_list = [Mp, 0, 0]

    weight_sum = zeros(Real, Threads.nthreads())
    prefactor = 1 / (2 * M)

    p = Progress(num_points; desc="Evaluating...")
    counter = Threads.Atomic{Int}(0)
    ProgressMeter.update!(p, counter[])
    Threads.@threads for _ ∈ 1:num_points
        ρn, pi_list  = phase_space_point_Monte_Carlo(P, mass_list)

        weight_sum[Threads.threadid()] += ρn

        Threads.atomic_add!(counter, 1)
        ProgressMeter.update!(p, counter[])
    end

    return sum(weight_sum) * prefactor / num_points
end
