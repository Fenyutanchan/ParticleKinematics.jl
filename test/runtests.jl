# Copyright (c) 2024 Quan-feng WU <wuquanfeng@ihep.ac.cn>
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using FytcKinematics
using ProgressMeter
using Test

include("test_2-body.jl")
include("test_3-body.jl")

include("legacy/kinematics.jl")
(include ∘ joinpath)("legacy", "test_2-body.jl")
(include ∘ joinpath)("legacy", "test_3-body.jl")

@testset "Legacy Kinematics Test" begin
    @info "legacy kinematics test"
    for _ ∈ 1:100
        M = rand() * 1000
        m_list = [rand() for _ ∈ 1:rand(2:10)]
        while sum(m_list) > M
            m_list = [rand() for _ ∈ 1:rand(2:10)]
        end

        pM = rand(3) * 1000
        EM = sqrt(transpose(pM) * pM + M^2)
        P = [EM, pM...]

        @info "M: $M, m_list: $m_list"
        ρn, pi_list = Φ(P, m_list)
        @test sum(pi_list) ≈ P
        @test calc_mass.(pi_list) ≈ m_list rtol=1e-6
    end
end

@testset "Kinematics Test" begin
    @info "kinematics test"
    for _ ∈ 1:100
        M = rand() * 1000
        m_list = [rand() for _ ∈ 1:rand(2:10)]
        while sum(m_list) > M
            m_list = [rand() for _ ∈ 1:rand(2:10)]
        end

        pM = rand(3) * 1000
        EM = sqrt(transpose(pM) * pM + M^2)
        P = [EM, pM...]
        
        @info "M: $M, m_list: $m_list"
        ρn, pi_list = phase_space_point_Monte_Carlo(P, m_list)
        # @show P - sum(pi_list)
        # @show FytcKinematics.mass.(pi_list) - m_list
        @test sum(pi_list) ≈ P
        @test FytcKinematics.mass.(pi_list) ≈ m_list rtol=1e-6
    end
end

@testset "Two-Body Phase Space" begin
    num_points = 10^6
    @info "number of points: $num_points"

    target_result = 5 / (4 * pi * 10)

    legacy_result = MC_2_body_legacy(num_points)
    legacy_rel_error = abs((legacy_result - target_result) / target_result)

    result = MC_2_body(num_points)
    rel_error = abs((result - target_result) / target_result)

    @info "target result: $target_result"

    @info "legacy result: $legacy_result, relative error: $legacy_rel_error"
    @test legacy_rel_error < 1e-10

    @info "result: $result, relative error: $rel_error"
    @test rel_error < 1e-10
end

@testset "Three-Body Phase Space" begin
    num_points = 10^7
    @info "number of points: $num_points"

    legacy_result = MC_3_body_legacy(num_points)
    result = MC_3_body(num_points)

    rel_error = abs((result - legacy_result) / legacy_result)

    @info "legacy result: $legacy_result"
    @info "result: $result"
    @info "relative error: $rel_error"
    @test legacy_result ≈ result rtol=1e-2 # The relative error is large
end
