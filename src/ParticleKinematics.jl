# Copyright (c) 2024 Quan-feng WU <wuquanfeng@ihep.ac.cn>
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

module ParticleKinematics

using Dates
using LinearAlgebra
using Random
using SHA

export initial_state_flux

export phase_space_point_Monte_Carlo
export phase_space_point

export generate_single_four_momentum

export generate_value_from_semi_infinite_range
export generate_value_from_log10_range

export scalar_product

include("initial_state_flux.jl")
include("kinematics.jl")
include("phase_space.jl")
include("random_seed.jl")
include("single.jl")
include("utilities.jl")

end # module ParticleKinematics
