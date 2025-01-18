# Copyright (c) 2024 Quan-feng WU <wuquanfeng@ihep.ac.cn>
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

set_random_seed_from_time() = (Random.seed! ∘ Int ∘ mod)(
    parse(BigInt,
        (bytes2hex ∘ sha256 ∘ string ∘ now)(),
        base=16
    ),
    typemax(Int)
)
