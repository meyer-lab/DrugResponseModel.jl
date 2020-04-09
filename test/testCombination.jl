@testset "Combination tests" begin
    concs, populations, g1s, g2s = load(192)
    g0 = g1s[1, 1, 1] + g2s[1, 1, 1]
    p = [198.616, 1.04429, 0.0386625, 2.11079, 0.172845, 0.194753, 125.853, 1.28641, 2.67474, 0.0194261, 0.57024, 0.690045, 16.9456, 2.29013, 1.85685, 0.990835, 0.0339845, 0.243558, 4.31167, 4.2881, 0.743701, 0.45077, 0.72208, 0.36916, 85.2162, 0.818204, 0.063744, 0.999612, 0.0541879, 0.0775908, 2.53385, 0.492065, 0.511441, 52.5239, 16.4035, 33.3496, 13.6045]
    effects = DrugResponseModel.getODEparamsAll(p, concs)
    # combination test case
    p1 = effects[:, :, 1]
    p1[2, 1] = 1.0
    p1[2, 2:end] .= 0.8
    p2 = effects[:, :, 2]
    p2[2, 1] = 1.0
    p2[2, 2:end] .= 0.5
    combination = BlissCombination(p1, p2, 8)
    @assert(all(combination[2, 2:end, 2] .>= 0.3))
    @assert(all(combination[2, 2:end, 2] .<= 0.5))
end
