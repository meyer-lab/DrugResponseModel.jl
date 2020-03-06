@testset "Combination tests" begin
    concs, populations, g1s, g2s = load()
    g0 = g1s[1,1,1]+g2s[1,1,1];
    p = [50.2552, 0.866725, 2.41775e-5, 1.96785, 0.0185627, 0.0447356, 108.176, 1.16593, 1.3052, 2.86005e-5, 0.348876, 0.438606, 14.4648, 2.15448, 1.42977, 2.99996, 0.0216477, 0.196135, 4.03273, 3.93278, 1.86825e-5, 0.867674, 0.430368, 0.241053, 1.68177, 1.17228, 0.508321, 34.5786, 37.9849, 20.7341, 9.11408];
    effects = DrugResponseModel.getODEparamsAll(p, concs)
    # combination test case
    p1 = effects[:,:,1]
    p1[2,1] = 1.0
    p1[2,2:end] .= 0.8 
    p2 = effects[:,:,2]
    p2[2,1] = 1.0
    p2[2,2:end] .= 0.5
    combination = BlissCombination(p1, p2, 8)
    @assert(all(combination[2,2:end,2] .>= 0.3))
    @assert(all(combination[2,2:end,2] .<= 0.5))
end