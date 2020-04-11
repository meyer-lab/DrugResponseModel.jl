@testset "Fit All drug at once tests" begin
    concs, populations, g1s, g2s = load(192)
    p = ones(45)

    effects = DrugResponseModel.getODEparamsAll(p, concs)
    @time optimize_hillAll(concs, g1s, g2s; maxstep = 1E2)
    @assert(all(x -> !isnan(x), effects))
    @profile DrugResponseModel.cost(
        effects[:, 1, 1],
        g1s[:, 1, 1],
        g2s[:, 1, 1],
    )
end
