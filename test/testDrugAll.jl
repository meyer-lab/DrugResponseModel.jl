@testset "Fit All drug at once tests" begin
    concs, populations, g1s, g2s = load()
    p = [0.5*maximum(concs[:, 1]), 1.0, 1.0, 1.0, 1.0, 1.0, 0.1, 0.1, 0.50, 0.5*maximum(concs[:, 2]), 1.0, 1.0, 1.0, 1.0, 1.0, 0.1, 0.1, 0.50, 0.5*maximum(concs[:, 3]), 1.0, 1.0, 1.0, 1.0, 1.0, 0.1, 0.1, 0.50, 0.5*maximum(concs[:, 4]), 1.0, 1.0, 1.0, 1.0, 1.0, 0.1, 0.1, 0.50, 6, 10, 50, 50]
    effects = DrugResponseModel.getODEparamsAll(p, concs)
    @time optimize_hillAll(concs, g1s, g2s; maxstep = 1E2)
    @assert(all(x -> !isnan(x), effects))
    
end