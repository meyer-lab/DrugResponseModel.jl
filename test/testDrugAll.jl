@testset "Fit All drug at once tests" begin
    concs, populations, g1s, g2s = load(189, 1)
    p = ones(41)
    effects = DrugResponseModel.getODEparamsAll(p, concs)

    # TODO: Profile residHillAll

    # test the local optimization function
    params = DrugResponseModel.optim_all(concs, g1s, g2s)
end
