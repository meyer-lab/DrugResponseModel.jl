@testset "Fit All drug at once tests" begin
    concs, populations, g1s, g2s = load(189, 1)
    p = ones(76)
    effects = DrugResponseModel.getODEparams(p, concs)

    # test the local optimization function
    params = DrugResponseModel.optim_all(concs, g1s, g2s, maxiter = 2)
end
