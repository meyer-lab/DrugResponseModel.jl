@testset "Hill tests" begin
    conc, g2, g1 = load(189, 1)
    params = [70.0, 1.0, 0.6, 0.03, 1.4, 2.0, 2.0, 0.2, 1.05, 2.0, 0.003, 1.0e-9, 0.06, 0.1, 0.53]

    # Check that these at least run
    effects = getODEparams(params, conc)
    num = DrugResponseModel.numcells(effects[:, 3, 1], g1[1] + g2[1])
end
