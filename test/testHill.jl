@testset "Hill tests" begin
    conc, pop, g2, g1 = setup_data("lapatinib")
    params = [10.0, 1.0, 1.0, 1.0, 1.0, 1.0, 10.0, 1.0, 1.0, 30, 30, 10, 10]

    DrugResponseModel.residHill(params, conc, g1, g2)
    @time DrugResponseModel.residHill(params, conc, g1, g2)
    @profile DrugResponseModel.residHill(params, conc, g1, g2)

    Profile.print(noisefloor = 2.0)

    # Check that these at least run
    effects = getODEparams(params, conc)
    num = DrugResponseModel.numcells(effects[:, 3], g1[1] + g2[1], 100)
end
