@testset "Hill tests" begin

    conc, pop, g2, g1 = setup_data("lapatinib")
    params = [10.0, 1.0, 1.0, 1.0, 1.0, 1.0, 10.0, 1.0, 1.0, 30, 30, 10, 10]

    DrugResponseModel.residHill(params, conc, g1, g2)
    @time DrugResponseModel.residHill(params, conc, g1, g2)
    @profile DrugResponseModel.residHill(params, conc, g1, g2)

    Profile.print(noisefloor = 5.0)

    outt = DrugResponseModel.residHillG(params, conc, g1, g2)

    # Check that these at least run
    effects = getODEparams(params, conc)
end
