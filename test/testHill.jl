@testset "Hill tests" begin
	conc, pop, g2, g1, g1_0, g2_0 = setup_data("lapatinib")
	params = [10.0, 1.0, 1.0, 1.0, 1.0, 1.0, 10.0, 1.0, 1.0, 1.0, 1.0, 30, 30]

	DrugResponseModel.residHill(params, conc, g1, g2, g1_0, g2_0)
	@time DrugResponseModel.residHill(params, conc, g1, g2, g1_0, g2_0)
	@profile DrugResponseModel.residHill(params, conc, g1, g2, g1_0, g2_0)

	Profile.print(noisefloor=5.0)

    # Check that these at least run
    effects = getODEparams(params, conc)
end