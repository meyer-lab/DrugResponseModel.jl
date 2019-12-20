@testset "Hill tests" begin
	conc, pop, g2, g1, g1_0, g2_0 = setup_data("lapatinib")

    # run the optimization
    fitness, params = optimize_hill(50.0, 350.0, conc_l, g1, g2, g1_0, g2_0)
    println(fitness)

	@time optimize_hill(50.0, 350.0, conc, g1, g2, g1_0, g2_0)

	@profile optimize_hill(50.0, 350.0, conc, g1, g2, g1_0, g2_0)

	Profile.print(noisefloor=5.0)

    # Check that these at least run
    effects = getODEparams(params, conc)
	plot_parameters(conc, effects)
end