@testset "Hill tests" begin
	conc, pop, g2, g1, g1_0, g2_0 = setup_data("lapatinib")

	t = range(0.0; stop = 95.5, length = 192)
    guess = [100.0, 1.0, 1.0, 1.0, 100, 1.0, 1.0, 1.0, 100, 1.0, 1.0, 100.0, 1.0, 1.0, 0.5]
    # run the optimization
    fitness, params = optimize_hill(guess, conc, g1, g2, g1_0, g2_0, 50.0, 300.0, 1, 1)
    println(fitness)

	@time optimize_hill(guess, conc, g1, g2, g1_0, g2_0, 50.0, 300.0, 1, 1)

	@profile optimize_hill(guess, conc, g1, g2, g1_0, g2_0, 50.0, 300.0, 1, 1)

	Profile.print(noisefloor=5.0)

    # Check that these at least run
    effects = getODEparams(params, conc)
	plot_parameters(conc, effects)
end