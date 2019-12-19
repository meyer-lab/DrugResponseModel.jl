@testset "ODE tests" begin
	_, pop, g2, g1, g1_0, g2_0 = setup_data("lapatinib")

	t = range(0.0; stop = 95.5, length = 192)
	DrugResponseModel.predict([1.0, 1.1, 0.2, 0.3, 0.5, 20, 20], 1.0, 1.0, t)

	# ODE optimization and estimation of the parameters
	fitness, params_ode = ODEoptimizer(4, g1, g2, g1_0, g2_0)
	println(fitness)

	@time ODEoptimizer(4, g1, g2, g1_0, g2_0)

	@profile ODEoptimizer(4, g1, g2, g1_0, g2_0)

	Profile.print(noisefloor=5.0)

	# Check that these at least run
	ode_plotIt(params_ode, g1, g2, g1_0, g2_0, pop, 4, "", false)
end