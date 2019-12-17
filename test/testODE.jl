@testset "ODE tests" begin
	_, pop, g2, g1, g1_0, g2_0 = setup_data("lapatinib")

	# ODE optimization and estimation of the parameters
	fitness, params_ode = ODEoptimizer(ones(4)*0.5, 4, g1, g2, g1_0, g2_0, 1, 1)

	@time ODEoptimizer(ones(4)*0.5, 4, g1, g2, g1_0, g2_0, 1, 1)

	@profile ODEoptimizer(ones(4)*0.5, 4, g1, g2, g1_0, g2_0, 1, 1)

	Profile.print(noisefloor=5.0)

	# Check that these at least run
	ode_plotIt(params_ode, g1, g2, g1_0, g2_0, pop, 4, "", false, 1, 1)
end