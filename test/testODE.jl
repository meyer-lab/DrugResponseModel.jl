@testset "ODE tests" begin
	_, pop, g2, g1, g1_0, g2_0 = setup_data("lapatinib")

	t = range(0.0; stop = 95.5, length = 192)
	p = [1.0, 1.1, 0.2, 0.3, 0.5, 20, 20]
	DrugResponseModel.predict(p, 1.0, 1.0, t, 20, 20)

	# ODE optimization and estimation of the parameters
	fitness, params_ode = ODEoptimizer(4, g1, g2, g1_0, g2_0)
	println(fitness)

	@time DrugResponseModel.cost(p, g1_0[1], g2_0[1], g1[:, 1], g2[:, 1], 100, 100)

	for ii in 1:10
		@profile DrugResponseModel.cost(p, g1_0[1], g2_0[1], g1[:, 1], g2[:, 1], 100, 100)
	end

	Profile.print(noisefloor=2.0)
end