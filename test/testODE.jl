@testset "ODE tests" begin
	_, pop, g2, g1, g1_0, g2_0 = setup_data("lapatinib")
	# initial values
	p = [8.870525324, 8.492087169, 0.43447323, 7.67847790]

	params_ode = zeros(4,8)
	# ODE optimization and estimation of the parameters
	for i in 1:8
		_, params_ode[:, i] = ODEoptimizer(p, i, g1, g2, g1_0, g2_0, 1, 1)
		# to test the estimated parameters are still in the range
		@test all(upp .>= params_ode[:, i] .>= low)
	end

	# Check that these at least run
	ode_plotIt(params_ode[:, 8], g1, g2, g1_0, g2_0, pop, 8, "", false)
	ODEplot_all(params_ode, g1, g2, g1_0, g2_0, pop)
end