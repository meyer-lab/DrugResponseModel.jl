@testset "ODE tests" begin
	_, pop, g2, g1, g1_0, g2_0 = setup_data("lapatinib")
	# initial values
	p = [8.870525324, 8.492087169, 0.43447323, 7.67847790]

	low = 0.0001*ones(4)
	upp = 0.1*ones(4)
	params_ode = zeros(4,8)
	# ODE optimization and estimation of the parameters
	for i in 1:8
		params_ode[:,i] = ODEoptimizer(low, upp, p, i, g1, g2, g1_0, g2_0)
		# to test the estimated parameters are still in the range
		for j in 1:4
			@test upp[j] >= params_ode[j,i] >= low[j]
		end
	end

	# Check that these at least run
	params_odes = ODEoptimizer(low, upp, p, 5, g1, g2, g1_0, g2_0)
	ODEoptimizer(low, upp, p, 5, g1, g2, g1_0, g2_0)
	ode_plotIt(params_ode[:, 8], g1, g2, g1_0, g2_0, pop, 8, "", false)
	ODEplot_all(params_ode, g1, g2, g1_0, g2_0, pop)
    turingODE(5, params_ode, g1, g2, g1_0, g2_0)
end