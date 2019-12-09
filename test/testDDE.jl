@testset "DDE tests" begin
	##------------------ Import data -----------------------##
	conc, pop, g2, g1, g1_0, g2_0 = setup_data("lapatinib")

	@test size(pop[:, 1],1) == size(g2[:, 1],1)
	@test size(pop[:, 1],1) == size(g1[:, 1],1)
	@test size(g1_0, 1) == size(g2_0, 1)

	for i in 1:8
		pop[:, i] = remove_peaks(pop[:, i])
		g2[:, i] = remove_peaks(g2[:, i])
		g1[:, i] = remove_peaks(g1[:, i])
	end
	@test size(g1[:, 1],1) == size(g2[:, 1],1)

	# initial guess
	initial_guess  = [0.02798, 0.025502, 15.3481, 15.2881, 0.001, 0.001]

	# bounds 
	lower_bnd = [-6.0, -6.0, 1.0, 1.0, -10.0, -10.0]
	upper_bnd = [0.0, 0.0, 6.0, 6.0, 0.0, 0.0]

	# max number of steps
	maxSteps = 10000
	times = range(0.0; stop = 95.5, length = 192)

	# Add profile of fitness function?

	parameters = zeros(6, 8)
	# Estimating the parameters for all trials
	for j in 1:8
		best_fit, parameters[:, j] = optimization(g1, g2, g1_0, g2_0, initial_guess, j, lower_bnd, upper_bnd, maxSteps)
		println("trial number $j")
		# to test the estimated parameters are still in the range
		@test all(exp.(upper_bnd) .>= parameters[:, j] .>= exp.(lower_bnd))
		@test best_fit <= 9000
	end

	# profiling the plot function
	plotIt(parameters[:, 8], 6, "", :false, pop, g2, g1, g2_0, g1_0)
	plot_all(parameters, pop, g2, g1, g2_0, g1_0)
	plot_parameters(conc, parameters)
end