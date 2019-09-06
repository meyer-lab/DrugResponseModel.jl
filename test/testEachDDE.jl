using Test
using Profile
using CSV
using DrugResponseModel



#------------- DDE optimization and checking the error
_, g1_l, g2_l, g1_0_l, g2_0_l = setup_data("lapatinib")
_, g1_d, g2_d, g1_0_d, g2_0_d = setup_data("doxorubicin")
_, g1_g, g2_g, g1_0_g, g2_0_g = setup_data("gemcitabine")
_, g1_t, g2_t, g1_0_t, g2_0_t = setup_data("paclitaxel")

#  [EC50, b_steepness, alpha_min, alpha_max, beta_min, beta_max, tau1_mean, tau1_max, tau2_min, tau2_max, gamma1_max, gamma2_max]
# initial guess
initial_guess  = [0.02798, 0.025502, 21.3481, 10.2881, 0.0001, 0.0001]

# bounds 
lower_bnd = [-8.0, -8.0, 0.0, 0.0, -10.0, -10.0]
upper_bnd = [0.0, 0.0, 6.0, 6.0, 0.0, 0.0]
bound = collect(zip(lower_bnd, upper_bnd))

# max number of steps
maxSteps = 1e4
parameters = zeros(4, 8, 6)

println("############## start test optimize for each dde ###############")
for j in 1:8
    # lapatinib parameters
    println(" ##############lapatinib trial number $j ")
    fitness_l, parameters[1, j, :] = optimization(g1_l, g2_l, g1_0_l, g2_0_l, initial_guess, j, lower_bnd, upper_bnd, maxSteps)
    @test fitness_l <= 8000
    @test all(p -> p>0.0, exp.(parameters[1, :, j]))

    # doxorubicin parameters
    println("############## doxorubicin trial number $j ")
    fitness_d, parameters[2, j, :] = optimization(g1_d, g2_d, g1_0_d, g2_0_d, initial_guess, j, lower_bnd, upper_bnd, maxSteps)
    @test fitness_d <= 8000
    @test all(p -> p>0.0, exp.(parameters[2, :, j]))

    # gemcitabine parameters
    println("############## gemcitabine trial number $j ")
    fitness_g, parameters[3, j, :] = optimization(g1_g, g2_g, g1_0_g, g2_0_g, initial_guess, j, lower_bnd, upper_bnd, maxSteps)
    @test fitness_g <= 8000
    @test all(p -> p>0.0, exp.(parameters[3, :, j]))

    # paclitaxel parameters
    println("############# paclitaxel trial number $j ")
    fitness_t, parameters[4, j, :] = optimization(g1_t, g2_t, g1_0_t, g2_0_t, initial_guess, j, lower_bnd, upper_bnd, maxSteps)
    @test fitness_t <= 8000
    @test all(p -> p>0.0, exp.(parameters[4, :, j]))
end