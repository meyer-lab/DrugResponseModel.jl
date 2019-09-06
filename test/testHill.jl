using Test
using Profile
using CSV
using DrugResponseModel

# import G1, G2, population, and concentrations data
concentrations, pop, g2, g1, g2_0, g1_0 = setup_data("lapatinib")

# make sure there are 8 concentrations
@test length(concentrations) == 8

# testing the BlackBoxOptim 

# guess
guess = [100.0, 1.0 , 0.007, 0.01, 0.05, 0.045, 8.0286, 30.1100, 12.0, 8.0, 0.0035, 0.04]
# max number of iterations 
num_steps=100
# do the optimization
_, pt = optimize_hill(guess, concentrations, g1, g2, g1_0, g2_0, num_steps)

# check all the parameters to be positive
@test all(x -> x>0, pt)
@test length(pt) == 12

# make sure all the dde parameters are positive
dde_params = getDDEparams(pt, concentrations)
for i in 1:6
    temp = dde_params[i]
    @test all(x -> x>0, temp)
end

# profiling for Hill model
@profile optimize_hill(guess, concentrations, g1, g2, g1_0, g2_0, num_steps)
Profile.print(noisefloor=10.0)
