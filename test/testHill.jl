using Test
using Profile
using CSV
using DrugResponseModel
using BlackBoxOptim

println("#####################   hill tests begin ...")
# import G1, G2, population, and concentrations data
conc_l, popl, g2l, g1l, g2_0l, g1_0l = setup_data("lapatinib")

# make sure there are 8 concentrations
@test length(conc_l) == length(conc_d) == length(conc_g) == length(conc_t) == 8

# testing the BlackBoxOptim 

# guess
guess =  [125.0, 0.04, 0.007, 0.005, 0.007, 0.005, 30.0, 20.0, 0.003, 0.02]
# max number of iterations 
num_steps=50

# profiling for Hill model
println("profiling for residHill function  \n")
@profile residHill(guess, conc_l, g1l, g2l, g1_0l, g2_0l)
Profile.print(noisefloor=10.0)

residue(hillParams) = residHill(guess, conc_l, g1l, g2l, g1_0l, g2_0l)

# lower bound
low = [50.0, 0.01, 0.005, 0.04, 0.005, 0.01, 20.0, 5.0, 0.0001, 0.0001]
# upper bound
high = [250, 2.0, 0.02, 0.1, 0.2, 0.03, 40.0, 20.0, 0.05, 0.1]

println("profiling bboptimize in Hill ")
@profile bboptimize(residue; SearchRange=collect(zip(low, high)), MaxSteps=100, Method =:adaptive_de_rand_1_bin_radiuslimited)
Profile.print(noisefloor=10.0)

# do the optimization
println("### lapatinib ###")
best_fitL, pt_l = optimize_hill(guess, conc_l, g1l, g2l, g1_0l, g2_0l, num_steps)

Profile.init(delay=0.5)
# profiling the get_dde params
println("profiling the getDDEparams  \n")
@profile getDDEparams(pt_l, conc_l)
# check all the parameters to be positive
@test all(x -> x>0, pt_l)
@test length(pt_l) == 10

# make sure all the dde parameters are positive
dde_paramsL = getDDEparams(pt_l, conc_l)

# test the fitness of the model
@test best_fitL <= 5e5
