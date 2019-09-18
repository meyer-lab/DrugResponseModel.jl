using Test
using Profile
using CSV
using DrugResponseModel
using BlackBoxOptim

println("#####################   hill tests begin ...")
# import G1, G2, population, and concentrations data
conc_l, popl, g2l, g1l, g2_0l, g1_0l = setup_data("lapatinib")
conc_d, popd, g2d, g1d, g2_0d, g1_0d = setup_data("doxorubicin")
conc_g, popg, g2g, g1g, g2_0g, g1_0g = setup_data("gemcitabine")
conc_t, popt, g2t, g1t, g2_0t, g1_0t = setup_data("paclitaxel")
# make sure there are 8 concentrations
@test length(conc_l) == length(conc_d) == length(conc_g) == length(conc_t) == 8

# testing the BlackBoxOptim 

# guess
l_guess =  [125.0, 0.04, 0.007, 0.005, 0.007, 0.005, 30.0, 20.0, 0.003, 0.02]
d_guess =  [125.0, 0.04, 0.007, 0.005, 0.007, 0.005, 30.0, 20.0, 0.003, 0.02]
g_guess =  [125.0, 0.04, 0.007, 0.005, 0.007, 0.005, 30.0, 20.0, 0.003, 0.02]
t_guess =  [125.0, 0.04, 0.007, 0.005, 0.007, 0.005, 30.0, 20.0, 0.003, 0.02]
# max number of iterations 
num_steps=5

# profiling for Hill model
println("profiling for residHill function  \n")
@profile residHill(l_guess, conc_l, g1l, g2l, g1_0l, g2_0l)
Profile.print(noisefloor=10.0)

residue(hillParams) = residHill(l_guess, conc_l, g1l, g2l, g1_0l, g2_0l)

# lower bound
low = [50.0, 0.01, 0.005, 0.04, 0.005, 0.01, 20.0, 5.0, 0.0001, 0.0001]
# upper bound
high = [250, 2.0, 0.02, 0.1, 0.2, 0.03, 40.0, 20.0, 0.05, 0.1]

println("profiling bboptimize in Hill ")
@profile bboptimize(residue; SearchRange=collect(zip(low, high)), MaxSteps=100, Method =:adaptive_de_rand_1_bin_radiuslimited)
Profile.print(noisefloor=10.0)

@profile optimize_hill(l_guess, conc_l, g1l, g2l, g1_0l, g2_0l, num_steps)
Profile.print(noisefloor=10.0)

# do the optimization
println("### lapatinib ###")
best_fitL, pt_l = optimize_hill(l_guess, conc_l, g1l, g2l, g1_0l, g2_0l, num_steps)
println("### doxorubicin ###")
best_fitD, pt_d = optimize_hill(d_guess, conc_d, g1d, g2d, g1_0d, g2_0d, num_steps)
println("### gemcitabine ###")
best_fitG, pt_g = optimize_hill(g_guess, conc_g, g1g, g2g, g1_0g, g2_0g, num_steps)
println("### paclitaxel ###")
best_fitT, pt_t = optimize_hill(t_guess, conc_t, g1t, g2t, g1_0t, g2_0t, num_steps)

Profile.init(delay=0.5)
# profiling the get_dde params
println("profiling the getDDEparams  \n")
@profile getDDEparams(pt_l, conc_l)
# check all the parameters to be positive
@test all(x -> x>0, pt_l)
@test length(pt_l) == 10
@test all(x -> x>0, pt_d)
@test length(pt_d) == 10
@test all(x -> x>0, pt_g)
@test length(pt_g) == 10
@test all(x -> x>0, pt_t)
@test length(pt_t) == 10

# make sure all the dde parameters are positive
dde_paramsL = getDDEparams(pt_l, conc_l)
dde_paramsD = getDDEparams(pt_d, conc_d)
dde_paramsG = getDDEparams(pt_g, conc_g)
dde_paramsT = getDDEparams(pt_t, conc_t)

# test the fitness of the model
@test best_fitL <= 1e6
@test best_fitD <= 1e6
@test best_fitG <= 1e6
@test best_fitT <= 1e6

