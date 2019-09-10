using Test
using Profile
using CSV
using DrugResponseModel

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
l_guess = [100.0, 1.0 , 0.007, 0.01, 0.05, 0.045, 8.0286, 30.1100, 12.0, 8.0, 0.0035, 0.04]
d_guess = [150.8, 0.09, 0.06, 0.009, 0.008, 0.043, 38.94, 5.49, 19.9, 3.2, 0.0029, 0.022]
g_guess = [80.0, 0.097836, 0.00682578, 0.00958013, 0.00503788, 0.0433964, 38.941, 5.49712, 9.94827, 7.20807, 0.00294624, 0.0222821]
t_guess = [80.0, 0.097836, 0.00682578, 0.00958013, 0.00503788, 0.0433964, 38.941, 5.49712, 9.94827, 7.20807, 0.00294624, 0.0222821]
# max number of iterations 
num_steps=500

# do the optimization
println("### lapatinib ###")
best_fitL, pt_l = optimize_hill(l_guess, conc_l, g1l, g2l, g1_0l, g2_0l, num_steps)
println("### doxorubicin ###")
best_fitD, pt_d = optimize_hill(d_guess, conc_d, g1d, g2d, g1_0d, g2_0d, num_steps)
println("### gemcitabine ###")
best_fitG, pt_g = optimize_hill(g_guess, conc_g, g1g, g2g, g1_0g, g2_0g, num_steps)
println("### paclitaxel ###")
best_fitT, pt_t = optimize_hill(t_guess, conc_t, g1t, g2t, g1_0t, g2_0t, num_steps)

# check all the parameters to be positive
@test all(x -> x>0, pt_l)
@test length(pt_l) == 12
@test all(x -> x>0, pt_d)
@test length(pt_d) == 12
@test all(x -> x>0, pt_g)
@test length(pt_g) == 12
@test all(x -> x>0, pt_t)
@test length(pt_t) == 12

# make sure all the dde parameters are positive
dde_paramsL = getDDEparams(pt_l, conc_l)
dde_paramsD = getDDEparams(pt_d, conc_d)
dde_paramsG = getDDEparams(pt_g, conc_g)
dde_paramsT = getDDEparams(pt_t, conc_t)

for i in 1:6
    temp = dde_paramsL[i]
    @test all(x -> x>0, temp)
    temp = dde_paramsD[i]
    @test all(x -> x>0, temp)
    temp = dde_paramsG[i]
    @test all(x -> x>0, temp)
    temp = dde_paramsT[i]
    @test all(x -> x>0, temp)
end

# test the fitness of the model
@test best_fitL <= 9e5
@test best_fitD <= 9e5
@test best_fitG <= 9e5
@test best_fitT <= 9e5
# profiling for Hill model
@profile optimize_hill(guess, conc_l, g1l, g2l, g1_0l, g2_0l, num_steps)
Profile.print(noisefloor=10.0)
