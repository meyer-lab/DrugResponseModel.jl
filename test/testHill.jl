
println("#####################   hill tests begin ...")
# import G1, G2, population, and concentrations data
conc_l, popl, g2l, g1l, g2_0l, g1_0l = setup_data("lapatinib")

# make sure there are 8 concentrations
@test length(conc_l) == 8

# testing the BlackBoxOptim 

# guess
guess = [125.0, 0.04, 0.007, 0.005, 0.007, 0.005, 39.0, 30.0, 11.0, 20.0, 0.003, 0.02]
# max number of iterations 
num_steps=500

# profiling for Hill model
residHill(guess, conc_l, g1l, g2l, g1_0l, g2_0l)
Profile.clear()
println("profiling for residHill function  \n")
@profile residHill(guess, conc_l, g1l, g2l, g1_0l, g2_0l)
Profile.print(noisefloor=1.0)

# lower bound
low = [50.0, 0.01, 0.005, 0.04, 0.005, 0.01, 22.0, 20.0, 6.0, 5.0, 0.0001, 0.0001]
# upper bound
high = [250, 0.01, 0.02, 0.1, 0.2, 0.03, 40.0, 35.0, 15.0, 20.0, 0.05, 0.1]

# do the optimization
println("### lapatinib ###")
best_fitL, pt_l = optimize_hill(guess, conc_l, g1l, g2l, g1_0l, g2_0l, num_steps)

# check all the parameters to be positive
@test all(x -> x>0, pt_l)
@test length(pt_l) == 12

# make sure all the dde parameters are positive
dde_paramsL = getDDEparams(pt_l, conc_l)

# test the fitness of the model
@test best_fitL <= 5e5
