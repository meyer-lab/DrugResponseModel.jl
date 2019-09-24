# FactCheck and Base.Test are two functions to use for testing
using Test
using Profile
using DrugResponseModel

println("####################  DDE model tests begin ... ")
##------------------ Import data -----------------------##
_, pop, g2, g1, g1_0, g2_0 = setup_data("lapatinib")
_, popd, g2d, g1d, g1_0d, g2_0d = setup_data("doxorubicin")
_, popg, g2g, g1g, g1_0g, g2_0g = setup_data("gemcitabine")
_, popt, g2t, g1t, g1_0t, g2_0t = setup_data("paclitaxel")

@test size(pop[:, 1],1) == size(g2[:, 1],1)
@test size(pop[:, 1],1) == size(g1[:, 1],1)
@test size(g1_0, 1) == size(g2_0, 1)

for i in 1:8
    pop[:, i] = remove_peaks(pop[:, i])
    g2[:, i] = remove_peaks(g2[:, i])
    g1[:, i] = remove_peaks(g1[:, i])
end
@test size(g1[:, 1],1) == size(g2[:, 1],1)

##------------------ Simple tests for DDEsolve function -------------------------##

# initial guess
initial_guess  = [0.02798, 0.025502, 15.3481, 15.2881, 0.001, 0.001]

# bounds 
lower_bnd = [-6.0, -6.0, 1.0, 1.0, -10.0, -10.0]
upper_bnd = [0.0, 0.0, 6.0, 6.0, 0.0, 0.0]

# max number of steps
maxSteps = 1
j = 6
times = range(0.0; stop = 95.5, length = 192)

# profiling to DDEmodel
println("  \n profiling find_history  \n ")
@profile find_history(g1, g2)
Profile.print(noisefloor=10.0)
println("  \n profiling ddesolve function  \n")
@profile ddesolve(times, g1, g2, g1_0, g2_0, initial_guess, j)
Profile.print(noisefloor=10.0)
println("  \n profiling optimization function   \n")
@profile optimization(g1, g2, g1_0, g2_0, initial_guess, j, lower_bnd, upper_bnd, maxSteps)
Profile.print(noisefloor=10.0)

println("+++++++++++++++++ trials for lapatinib +++++++++++++++++++")
# Estimating the parameters for all trials
for j in 1:8
    best_fit, parameters = optimization(g1, g2, g1_0, g2_0, initial_guess, j, lower_bnd, upper_bnd, maxSteps)
    println("trial number $j")
    # to test the estimated parameters are still in the range
    @test length(parameters) == 6
    for i in 1:6
        @test exp.(upper_bnd[i]) >= parameters[i] >= exp.(lower_bnd[i])
        @test best_fit <= 9000
    end

end

println("+++++++++++++++++ trials for doxorubicin +++++++++++++++++++")
for j in 1:8
    best_fit, parameters = optimization(g1d, g2d, g1_0d, g2_0d, initial_guess, j, lower_bnd, upper_bnd, maxSteps)
    println("trial number $j")
    # to test the estimated parameters are still in the range
    @test length(parameters) == 6
    for i in 1:6
        @test exp.(upper_bnd[i]) >= parameters[i] >= exp.(lower_bnd[i])
        @test best_fit <= 9000
    end

end

println("+++++++++++++++++ trials for gemcitabine +++++++++++++++++++")
for j in 1:8
    best_fit, parameters = optimization(g1g, g2g, g1_0g, g2_0g, initial_guess, j, lower_bnd, upper_bnd, maxSteps)
    println("trial number $j")
    # to test the estimated parameters are still in the range
    @test length(parameters) == 6
    for i in 1:6
        @test exp.(upper_bnd[i]) >= parameters[i] >= exp.(lower_bnd[i])
        @test best_fit <= 9000
    end

end

println("+++++++++++++++++ trials for paclitaxel +++++++++++++++++++")
for j in 1:8
    best_fit, parameters = optimization(g1t, g2t, g1_0t, g2_0t, initial_guess, j, lower_bnd, upper_bnd, maxSteps)
    println("trial number $j")
    # to test the estimated parameters are still in the range
    @test length(parameters) == 6
    for i in 1:6
        @test exp.(upper_bnd[i]) >= parameters[i] >= exp.(lower_bnd[i])
        @test best_fit <= 9000
    end

end

best_fit, parameters = optimization(g1t, g2t, g1_0t, g2_0t, initial_guess, 6, lower_bnd, upper_bnd, maxSteps)
# profiling the plot function
@profile plotIt(parameters, 6, "", :false)
