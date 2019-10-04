# FactCheck and Base.Test are two functions to use for testing
using Test
using Profile
using DrugResponseModel

println("####################  DDE model tests begin ... ")
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

##------------------ Simple tests for DDEsolve function -------------------------##

# initial guess
initial_guess  = [0.02798, 0.025502, 15.3481, 15.2881, 0.001, 0.001]

# bounds 
lower_bnd = [-6.0, -6.0, 1.0, 1.0, -10.0, -10.0]
upper_bnd = [0.0, 0.0, 6.0, 6.0, 0.0, 0.0]

# max number of steps
maxSteps = 10000
times = range(0.0; stop = 95.5, length = 192)

# profiling to DDEmodel
println("  \n profiling find_history  \n ")
@profile find_history(g1, g2)
Profile.print(noisefloor=10.0)
println("  \n profiling ddesolve function  \n")
@profile ddesolve(times, g1, g2, g1_0, g2_0, initial_guess, 6)
Profile.print(noisefloor=10.0)
println("  \n profiling optimization function   \n")
@profile optimization(g1, g2, g1_0, g2_0, initial_guess, 6, lower_bnd, upper_bnd, maxSteps)
Profile.print(noisefloor=10.0)

parameters = zeros(6, 8)
println("+++++++++++++++++ trials for lapatinib +++++++++++++++++++")
# Estimating the parameters for all trials
for j in 1:8
    best_fit, parameters[:, j] = optimization(g1, g2, g1_0, g2_0, initial_guess, j, lower_bnd, upper_bnd, maxSteps)
    println("trial number $j")
    # to test the estimated parameters are still in the range
    for i in 1:6
        @test exp.(upper_bnd[i]) >= parameters[i,j] >= exp.(lower_bnd[i])
        @test best_fit <= 9000
    end

end


best_fit, parameter = optimization(g1, g2, g1_0, g2_0, initial_guess, 6, lower_bnd, upper_bnd, maxSteps)
# profiling the plot function
plotIt(parameter, 6, "", :false, pop, g2, g1, g2_0, g1_0)
plot_all(parameters, pop, g2, g1, g2_0, g1_0)
plot_parameters(conc, parameters)
