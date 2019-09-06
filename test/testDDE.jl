# FactCheck and Base.Test are two functions to use for testing
using Test
using Profile
using DrugResponseModel

##------------------ Import data -----------------------##
_, pop, g2, g1, g1_0, g2_0 = setup_data("doxorubicin")

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
# j is the number of the column we are using from the data (# of trial)
j = 6

# initial guess
initial_guess  = [0.02798, 0.025502, 21.3481, 10.2881, 0.0001, 0.0001]

# bounds 
lower_bnd = [-6.0, -6.0, 0.0, 0.0, -10.0, -10.0]
upper_bnd = [0.0, 0.0, 6.0, 6.0, 0.0, 0.0]

# max number of steps
maxSteps = 5

# Estimating the parameters for all trials
for j in 1:8
    best_fit, parameters = optimization(g1, g2, g1_0, g2_0, initial_guess, j, lower_bnd, upper_bnd, maxSteps)
    println("trial number $j")
    # to test the estimated parameters are still in the range
    @test length(parameters) == 6
    for i in 1:6
        @test upper_bnd[i] >= parameters[i] >= lower_bnd[i]
    end

end

# profiling to DDEmodel
@profile optimization(g1, g2, g1_0, g2_0, initial_guess, j, lower_bnd, upper_bnd, maxSteps)
Profile.print(noisefloor=10.0)
