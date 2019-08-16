# FactCheck and Base.Test are two functions to use for testing
using Test
using Profile
using DrugResponseModel


##------------------ Simple tests for get_input function -----------------------##
pop, g2, g1, g1_0, g2_0 = get_data("..//data//lap.csv", "..//data//lap_pop.csv")

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
lower_bnd = [-6.0, -6.0, 2.0, 2.0, -10.0, -10.0]
upper_bnd = [0.0, 0.0, 6.0, 6.0, 0.0, 0.0]
bound = collect(zip(lower_bnd, upper_bnd))

# Estimating the parameters for trial j
params = optimization(g1, g2, g1_0, g2_0, initial_guess, j, bound)

# to test the estimated parameters are still in the range
@test size(params) == 6
for i in 1:6
    @test upper_bnd[i] >= log.params[i] >= lower_bnd[i]
end

# profiling to DDEmodel
@profile optimization(g1, g2, g1_0, g2_0, initial_guess, j, bound)
Profile.print()

