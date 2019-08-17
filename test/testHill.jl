using Test
using Profile
using CSV
using DrugResponseModel


# import drug concentrations
param_lap_dde = CSV.read(joinpath("..", "data", "params_lap_DDE.csv"))
concentrations = permutedims(Vector(param_lap_dde[8,2:end]));

# make sure there are 8 concentrations
@test length(concentrations) == 8

# import G1, G2, and population data
pop, g2, g1, g2_0, g1_0 = get_data(joinpath("..", "data", "lap.csv"),
                                   joinpath("..", "data", "lap_pop.csv"));

# removing peaks from the data
for i in 1:length(concentrations)
    pop[:, i] = remove_peaks(pop[:, i])
    g2[:, i] = remove_peaks(g2[:, i])
    g1[:, i] = remove_peaks(g1[:, i])
end

# testing the BlackBoxOptim 

# lower bound
low = [50.0, 0.001, 0.2, 0.01, 0.01, 0.001, 0.01, 20.0, 0.1, 0.01, 19.0, 0.01, 0.01, 0.001, 0.01, 0.001, 0.01]
# upper bound
high = [250.0, 0.01, 0.4, 0.04, 0.05, 0.02, 0.04, 35.0, 3.0, 0.04, 23.0, 2.0, 0.4, 0.01, 0.04, 0.04, 0.04]
# guess
guess = [50.0,0.01,0.2,0.01,0.021,0.02,0.0101,20.0286,0.1100,0.0125,21.2519,0.024374,0.0825,0.01,0.01857,0.04,0.0102]
# max number of iterations 
num_steps=10
# do the optimization
parameterrs = optimize_hill(guess, concentrations, g1, g2, g1_0, g2_0, low, high, num_steps)
# check all the parameters to be positive
@test all(x->x>=0.00000001, parameterrs)
@test length(parameterrs) == 17

# profiling for Hill model
@profile optimize_hill(guess, concentrations, g1, g2, g1_0, g2_0, low, high, num_steps)
Profile.print()
