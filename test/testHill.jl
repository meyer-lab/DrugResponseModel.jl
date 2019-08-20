using Test
using Profile
using CSV
using DrugResponseModel


# import drug concentrations
param_lap_dde = CSV.read(joinpath("..", "data", "params_lap_DDE.csv"))
concentrations = permutedims(Vector(param_lap_dde[8,2:end]));

# make sure there are 8 concentrations
@assert length(concentrations) == 8

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

# guess
guess = [100.0, 1.0 , 0.007, 0.01, 0.05, 0.045, 8.0286, 30.1100, 12.0, 8.0, 0.0035, 0.04]
# max number of iterations 
num_steps=10
# do the optimization
_, pt = optimize_hill(guess, concentrations, g1, g2, g1_0, g2_0, num_steps)

# check all the parameters to be positive
@assert all(x -> x>0, pt)
@assert length(pt) == 12

# profiling for Hill model
@profile optimize_hill(guess, concentrations, g1, g2, g1_0, g2_0, num_steps)
Profile.print()
