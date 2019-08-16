using Test
using Profile
using DrugResponseModel


# import drug concentrations
param_lap_dde = CSV.read(".//figures//Dox//params_dox_DDE.csv")
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
@testset "Hill model test." begin
    # lower bound
    low = [50.0, 0.001, 0.2, 0.01, 0.01, 0.001, 0.01, 20.0, 0.1, 0.01, 19.0, 0.01, 0.01, 0.001, 0.01, 0.001, 0.01]
    # upper bound
    high = [250.0, 0.01, 0.4, 0.04, 0.05, 0.02, 0.04, 35.0, 3.0, 0.04, 23.0, 2.0, 0.4, 0.01, 0.04, 0.04, 0.04]
    params = optimize_hill(low, high)
    # check all the parameters to be positive
    @test all(params .>= 0.0)
    @test length(params == 17)
end


# profiling
low = [50.0, 0.001, 0.2, 0.01, 0.01, 0.001, 0.01, 20.0, 0.1, 0.01, 19.0, 0.01, 0.01, 0.001, 0.01, 0.001, 0.01]
# upper bound
high = [250.0, 0.01, 0.4, 0.04, 0.05, 0.02, 0.04, 35.0, 3.0, 0.04, 23.0, 2.0, 0.4, 0.01, 0.04, 0.04, 0.04]

# profiling for Hill model
@profile optimize_hill(low, high)
Profile.print()
