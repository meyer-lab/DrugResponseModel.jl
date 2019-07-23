# FactCheck and Base.Test are two functions to use for testing
using Test

include("importData.jl")
include("DDEmodel.jl")


##------------------ Simple tests for get_input function -----------------------##
pop, g2, g1, g1_0, g2_0 = get_data("..//data//lap.csv", "..//data//lap_pop.csv")

@test size(pop[:, 1],1) == size(g2[:, 1],1)
@test size(pop[:, 1],1) == size(g1[:, 1],1)
@test size(g1_0, 1) == size(g2_0, 1)


##------------------ Simple tests for DDEsolve function -------------------------##
# i is the number of the column we are using from the data (# of trial)
i = 6

# initial guess
p  = [0.02798, 0.025502, 21.3481, 10.2881, 0.0001, 0.0001]

# setting lowest delay for tau1 to be half an hour and for tau2 to be 3 hours.
low = [0.015, 0.003, 3.0, 7.0, 0.0001, 0.0001]
upp = [0.075, 0.075, 30.0, 100.0, 0.05, 0.05]

# Estimating the parameters for trial i
params = optimIt(p, low, upp, i, g1, g2)

# to test the estimated parameters are still in the range
@test length(params) == 6
for i in 1:6
    @test upp[i] >= params[i] >= low[i]
end

plotIt(params, g1, g2, g1_0, g2_0, pop, i, "Lapatinib")
