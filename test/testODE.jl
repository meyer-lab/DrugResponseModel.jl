using Test
using Profile
using DrugResponseModel

println("####################  ODE model tests begin ... ")

_, pop, g2, g1, g1_0, g2_0 = setup_data("lapatinib")
# initial values
p = [8.870525324, 8.492087169, 0.43447323, 7.67847790]

# setting lowest delay for tau1 to be half an hour and for tau2 to be 3 hours.
low = 0.0001*ones(4)
upp = 0.1*ones(4)

# ODE optimization and estimation of the parameters
for i in 1:8
    params_ode = ODEoptimizer(low, upp, p, i, g1, g2, g1_0, g2_0)
    # to test the estimated parameters are still in the range
    @test length(params_ode) ==4
    for i in 1:4
        @test upp[i] >= params_ode[i] >= low[i]
    end
end
params_ode = ODEoptimizer(low, upp, p, 5, g1, g2, g1_0, g2_0)
@profile ODEoptimizer(low, upp, p, 5, g1, g2, g1_0, g2_0)
@profile ode_plotIt(params_ode, g1, g2, g1_0, g2_0, pop, 8, "")

