@testset "ODE tests" begin
    _, pop, g2, g1, g1_0, g2_0 = setup_data("lapatinib")

    t = range(0.0; stop = 95.5, length = 192)
    p = [1.0, 1.1, 0.2, 0.3, 0.5, 20, 20, 10, 10]
    DrugResponseModel.predict(p, 1.0, 1.0, t, 20, 20, 10, 10)

    # ODE optimization and estimation of the parameters
    fitness, params_ode = ODEoptimizer(4, g1, g2, g1_0, g2_0)
    @test fitness < 800.0
end
