@testset "ODE tests" begin
    _, pop, g2, g1 = setup_data("lapatinib")

    t = range(0.0; stop = 95.5, length = 192)
    p = [1.0, 1.1, 0.2, 0.3, 0.5, 20, 20, 10, 10]
    DrugResponseModel.predict(p, 1.0, t, 20, 20, 10, 10)

    # ODE optimization and estimation of the parameters
    fitness, params_ode = ODEoptimizer(4, g1, g2)
    @test fitness < 800.0

    @time DrugResponseModel.cost(p, g1[:, 1], g2[:, 1], 100, 100, 20, 20)
end
