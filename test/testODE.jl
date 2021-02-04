@testset "ODE tests" begin
    _, g2, g1 = setup_data("Lapatinib1")

    t = LinRange(0.0, 95.5, size(g1, 1))
    p = [1.0, 1.1, 0.2, 0.3, 0.5, 20, 20, 10, 10]
    DrugResponseModel.predict(p, g1[1, 1] + g2[1, 1], t, g1[:, 1], g2[:, 1])

    @time DrugResponseModel.predict(p, g1[1, 1] + g2[1, 1], t, g1[:, 1], g2[:, 1])
end
