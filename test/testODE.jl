@testset "ODE tests" begin
    _, g2, g1 = setup_data("Lapatinib1")

    t = LinRange(0.0, 95.5, size(g1, 1))
    p2 = [100.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.65];
    DrugResponseModel.predict(p, g1[1, 1] + g2[1, 1], t, g1[:, 1], g2[:, 1])
    
    @time DrugResponseModel.predict(p, g1[1, 1] + g2[1, 1], t, g1[:, 1], g2[:, 1])
end
