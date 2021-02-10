@testset "ODE tests" begin
    _, g2, g1 = setup_data("Lapatinib1")

    t = LinRange(0.0, 95.5, size(g1, 1))
    p2 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.65]
    pC = [1.0, 1.0, 0.4, 1.0, 1.0, 0.3, 1.0, 1.0, 0.65]
    predict(p2, pC, t, g1[:, 1], g2[:, 1])

    @time predict(p2, pC, t, g1[:, 1], g2[:, 1])
end
