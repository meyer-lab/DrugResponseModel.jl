@testset "Hill tests" begin
    conc, g2, g1 = load(189, 1)
    params = [70.0, 1.0, 0.6, 0.03, 1.4, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 0.2, 1.05, 2.0, 0.003, 1.0e-9, 0.06, 0.01, 0.1, 0.1]

    # Check that these at least run
    effects = getODEparams(params, conc)
end
