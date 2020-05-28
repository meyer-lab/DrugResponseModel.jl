@testset "Hill tests" begin
    conc, pop, g2, g1 = setup_data("Lapatinib1")
    params = [10.0, 1.0, 1.0, 1.0, 1.0, 1.0, 10.0, 1.0, 1.0, 30, 30, 10, 10]

    # Check that these at least run
    effects = getODEparams(params, conc)
    num = DrugResponseModel.numcells(effects[:, 3], g1[1] + g2[1])
end

@testset "Test function appropriately works for different array dimensions" begin
    conc, pop, g2, g1 = setup_data("Lapatinib1")
    p1 = [10.0, 1.0, 1.0, 1.0, 1.0, 1.0, 10.0, 1.0, 1.0, 30, 30, 0, 0]
    eff1 = getODEparams(p1, conc)
end
