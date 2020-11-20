@testset "Sensitivity tests" begin
    concs, _, _, _ = load(189, 1);
    p = [97.7434, 1.41726, 0.0501761, 0.930205, 0.111494, 0.123604, 0.552216, 
        92.9045, 1.40024, 0.996396, 0.0513087, 0.521831, 0.6535, 0.566578, 
        15.5317, 2.3689, 0.592433, 0.999986, 0.0283363, 0.286975, 0.503328, 
        3.96929, 4.62768, 0.0512281, 0.307528, 0.549714, 0.378717, 0.50959, 
        63.4248, 0.976052, 0.16582, 0.740009, 0.0572609, 0.0776912, 0.534201, 
        0.734513, 0.375555, 16.8387, 12.3945, 30.176, 14.5352]
    g0 = 20.0
    effs = getODEparamsAll(p, concs);
    combin = DrugResponseModel.calc_cellNum(effs[:, 2, 1], effs[:, 3, 2], g0)
    # Check that these at least run
    der = DrugResponseModel.get_derivative(effs[:, 2, 1], effs[:, 3, 2], effs, combin, concs, 1, 2, g0)

end