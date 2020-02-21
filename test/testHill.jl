@testset "Hill tests" begin
    conc, pop, g2, g1 = setup_data("lapatinib")
#     params = [10.0, 1.0, 1.0, 1.0, 1.0, 1.0, 10.0, 1.0, 1.0, 30, 30, 10, 10]
    params = [60.0257, 0.925401, 6.87587e-9, 0.792481, 0.368774, 0.750335, 0.00827934, 0.0204481, 0.536002, 20.7355, 13.4557, 0.530622, 2.83636]

    DrugResponseModel.residHill(params, conc, g1, g2)
    @time DrugResponseModel.residHill(params, conc, g1, g2)
    @profile DrugResponseModel.residHill(params, conc, g1, g2)

    Profile.print(noisefloor = 5.0)

    outt = DrugResponseModel.residHillG(params, conc, g1, g2)

    # Check cell number is not negative
    effects = getODEparams(params, conc)
    for n=1:7
        num = numcells(effects[:, n], g1[1]+g2[1], 100)
    end
    dif = zeros(4, length(conc) - 1)
    # assert all the derivations of death rate are negative.
    for i = 2:length(conc)
        dif[:, i] = diffCell(effects[:, i], g1[1]+g2[1], 100)[1:4]
        @assert(all(x -> x<=0, dif[3:4, i]), "positive gradient for death rates due to parameter set: $params")
    end
end
