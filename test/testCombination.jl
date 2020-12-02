@testset "Combination tests where g1 prog. rate increased for one drug and decreases for the other drug and g2 prog. rates where in both drugs the rate is decreasing compared to control." begin
    p1 = ones(9, 8)
    p2 = ones(9, 8)

    # drugA
    p1[1, 1] = 0.6          # g1 prog. rate in control
    p1[1, 2:end] .= 0.8     # g1 prog. rate in all other concentrations
    p1[2, 1] = 1.0          # g2 prog. rate in control
    p1[2, 2:end] .= 0.8     # g2 prog. rate in all other concentrations
    # drugB
    p2[2, 1] = 1.0          # g2 prog. rate in control
    p2[2, 2:end] .= 0.5     # g2 prog. rate in all other concentrations
    p2[1, 1] = 0.6          # g1 prog. rate in control
    p2[1, 2:end] .= 0.3     # g1 prog. rate in all other concentrations
    combination = DrugResponseModel.AllBliss_params(p1, p2)
    @assert(all(combination[1:2, 2:end, 2] .>= 0.3))
    @assert(all(combination[1:2, 2:end, 2] .<= 0.5))
end

@testset "Combination tests from estimated parameters to converting to ODE parameters where for both drugs, rates are decreasing and one reaches to zero." begin
    concs, _, _, _ = load(189, 1)
    gem_before = [10.0, 0.9, 0.9, 1.8, 0.2, 0.5, 0.00593379, 0.110279, 0.5, 10.0, 10.0, 10.0, 10.0]
    dox_before = [100.0, 0.1, 0.04, 0.8, 0.16, 0.0, 0.0720467, 0.14468, 0.5, 10.0, 10.0, 10.0, 10.0]
    p1 = DrugResponseModel.getODEparams(gem_before, concs[:, 3])
    p2 = DrugResponseModel.getODEparams(dox_before, concs[:, 2])

    cmb = DrugResponseModel.Bliss_params_unit(p1[:, 2], p2[:, 3], hcat(p1[:, 1], p2[:, 1]))
    @assert(all(cmb[1, :, end] .>= 0.4))
    @assert(all(cmb[1, :, end] .<= 0.5))
end
