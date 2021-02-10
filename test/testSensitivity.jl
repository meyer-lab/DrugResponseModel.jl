@testset "Sensitivity tests" begin
    concs, _, _, _ = load(189, 1)
    p = ones(59)
    effs = getODEparams(p, concs)
    # Check that these at least run
    der = DrugResponseModel.get_derivative(p, 1, 2, concs, effs[:, 1, 1], 1)
    @assert(all(der .>= -1000.0))
    @assert(all(der .<= 1000.0))

end
