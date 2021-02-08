@testset "Sensitivity tests" begin
    concs, _, _, _ = load(189, 1)
    p = ones(59)
    g0 = 20.0
    effs = getODEparams(p, concs, 5)
    # Check that these at least run
    der = DrugResponseModel.get_derivative(p, 1, 2, concs, g0, 1)
    @assert(all(der .>= -1000.0))
    @assert(all(der .<= 1000.0))

end
