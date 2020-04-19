@testset "Fit All drug at once tests" begin
    concs, populations, g1s, g2s = load(189, 1)
    p = ones(45)

    effects = DrugResponseModel.getODEparamsAll(p, concs)
    @time optimize_hillAll(concs, g1s, g2s; maxstep = 1E2)
    @assert(all(x -> !isnan(x), effects))

    # TODO: Profile residHillAll

    # test the local optimization function
    params = optim_all(concs, g1s, g2s)
    @assert(all(x -> !isnan(x), params))
    @profile optim_all(concs, g1s, g2s)

    Profile.print(noisefloor = 2.0)
end
