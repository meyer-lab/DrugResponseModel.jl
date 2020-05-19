@testset "Fit All drug at once tests" begin
    concs, populations, g1s, g2s = load(189, 1)
    p = ones(45)

    effects = DrugResponseModel.getODEparamsAll(p, concs)
    #@time optimize_hillAll(concs, g1s, g2s; maxstep = 1E2)
    #@assert(all(x -> !isnan(x), effects))

    # TODO: Profile residHillAll

    # test the local optimization function
    initial = [52.1243, 1.2856, 0.0026682, 2.9096, 0.0250893, 0.0623134, 
        154.803, 1.10582, 0.304235, 0.295368, 0.411795, 0.613183,
        82.5301, 1.05063, 0.391899, 0.919925, 0.0730874, 0.338454,
        3.21729, 3.41872, 0.0841223, 1.52284, 0.281485, 0.207797, 
        60.2809, 0.913883, 0.00228251, 2.96638, 0.0583304, 0.0770743, 
        1.46021, 2.18259, 0.547192, 36.894, 70.8683, 17.2029, 10.6778];
    params = optim_all(concs, g1s, g2s, initial)
    @assert(all(x -> !isnan(x), params))
end
