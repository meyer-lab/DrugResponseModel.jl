""" In this file we fit all the drugs att once. """

function getODEparamsAll(p::Array{Float64, 1}, concentrations::Array{Float64, 2})
    effects = zeros(9, length(concentrations[:, 1]), 5)

    k = 1
    # Scaled drug effect
    for i = 1:5
        xx = 1.0 ./ (1.0 .+ (p[k] ./ concentrations[:, i]) .^ p[k + 1])

        effects[1, :, i] = p[31] .+ (p[k + 2] - p[31]) .* xx
        effects[2, :, i] = p[32] .+ (p[k + 3] - p[32]) .* xx
        effects[3, :, i] = p[k + 4] .* xx
        effects[4, :, i] = p[k + 5] .* xx
        k += 6
    end
    effects[5, :, :] .= p[33] #percentage in G1
    effects[6, :, :] .= p[34] #nG1
    effects[7, :, :] .= p[35] #nG2
    effects[8, :, :] .= p[36] #nD1
    effects[9, :, :] .= p[37] #nD2

    return effects
end

function residHillAll(hillParams::Array{Float64, 1}, concentrations::Array{Float64, 2}, g1::Array{Float64, 3}, g2::Array{Float64, 3})
    res = Atomic{eltype(hillParams)}(0.0)

    # Solve for all drugs
    t = 1
    for j = 1:5
        hill = [hillParams[t], hillParams[31], hillParams[t+2], hillParams[t+1], hillParams[32], hillParams[t+3], hillParams[t+4], hillParams[t+5], hillParams[33], hillParams[34], hillParams[35], hillParams[36], hillParams[37]]
        t += 6
        @threads for ii = 1:length(concentrations[:, j])
            atomic_add!(res, residHill(hill, conc_l, g1, g2))
        end
    end

    return res[]
end

""" Hill optimization function for all drugs. """
function optimize_hillAll(concs::Array{Float64, 2}, g1::Array{Float64, 3}, g2::Array{Float64, 3}; maxstep = 1E5)
    hillCostAll(hillParams) = residHillAll(hillParams, concs, g1, g2)

    # The parameters used here in order:
    #(:Lap_EC50, :Lap_steepness, :Lap_maxG1ProgRate, :Lap_maxG2ProgRate, :Lap_maxDeathG1Rate, :Lap_maxDeathG2Rate, :Dox_EC50, :Dox_steepness, :Dox_maxG1ProgRate, :Dox_maxG2ProgRate, :Dox_maxDeathG1Rate, :Dox_maxDeathG2Rate, :Gem_EC50, :Gem_steepness, :Gem_maxG1ProgRate, :Gem_maxG2ProgRate, :Gem_maxDeathG1Rate, :Gem_maxDeathG2Rate, :Tax_EC50, :Tax_steepness, :Tax_maxG1ProgRate, :Tax_maxG2ProgRate, :Tax_maxDeathG1Rate, :Tax_maxDeathG2Rate, :pal_EC50, :pal_steepness, :pal_maxG1ProgRate, :pal_maxG2ProgRate, :pal_maxDeathG1Rate, :pal_maxDeathG2Rate, :G1ProgRateControl, :G2ProgRateControl, :percG1, :nG1, :nG2, :nD1, :nD2)
    lowPiece = [0.01, 1e-9, 1e-9, 0.0, 0.0]
    low = vcat(
        minimum(concs[:, 1]),
        lowPiece,
        minimum(concs[:, 2]),
        lowPiece,
        minimum(concs[:, 3]),
        lowPiece,
        minimum(concs[:, 4]),
        lowPiece,
        minimum(concs[:, 5]),
        lowPiece,
        1e-9,
        1e-9,
        0.45,
        2,
        10,
        0,
        0,
    )
    highPiece = [10.0, 3.0, 3.0, 1.0, 1.0]
    high = vcat(
        maximum(concs[:, 1]),
        highPiece,
        maximum(concs[:, 2]),
        highPiece,
        maximum(concs[:, 3]),
        highPiece,
        maximum(concs[:, 4]),
        highPiece,
        maximum(concs[:, 5]),
        highPiece,
        3.0,
        3.0,
        0.55,
        60,
        180,
        50,
        50,
    )

    results_ode = bboptimize(
        hillCostAll;
        SearchRange = collect(zip(low, high)),
        NumDimensions = length(low),
        TraceMode = :verbose,
        TraceInterval = 100,
        MaxSteps = maxstep,
    )

    return best_fitness(results_ode), best_candidate(results_ode)
end

""" This function calculates the mean and std of the data for the three replicates. """
function find_mean_std_gs(g1s1, g1s2, g1s3, g2s1, g2s2, g2s3)
    
    meang1 = ones(189, 8, 5)
    meang2 = ones(189, 8, 5)
    stdg1 = ones(189, 8, 5)
    stdg2 = ones(189, 8, 5)

    for k=1:189
        meang1[k, :, :] .= mean([g1s1[k, :, :] , g1s2[k, :, :] , g1s3[k, :, :]])
        meang2[k, :, :] .= mean([g2s1[k, :, :] , g2s2[k, :, :] , g2s3[k, :, :]])
        stdg1[k, :, :] .= std([g1s1[k, :, :] , g1s2[k, :, :] , g1s3[k, :, :]])
        stdg2[k, :, :] .= std([g2s1[k, :, :] , g2s2[k, :, :] , g2s3[k, :, :]])
    end

    return meang1, meang2, stdg1, stdg2
end

function find_mean_std_simul(rep1, rep2, rep3, concs, g0)
    t = LinRange(0.0, 95.0, 189)
    p1 = getODEparamsAll(rep1, concs) 
    p2 = getODEparamsAll(rep2, concs)
    p3 = getODEparamsAll(rep3, concs)
    j=1; # drug number
    G1_1 = ones(189,8,5)
    G2_1 = ones(189,8,5)
    G1_2 = ones(189,8,5)
    G2_2 = ones(189,8,5)
    G1_3 = ones(189,8,5)
    G2_3 = ones(189,8,5)

    for j=1:5
        for i=1:8 # concentration number
            G1_1[:, i, j], G2_1[:, i, j], _ = predict(p1[:, i, j], g0, t)
            G1_2[:, i, j], G2_2[:, i, j], _ = predict(p2[:, i, j], g0, t)
            G1_3[:, i, j], G2_3[:, i, j], _ = predict(p3[:, i, j], g0, t)
        end
    end
    meanG1, meanG2, stdG1, stdG2 = find_mean_std_gs(G1_1, G1_2, G1_3, G2_1, G2_2, G2_3);
    return meanG1, meanG2, stdG1, stdG2
end

""" This function takes in the estimated Hill params for the 3 replicates, and outputs average and standard deviations of ODE model parameters"""
function avgRepsParams(rep1, rep2, rep3, concs)
    effs1 = getODEparamsAll(rep1, concs)
    effs2 = getODEparamsAll(rep2, concs)
    effs3 = getODEparamsAll(rep3, concs)
    # each of the following are 9x8 matrices
    lapat1 = effs1[:,:,1];
    dox1 = effs1[:,:,2];
    gemc1 = effs1[:,:,3];
    pac1 = effs1[:,:,4];
    pal1 = effs1[:,:,5];

    lapat2 = effs2[:,:,1];
    dox2 = effs2[:,:,2];
    gemc2 = effs2[:,:,3];
    pac2 = effs2[:,:,4];
    pal2 = effs2[:,:,5];

    lapat3 = effs3[:,:,1];
    dox3 = effs3[:,:,2];
    gemc3 = effs3[:,:,3];
    pac3 = effs3[:,:,4];
    pal3 = effs3[:,:,5];

    lapat = ones(9,8)
    dox = ones(9,8)
    gemc = ones(9,8)
    pac = ones(9,8)
    pal = ones(9,8)
    lapatstd = ones(9,8)
    doxstd = ones(9,8)
    gemcstd = ones(9,8)
    pacstd = ones(9,8)
    palstd = ones(9,8)

    for i=1:9
        lapat[i,:] .= mean([lapat1[i,:], lapat2[i,:], lapat[i,:]])
        dox[i,:] .= mean([dox1[i,:], dox2[i,:], dox3[i,:]])
        gemc[i,:] .= mean([gemc1[i,:], gemc2[i,:], gemc3[i,:]])
        pac[i,:] .= mean([pac1[i,:], pac2[i,:], pac3[i,:]])
        pal[i,:] .= mean([pal1[i,:], pal2[i,:], pal3[i,:]])

        lapatstd[i,:] .= std([lapat1[i,:], lapat2[i,:], lapat[i,:]])
        doxstd[i,:] .= std([dox1[i,:], dox2[i,:], dox3[i,:]])
        gemcstd[i,:] .= std([gemc1[i,:], gemc2[i,:], gemc3[i,:]])
        pacstd[i,:] .= std([pac1[i,:], pac2[i,:], pac3[i,:]])
        palstd[i,:] .= std([pal1[i,:], pal2[i,:], pal3[i,:]])
    end
    avgs = ones(9,8,5)
    stds = ones(9,8,5)
    avgs[:,:,1] = lapat; avgs[:,:,2] = dox; avgs[:,:,3] = gemc; avgs[:,:,4] = pac; avgs[:,:,5] = pal;
    stds[:,:,1] = lapatstd; stds[:,:,2] = doxstd; stds[:,:,3]=gemcstd; stds[:,:,4] = pacstd; stds[:,:,5]=palstd;
    return avgs, stds
end