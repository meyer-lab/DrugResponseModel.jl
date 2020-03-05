""" In this file we fit all the drugs att once. """

function getODEparamsAll(p::Array{Float64, 1}, concentrations::Array{Float64, 2})
    effects = zeros(9, length(concentrations[:, 1]), 4)

    k = 1
    # Scaled drug effect
    for i = 1:4
        xx = 1.0 ./ (1.0 .+ (p[k] ./ concentrations[:, i]) .^ p[k + 1])

        effects[1, :, i] = p[25] .+ (p[k + 2] - p[25]) .* xx
        effects[2, :, i] = p[26] .+ (p[k + 3] - p[26]) .* xx
        effects[3, :, i] = p[k + 4] .* xx
        effects[4, :, i] = p[k + 5] .* xx
        k += 6
    end
    effects[5, :, :] .= p[27]
    effects[6, :, :] .= floor(p[28]) #nG1
    effects[7, :, :] .= floor(p[29]) #nG2
    effects[8, :, :] .= floor(p[30]) #nD1
    effects[9, :, :] .= floor(p[31]) #nD2

    return effects
end

function residHillAll(hillParams::Array{Float64, 1}, concentrations::Array{Float64, 2}, g1::Array{Float64, 3}, g2::Array{Float64, 3})
    res = Atomic{eltype(hillParams)}(0.0)
    params = getODEparamsAll(hillParams, concentrations)

    # Solve for all drugs
    for j = 1:4
        @threads for ii = 1:length(concentrations[:, j])
            atomic_add!(
                res,
                cost(
                    params[:, ii, j],
                    g1[:, ii, j],
                    g2[:, ii, j],
                    Int(floor(params[6, ii, j])),
                    Int(floor(params[7, ii, j])),
                    Int(floor(params[8, ii, j])),
                    Int(floor(params[9, ii, j])),
                ),
            )
        end
    end

    return res[]
end

""" Hill optimization function for all drugs. """
function optimize_hillAll(concs::Array{Float64, 2}, g1::Array{Float64, 3}, g2::Array{Float64, 3}; maxstep = 1E5)
    hillCostAll(hillParams) = residHillAll(hillParams, concs, g1, g2)

    # The parameters used here in order:
     #(:Lap_EC50, :Lap_steepness, :Lap_maxG1ProgRate, :Lap_maxG2ProgRate, :Lap_maxDeathG1Rate, :Lap_maxDeathG2Rate, :Dox_EC50, :Dox_steepness, :Dox_maxG1ProgRate, :Dox_maxG2ProgRate, :Dox_maxDeathG1Rate, :Dox_maxDeathG2Rate, :Gem_EC50, :Gem_steepness, :Gem_maxG1ProgRate, :Gem_maxG2ProgRate, :Gem_maxDeathG1Rate, :Gem_maxDeathG2Rate, :Tax_EC50, :Tax_steepness, :Tax_maxG1ProgRate, :Tax_maxG2ProgRate, :Tax_maxDeathG1Rate, :Tax_maxDeathG2Rate, :G1ProgRateControl, :G2ProgRateControl, :percG1, :nG1, :nG2, :nD1, :nD2)
    low = [
        minimum(concs[:, 1]),
        0.01,
        1e-9,
        1e-9,
        0.0,
        0.0,
        minimum(concs[:, 2]),
        0.01,
        1e-9,
        1e-9,
        0.0,
        0.0,
        minimum(concs[:, 3]),
        0.01,
        1e-9,
        1e-9,
        0.0,
        0.0,
        minimum(concs[:, 4]),
        0.01,
        1e-9,
        1e-9,
        0.0,
        0.0,
        1e-9,
        1e-9,
        0.45,
        2,
        10,
        0,
        0,
    ]
    high = [
        maximum(concs[:, 1]),
        10.0,
        3.0,
        3.0,
        1.0,
        1.0,
        maximum(concs[:, 2]),
        10.0,
        3.0,
        3.0,
        1.0,
        1.0,
        maximum(concs[:, 3]),
        10.0,
        3.0,
        1.0,
        1.0,
        1.0,
        maximum(concs[:, 4]),
        10.0,
        3.0,
        3.0,
        1.0,
        1.0,
        3.0,
        3.0,
        0.55,
        60,
        180,
        50,
        50,
    ]

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


""" Combination functions. """
function ParamForBliss(p,n)
    """ To calculate Bliss independence drug effect
    we assume delays are constant, death rates are additive,
    and will keep the alpha and beta intact."""
    par = zeros(4,n)
    par[1,:] = p[1,1] .- p[1,:] # alpha stays the same
    par[2,:] = p[2,1] .- p[2,:] # beta stays the same
    par[3,:] = p[3,1] .- p[3,:] # death rate in G1
    par[4,:] = p[4,1] .- p[4,:] # death rate in G2
    return par
end

function BlissCombination(p1::Matrix{Float64}, p2::Matrix{Float64}, n::Int)
    """ A function to calculate Bliss independence for drug combination assuming
    the two drugs hit different pathways and they effect independently. """

    param1 = ParamForBliss(p1,n)
    param2 = ParamForBliss(p2,n)
    """ For 8x8 combination of drug concentrations for G1 progression rate, G2 progression rate, and death rates in G1 and G2, respectively. """
    combined = zeros(n,n,4)
    for j in 1:n
        for k in 1:n
            combined[j,k,1:2] .= -(param1[1:2,j] .+ param2[1:2,k] .- param1[1:2,j] .* param2[1:2,k]) .+ p1[1:2,1,1]
            combined[j,k,3:4] .= -(param1[3:4,j] .+ param2[3:4,k])
            end
        end
    return combined
end

""" To output the full ODE params for plotting the cell number. """
function fullCombinationParam(origP1, origP2, origFullParam,n)
    """ Here we assume the base is origP1, and we just want to get the params of EC50 from origP2. """
    combined = BlissCombination(origP1, origP2,n)
    fullparam = zeros(9,n,n)
    fullparam[5:9, :, :] .= origFullParam[5:9, 1, 1]
    for i = 1:4
        fullparam[i, :, :] .= combined[:, :, i]
    end
    return fullparam
end

""" Function unit to plot drug effects before and after combination. """
function plotunitCombin(conc, gemc, titles, combin)
    concs = log.(conc)
    plot(concs, gemc, ylabel=titles, label = "taxol alone", legendfontsize = 7, lw = 3, fg_legend = :transparent, shape=:circle, color=:purple)
    plot!(concs, combin, label = "taxol w/ 5nM gemc.", lw=3, shape=:circle, color=:green) 
end

""" Function to plot all of the drug effects before and after drug combination. """
function plotEffectsCombin(concs, gemc, combin)
    titles = ["G1 prog. rate", "G2 prog. rate", "G1 death rate", "G2 death rate"]
    pl = [plotunitCombin(concs[:, 4], gemc[i, :], titles[i], combin[:, 4, i]) for i = 1:4]
    plot(pl..., layout = (2, 2))
    plot!(size = (800, 500), margin = 0.4cm, dpi = 150)
end

function plotNumcells(drugB, combination, concDrugB, g0, n)
    numscomb = zeros(n)
    nums = zeros(n)
    for j =1:n
        numscomb[j] = numcells(combination[:, j], g0, 96)
        nums[j] = numcells(drugB[:, j], g0, 96)
    end
    plot(log.(concDrugB), numscomb, label = "pac + gemc", legendfontsize = 7, lw = 3, fg_legend = :transparent, shape=:circle, color=:purple)
    plot!(log.(concDrugB), nums, label="pac", lw=3, xlabel="log drug concentration", ylabel="cell #", shape=:circle, color=:green)
    plot!(dpi=150)
end

""" Plotting cell# versus concentration for2 drugs """
function combin2drugs(d1, d2, concd1, concd2, named1, named2, effs, g0)  
    combin = fullCombinationParam(d1, d2, effs,8);
    n=8;

    numscomb = zeros(8,8)
    nums = zeros(n)
    for j =1:n
        nums[j] = numcells(d1[:, j], g0, 96)
        for m=1:8
            numscomb[j,m] = numcells(combin[:, j, m], g0, 96)
        end
    end
    plot(log.(concd1), nums, label=string(named1), lw=3, xlabel="log drug concentration", ylabel="cell #", shape=:circle, color=:green)
    for k = 2:8
        plot!(log.(concd1), numscomb[:, k], label = string(named1, " +", Int(concd2[k,1]), "nM ", named2), legendfontsize = 7, lw = 3, fg_legend = :transparent, shape=:circle, color=:purple, alpha = (1-0.1*k))
    end
end