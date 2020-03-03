""" In this file we fit all the drugs att once. """

function getODEparamsAll(p::Array{Float64, 1}, concentrations::Array{Float64, 2})
    effects = effects = zeros(9, 8, 4)

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
        3.0,
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
function ParamForBliss(p)
    """ To calculate Bliss independence drug effect
    we assume delays are constant, death rates are additive,
    and will keep the alpha and beta intact."""
    par = zeros(4,8)
    par[1,:] = p[1,:] # alpha stays the same
    par[2,:] = p[2,:] # beta stays the same
    par[3,:] = p[3,:] # death rate in G1
    par[4,:] = p[4,:] # death rate in G2
    @assert(all(par .>= 0.0), " the drug effects <= 0")
    return par
end

function BlissCombination(p1::Matrix{Float64}, p2::Matrix{Float64})
    """ A function to calculate Bliss independence for drug combination assuming
    the two drugs hit different pathways and they effect independently. """

    param1 = ParamForBliss(p1)
    param2 = ParamForBliss(p2)
    """ For 8x8 combination of drug concentrations for G1 progression rate, G2 progression rate, and death rates in G1 and G2, respectively. """
    combined = zeros(8,8,4)
    for j in 1:8
        for k in 1:8
            combined[j,k,1:2] .= param1[1:2,j] .+ param2[1:2,k] .- param1[1:2,j] .* param2[1:2,k]
            combined[j,k,3:4] .= param1[3:4,j] .+ param2[3:4,k]
            end
        end
    return combined
end

""" To output the full ODE params for plotting the cell number. """
function fullCombinationParam(combined, origFullParam)
    fullparam = zeros(9,8,8)

    fullparam[5:9, :, :] .= origFullParam[5:9, 1, 1]
    for i=1:4
        fullparam[i, :, :] .= combined[:, :, i]
    end
    return fullparam
end

function plotCombinODE(params, g0, title, ymax)
    t = LinRange(0.0, 120, 200)
    G1, G2 = predict(params, g0, t, Int(floor(params[6])), Int(floor(params[7])), Int(floor(params[8])), Int(floor(params[9])))

    plot(t,
        G1,
        label = "G1 est",
        xlabel = "time [hours]",
        ylabel = "# of cells",
        xguidefontsize = 8,
        yguidefontsize = 8,
        lw = 2.0,
        alpha = 0.6,
        color = :green,
    )
    plot!(t, G2, label = "G2 est", legend = :topleft, legendfontsize = 4, fg_legend = :transparent, lw = 2.0, alpha = 0.6, color = :sienna)
    plot!(t, G1 .+ G2, label = "total est", dpi = 150, lw = 2.0, alpha = 0.6, color = :hotpink)
    plot!(annotation = [(60, ymax, text(title, 8))])
    ylims!((0.0, ymax))
end

function combinplot_all(params_ode, g0, conc::Array{Float64, 1})
    # plotting the fitted curves
    rl = [plotCombinODE(params_ode[:, i], g0, string(conc[i], " nM"), 80.0) for i = 1:4]
    r2 = [plotCombinODE(params_ode[:, i], g0, string(conc[i], " nM"), 40.0) for i = 5:7]
    r8 = plotCombinODE(params_ode[:, 8], g0, string(conc[8], " nM"), 40.0)
    plot(rl..., r2..., r8, layout = (2, 4))
    plot!(size = (900, 400), margin = 0.4cm, dpi = 200)
end

function plotunitCombin(conc, gemc, titles, combin)
    concs = log.(conc)
    plot(concs, gemc, ylabel=titles, label = "taxol alone", legendfontsize = 7, lw = 3, fg_legend = :transparent)
    plot!(concs, combin, label = "taxol w/ 5nM gemc.", lw=3) 
end
function plotEffectsCombin(concs, gemc, combin)
    titles = ["G1 prog. rate", "G2 prog. rate", "G1 death rate", "G2 death rate"]
    pl = [plotunitCombin(concs[:, 4], gemc[i, :], titles[i], combin[:, 4, i]) for i = 1:4]
    plot(pl..., layout = (2,2))
    plot!(size = (800, 500), margin = 0.4cm, dpi=150)
end