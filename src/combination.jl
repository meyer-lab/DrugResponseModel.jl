""" This file contains all the functions related to Bliss and temporal combinations. """
function ParamForBliss(p::Matrix{Float64}, n::Int)
    """ To calculate Bliss independence drug effect
    we assume delays are constant, death rates are additive,
    and will keep the alpha and beta intact."""
    par = zeros(4, n)
    par[1, :] = p[1, 1] .- p[1, :] # alpha stays the same
    par[2, :] = p[2, 1] .- p[2, :] # beta stays the same
    par[3, :] = p[3, 1] .- p[3, :] # death rate in G1
    par[4, :] = p[4, 1] .- p[4, :] # death rate in G2
    return par
end

function BlissCombination(p1::Array{Float64, 2}, p2::Array{Float64, 2}, n::Int)
    """ A function to calculate Bliss independence for drug combination assuming
    the two drugs hit different pathways and they effect independently. """

    param1 = ParamForBliss(p1, n)
    param2 = ParamForBliss(p2, n)
    """ For 8x8 combination of drug concentrations for G1 progression rate, G2 progression rate, and death rates in G1 and G2, respectively. """
    combined = zeros(n, n, 4)
    for j = 1:n
        for k = 1:n
            combined[j, k, 1:2] .= -(param1[1:2, j] .+ param2[1:2, k] .- param1[1:2, j] .* param2[1:2, k]) .+ p1[1:2, 1, 1]
            combined[j, k, 3:4] .= -(param1[3:4, j] .+ param2[3:4, k])
        end
    end
    @assert(all(combined .>= 0.0))
    return combined
end

""" To output the full ODE params for plotting the cell number. """
function fullCombinationParam(origP1::Array{Float64, 2}, origP2::Array{Float64, 2}, origFullParam::Array{Float64, 3}, n::Int)
    """ Here we assume the base is origP1, and we just want to get the params of EC50 from origP2. """
    combined = BlissCombination(origP1, origP2, n)
    fullparam = zeros(9, n, n)
    fullparam[5:9, :, :] .= origFullParam[5:9, 1, 1]
    for i = 1:4
        fullparam[i, :, :] .= combined[:, :, i]
    end
    return fullparam
end

""" Function unit to plot drug effects before and after combination. """
function plotunitCombin(conc::Array{Float64, 1}, gemc::Array{Float64, 1}, titles, combin::Array{Float64, 1})
    concs = log.(conc)
    plot(concs, gemc, ylabel = titles, label = "taxol alone", lw = 3, fg_legend = :transparent, shape = :circle, color = :purple)
    plot!(concs, combin, label = "taxol w/ 5nM gemc.", lw = 3, shape = :circle, color = :green)
end

""" Function to plot all of the drug effects before and after drug combination. """
function plotEffectsCombin(concs::Array{Float64, 2}, gemc::Array{Float64, 2}, combin::Array{Float64, 3})
    titles = ["G1 prog. rate", "G2 prog. rate", "G1 death rate", "G2 death rate"]
    pl = [plotunitCombin(concs[:, 4], gemc[i, :], titles[i], combin[:, 4, i]) for i = 1:4]
    plot(pl..., layout = (2, 2))
end

function plotNumcells(drugB::Array{Float64, 2}, combination::Array{Float64, 2}, concDrugB::Array{Float64, 1}, g0::Float64, n::Int)
    numscomb = zeros(n)
    nums = zeros(n)
    for j = 1:n
        numscomb[j] = numcells(combination[:, j], g0, 96)
        nums[j] = numcells(drugB[:, j], g0, 96)
    end
    plot(log.(concDrugB), numscomb, label = "pac + gemc", lw = 3, fg_legend = :transparent, shape = :circle, color = :purple)
    plot!(log.(concDrugB), nums, label = "pac", lw = 3, xlabel = "log drug concentration", ylabel = "cell #", shape = :circle, color = :green)
end

function helperPlot(concd1, named1, concd2, named2, numscomb)
    p = plot(
        log.(concd1),
        numscomb[:, 1],
        label = string(named1),
        lw = 3,
        xlabel = "log drug concentration",
        ylabel = "cell # at t = 96 hrs",
        shape = :circle,
        color = :green,
    )
    for k = 2:8
        plot!(
            log.(concd1),
            numscomb[:, k],
            label = string(named1, " +", concd2[k, 1], "nM ", named2),
            lw = 3,
            fg_legend = :transparent,
            shape = :circle,
            color = :purple,
            alpha = (1 - 0.1 * k),
            show = true,
        )
    end
    p
end
""" Plotting cell# versus concentration for2 drugs """
function combin2drugs(
    d1::Array{Float64, 2},
    d2::Array{Float64, 2},
    concd1::Array{Float64, 1},
    concd2::Array{Float64, 1},
    named1::String,
    named2::String,
    effs::Array{Float64, 3},
    g0::Float64,
)
    n = 8
    combin = fullCombinationParam(d1, d2, effs, n)

    numscomb = zeros(n, n)
    for j = 1:n
        for m = 1:n
            numscomb[j, m] = numcells(combin[:, j, m], g0, 96)
        end
    end
    helperPlot(concd1, named1, concd2, named2, numscomb)
end

""" In this function, we apply the Bliss synergy to the cell numbers. 
This is to compare the traditional way of representing the combination effect, compare to the way we do in our model."""
function blissCellNum(g1s, g2s; T = 96, n = 8)
    num = zeros(n, 4)
    # for no specific reason, I chose lapatinib's control trial to be the base case for converting.
    base = g1s[T, 1, 1] + g2s[T, 1, 1]
    for i = 1:4
        # num is a 8 x 4 matrix, holding cell numbers for 4 drugs, in 8 concenntration, for a specific time point.
        num[:, i] = 1.0 .- ((g1s[T, :, i] + g2s[T, :, i]) ./ base)
    end
    combined = zeros(n, n, 6)
    for j = 1:n
        combined[j, :, 1] = -(num[:, 1] .+ num[j, 2] .- (num[:, 1] .* num[j, 2]) .- 1.0) .* base # lap w/ dox
        combined[j, :, 2] = -(num[:, 1] .+ num[j, 3] .- (num[:, 1] .* num[j, 3]) .- 1.0) .* base # lap w/ gem
        combined[j, :, 3] = -(num[:, 1] .+ num[j, 4] .- (num[:, 1] .* num[j, 4]) .- 1.0) .* base # lap w/ pac
        combined[j, :, 4] = -(num[:, 2] .+ num[j, 3] .- (num[:, 2] .* num[j, 3]) .- 1.0) .* base # dox w/ gem
        combined[j, :, 5] = -(num[:, 2] .+ num[j, 4] .- (num[:, 2] .* num[j, 4]) .- 1.0) .* base # dox w/ pac
        combined[j, :, 6] = -(num[:, 3] .+ num[j, 4] .- (num[:, 3] .* num[j, 4]) .- 1.0) .* base # gem w/ pac
    end
    @assert(all(combined .>= 0.0))
    return combined
end

""" Function for calculating temporal combination of two drugs. """
function temporal_combination(params1, params2, g0)
    t1 = LinRange(0.0, 60.0, 100)

    g1L, g2L, vecL = predict(params1, g0, t1, Int(floor(params1[6])), Int(floor(params1[7])), Int(floor(params1[8])), Int(floor(params1[9])))
    g1G, g2G, _ = predict(params2, vec(vecL), t1, Int(floor(params2[6])), Int(floor(params2[7])), Int(floor(params2[8])), Int(floor(params2[9])))

    return vcat(g1L, g1G), vcat(g2L, g2G)
end

function helperPlotCombin(G1, G2, g0, title::String, legend::Any, ymax)
    t_new = LinRange(0.0, 120, 200)
    plot(
        t_new,
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
    plot!(t_new, G2, label = "G2 est", legend = legend, legendfontsize = 4, fg_legend = :transparent, lw = 2.0, alpha = 0.6, color = :sienna)
    plot!(t_new, G1 .+ G2, label = "total est", dpi = 150, lw = 2.0, alpha = 0.6, color = :hotpink)
    plot!(annotation = [(60, ymax, text(title, 8))])
    ylims!((0.0, ymax))
end

""" Function to plot temporal combinations of two drugs. """
function plotTemporalCombin(params1, params2, g1s, g2s, pop, concl, concg)
    # This is right now specificly for lapatinib and gemcitabine
    i = 6 # 6th concentration of lapatinib
    j = 5 # 6th concentration of gemcitabine
    G1_1, G2_1 = temporal_combination(params1, params2, g1s[1,1,1]+g2s[1,1,1])
    G1_2, G2_2 = temporal_combination(params2, params1, g1s[1,1,1]+g2s[1,1,1])
    p1 = ode_plotIt(params1, g1s[:, :, 1], g2s[:, :, 1], pop[1], i, string(concl[i]," nM lap."), false, 70.0)
    p2 = ode_plotIt(params2, g1s[:, :, 3], g2s[:, :, 3], pop[3], j, string(concg[j]," nM Gemc."), true, 70.0)
    p3 = helperPlotCombin(G1_1, G2_1, g1s[1,1,1]+g2s[1,1,1], string(concl[i], " nM Lap+ ", concg[j], "nM Gemc"), false, 70.0) # first lapatinib, then gemcitabine
    p4 = helperPlotCombin(G1_2, G2_2, g1s[1,1,1]+g2s[1,1,1], string(concg[j], " nM Lap+ ", concl[i], "nM Gemc"), false, 70.0) # first gemcitabine then lapatinib
    plot(p1, p2, p3, p4, layout=(2,2))
end