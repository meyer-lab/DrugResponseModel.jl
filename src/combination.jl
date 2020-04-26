""" This file contains all the functions related to Bliss and temporal combinations. """

function BlissCombination(p1::Array{Float64, 2}, p2::Array{Float64, 2}, n::Int)
    """ A function to calculate Bliss independence for drug combination assuming
    the two drugs hit different pathways and they effect independently. """

    param1 = zeros(4, 8)
    param2 = zeros(4, 8)
    param1[1, :] .= p1[1, 1] .- p1[1, :]
    param1[2, :] .= p1[2, 1] .- p1[2, :]
    param1[3, :] .= p1[3, :]
    param1[4, :] .= p1[4, :]
    param2[1, :] .= p2[1, 1] .- p2[1, :]
    param2[2, :] .= p2[2, 1] .- p2[2, :]
    param2[3, :] .= p2[3, :]
    param2[4, :] .= p2[4, :]

    """ For 8x8 combination of drug concentrations for G1 progression rate, G2 progression rate, and death rates in G1 and G2, respectively. """
    combined = zeros(n, n, 4)
    for j = 1:n
        for k = 1:n
            combined[j, k, 1:2] .= -(param1[1:2, j] .+ param2[1:2, k] .- param1[1:2, j] .* param2[1:2, k]) .+ p1[1:2, 1, 1]
            combined[j, k, 3:4] .= param1[3:4, j] .+ param2[3:4, k]
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

function helperPlot(concd1, named1, concd2, named2, numscomb, legend, title, ymin, ymax)
    p = plot(
        log.(concd1),
        numscomb[:, 1],
        label = string(named1),
        lw = 3,
        xlabel = "log drug concentration",
        ylabel = "cell # at t = 96 hrs",
        shape = :circle,
        color = :green,
        title = title,
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
            legend = legend,
            show = true,
        )
    end
    ylims!((ymin, ymax))
    p
end


""" In this function, we apply the Bliss synergy to the cell numbers. 
This is to compare the traditional way of representing the combination effect, compare to the way we do in our model."""
function blissCellNum(g1s, g2s; T = 96, n = 8)
    num = zeros(n, 5)
    # for no specific reason, I chose lapatinib's control trial to be the base case for converting.
    base = g1s[T, 1, 1] + g2s[T, 1, 1]
    for i = 1:5
        # num is a 8 x 4 matrix, holding cell numbers for 4 drugs, in 8 concenntration, for a specific time point.
        num[:, i] = 1.0 .- ((g1s[T, :, i] + g2s[T, :, i]) ./ base)
    end
    combined = zeros(n, n, 10)
    for j = 1:n
        combined[j, :, 1] = -(num[:, 1] .+ num[j, 2] .- (num[:, 1] .* num[j, 2]) .- 1.0) .* base # lap w/ dox; meaning dox changes with rows and lap changes with columns
        combined[j, :, 2] = -(num[:, 1] .+ num[j, 3] .- (num[:, 1] .* num[j, 3]) .- 1.0) .* base # lap w/ gem
        combined[j, :, 3] = -(num[:, 1] .+ num[j, 4] .- (num[:, 1] .* num[j, 4]) .- 1.0) .* base # lap w/ pac
        combined[j, :, 4] = -(num[:, 1] .+ num[j, 5] .- (num[:, 1] .* num[j, 5]) .- 1.0) .* base # lap w/ palb
        combined[j, :, 5] = -(num[:, 2] .+ num[j, 3] .- (num[:, 2] .* num[j, 3]) .- 1.0) .* base # dox w/ gem
        combined[j, :, 6] = -(num[:, 2] .+ num[j, 4] .- (num[:, 2] .* num[j, 4]) .- 1.0) .* base # dox w/ pac
        combined[j, :, 7] = -(num[:, 2] .+ num[j, 5] .- (num[:, 2] .* num[j, 5]) .- 1.0) .* base # dox w/ palb
        combined[j, :, 8] = -(num[:, 3] .+ num[j, 4] .- (num[:, 3] .* num[j, 4]) .- 1.0) .* base # gem w/ pac
        combined[j, :, 9] = -(num[:, 3] .+ num[j, 5] .- (num[:, 3] .* num[j, 5]) .- 1.0) .* base # gem w/ palb
        combined[j, :, 10] = -(num[:, 4] .+ num[j, 5] .- (num[:, 4] .* num[j, 5]) .- 1.0) .* base # pac w/ palb
    end
    @assert(all(combined .>= 0.0))
    return combined
end

""" Function for calculating temporal combination of two drugs. """
function temporal_combination(params1, params2, g0::Float64, max1::Float64, max2::Float64)
    t1 = LinRange(0.0, max1, Int(max1*2))
    t2 = LinRange(0.0, max2, Int(max2*2))

    g1L, g2L, vecL = predict(params1, g0, t1)
    g1G, g2G, _ = predict(params2, vec(vecL), t2)

    return vcat(g1L, g1G), vcat(g2L, g2G)
end

function helperPlotCombin(G1, G2, g0::Float64, title::String, legend::Any, ymax::Float64)
    t_new = LinRange(0.0, 90.0, 180)
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
    plot!(t_new, G1 .+ G2, label = "total est", lw = 2.0, alpha = 0.6, color = :hotpink)
    plot!(annotation = [(100, ymax, text(title, 8))])
    ylims!((0.0, ymax))
end

""" Function to plot temporal combinations of two drugs. """
function plotTemporalCombin(params1, params2, g1s, g2s, pop, concl, concg, legend, i, j, named1, named2, k1, k2)
    # This is right now specificly for lapatinib and doxorubicin
    # ith concentration of lapatinib
    # jth concentration of doxorubicin
    diff = find_combin_order(params1, params2, g1s, g2s)
    max1 = max2 = 45.0
    G1_1, G2_1 = temporal_combination(params1, params2, g1s[1, 1, 1] + g2s[1, 1, 1], max1, max2)
    G1_2, G2_2 = temporal_combination(params2, params1, g1s[1, 1, 1] + g2s[1, 1, 1], max1, max2)
    p1 = ode_plotIt(params1, g1s[:, :, k1], g2s[:, :, k1], pop[:, :, k1], i, string(concl[i], " nM ", named1), false, 70.0; tnew=90)
    p2 = ode_plotIt(params2, g1s[:, :, k2], g2s[:, :, k2], pop[:, :, k2], j, string(concg[j], " nM ", named2), false, 70.0; tnew=90)
    p3 = helperPlotCombin(G1_1, G2_1, g1s[1, 1, 1] + g2s[1, 1, 1], string(concl[i], " nM ", named1, "+", concg[j], "nM ", named2), legend, 70.0) # first lapatinib, then gemcitabine
    p4 = helperPlotCombin(G1_2, G2_2, g1s[1, 1, 1] + g2s[1, 1, 1], string(concg[j], " nM ", named2, "+", concl[i], "nM ", named1), false, 70.0) # first gemcitabine then lapatinib
    plot(p1, p2, p3, p4, layout = (2, 2))
end

""" To find IC50 for each drug, separately."""
function find_IC50(population)
    lap = Array(population[189, :, 1])
    dox = Array(population[189, :, 2])
    gem = Array(population[189, :, 3])
    tax = Array(population[189, :, 4])
    pal = Array(population[189, :, 5])
    IC50_lap = argmin(abs.(0.5 * lap[1] .- lap)) #6
    IC50_dox = argmin(abs.(0.5 * dox[1] .- dox)) #3
    IC50_gem = argmin(abs.(0.5 * gem[1] .- gem)) #6
    IC50_tax = argmin(abs.(0.5 * tax[1] .- tax)) #4
    IC50_pal = argmin(abs.(0.5 * pal[1] .- pal)) #5
    return (IC50_lap, IC50_dox, IC50_gem, IC50_tax, IC50_pal)
end


function heatmap_combination(d1, d2, cellNum, i1, i2, d1name, d2name, effs, concs, g0)
    n = 8 # the number of concentrations we have
    combin = fullCombinationParam(d1, d2, effs, n)

    numscomb = zeros(n, n)
    for j = 1:n
        for m = 1:n
            numscomb[j, m] = numcells(combin[:, j, m], g0, 96)
        end
    end

    diffs = numscomb ./ cellNum # model prediction / reference
    concs[1, :] .= 0.6
    heatmap(
        string.(round.(log.(concs[:, i2]), digits = 1)),
        string.(round.(log.(concs[:, i1]), digits = 1)),
        diffs,
        xlabel = string(d2name, " log[nM]"),
        ylabel = string(d1name, " log [nM]"),
        title = "cell number fold diff",
        clim = (0.0, 2.0),
    )
end

""" find those conbiations in which order of treatment matters. """
function find_combin_order(params1, params2, g1s, g2s)

    diff = zeros(17, 17)
    i=1
    j=1
    for max1 = 10.0:5:90.0 # 5-hour interval, from 10 hours to 90 hours
        for max2 = 10.0:5:90.0
            G1_1, G2_1 = temporal_combination(params1, params2, g1s[1, 1, 1] + g2s[1, 1, 1], max1, max2)
            G1_2, G2_2 = temporal_combination(params2, params1, g1s[1, 1, 1] + g2s[1, 1, 1], max1, max2)
            diff[i,j] = abs((G1_1[end] + G2_1[end]) - (G1_2[end] + G2_2[end]))
            j += 1
        end
        j = 1
        i += 1
    end
    return diff
end

""" Plot the heatmap to describe the difference between the order of treatments. """
function plot_order_temporalCombin(params1, params2, g1s, g2s, named1, named2)
    diffs = DrugResponseModel.find_combin_order(params1, params2, g1s, g2s)
    heatmap(string.(10.0:5:90.0),
        string.(10.0:5:90.0),
        diffs,
        xlabel = string("time max2 [hr]"),
        ylabel = string("time max1 [hr]"),
        title = string("cell number diff", named1, "->", named2, "-", named2, "->", named1),
        clim = (0.0, 15.0),
    )
end