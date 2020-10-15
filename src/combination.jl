""" This file contains all the functions related to Bliss and temporal combinations. """

function BlissCombination(p1::Array{Float64, 2}, p2::Array{Float64, 2}, n::Int)
    """ A function to calculate Bliss independence for drug combination assuming
    the two drugs hit different pathways and they effect independently. """

    param1 = zeros(4, 8)
    param2 = zeros(4, 8)
    param1[1, :] .= 1.0 .- (p1[1, :] ./ p1[1, 1])
    param1[2, :] .= 1.0 .- (p1[2, :] ./ p1[2, 1])
    param1[3, :] .= p1[3, :]
    param1[4, :] .= p1[4, :]
    param2[1, :] .= 1.0 .- (p2[1, :] ./ p2[1, 1])
    param2[2, :] .= 1.0 .- (p2[2, :] ./ p2[2, 1])
    param2[3, :] .= p2[3, :]
    param2[4, :] .= p2[4, :]

    """ For 8x8 combination of drug concentrations for G1 progression rate, G2 progression rate, and death rates in G1 and G2, respectively. """
    combined = zeros(n, n, 4)
    for j = 1:n
        for k = 1:n
            combined[j, k, 1:2] .= (1.0 .- (param1[1:2, j] .+ param2[1:2, k] .- param1[1:2, j] .* param2[1:2, k])) .* p1[1:2, 1, 1]
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
        # num is a 8 x 5 matrix, holding cell numbers for 5 drugs, in 8 concenntration, for a specific time point.
        num[:, i] = 1.0 .- ((g1s[T, :, i] + g2s[T, :, i]) ./ base)
    end
    combined = zeros(n, n, 10)
    for j = 1:n
        combined[:, j, 1] = -(num[:, 1] .+ num[j, 2] .- (num[:, 1] .* num[j, 2]) .- 1.0) .* base # lap w/ dox; meaning dox changes with rows and lap changes with columns
        combined[:, j, 2] = -(num[:, 1] .+ num[j, 3] .- (num[:, 1] .* num[j, 3]) .- 1.0) .* base # lap w/ gem
        combined[:, j, 3] = -(num[:, 1] .+ num[j, 4] .- (num[:, 1] .* num[j, 4]) .- 1.0) .* base # lap w/ pac
        combined[:, j, 4] = -(num[:, 1] .+ num[j, 5] .- (num[:, 1] .* num[j, 5]) .- 1.0) .* base # lap w/ palb
        combined[:, j, 5] = -(num[:, 2] .+ num[j, 3] .- (num[:, 2] .* num[j, 3]) .- 1.0) .* base # dox w/ gem
        combined[:, j, 6] = -(num[:, 2] .+ num[j, 4] .- (num[:, 2] .* num[j, 4]) .- 1.0) .* base # dox w/ pac
        combined[:, j, 7] = -(num[:, 2] .+ num[j, 5] .- (num[:, 2] .* num[j, 5]) .- 1.0) .* base # dox w/ palb
        combined[:, j, 8] = -(num[:, 3] .+ num[j, 4] .- (num[:, 3] .* num[j, 4]) .- 1.0) .* base # gem w/ pac
        combined[:, j, 9] = -(num[:, 3] .+ num[j, 5] .- (num[:, 3] .* num[j, 5]) .- 1.0) .* base # gem w/ palb
        combined[:, j, 10] = -(num[:, 4] .+ num[j, 5] .- (num[:, 4] .* num[j, 5]) .- 1.0) .* base # pac w/ palb
    end
    @assert(all(combined .>= 0.0))
    return combined
end

""" Function for calculating temporal combination of two drugs. """
function temporal_combination(params1, params2, g0::Float64, max1::Float64, max2::Float64)
    t1 = LinRange(0.0, max1, Int(max1 * 2))
    t2 = LinRange(0.0, max2, Int(max2 * 2))

    g1L, g2L, vecL = predict(params1, g0, t1)
    g1G, g2G, _ = predict(params2, vec(vecL), t2)

    return vcat(g1L, g1G), vcat(g2L, g2G)
end

function helperPlotCombin(G1, G2, g0::Float64, title::String, legend::Any, ymax::Float64)
    t_new = LinRange(0.0, length(G1), Int(length(G1)))
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

""" Function to plot temporal combinations of two drugs. When you are given the specific concentrations for drugA and drugB,
This returns the temporal combination with particular order that has the most difference between drugA then drug B, or drugB then drugA."""
function plotTemporalCombin(params1, params2, g1s, g2s, pop, concl, concg, legend, i, j, named1, named2, k1, k2)
    # Let's say this is for lapatinib and doxorubicin
    # ith concentration of lapatinib
    # jth concentration of doxorubicin
    diff = find_combin_order(params1, params2, g1s, g2s)
    rowmaxind = zeros(length(diff[1, :]))
    rowmaxval = zeros(length(diff[1, :]))
    for j = 1:length(diff[1, :])
        rowmaxval[j], rowmaxind[j] = findmax(abs.(diff[j, :])) # specifies column
    end
    maxdiff, index = findmax(rowmaxval) # specifies row
    tim = 10.0:5:90.0
    max1 = tim[index]
    max2 = tim[Int(rowmaxind[index])]
    t_new = LinRange(0.0, max1 + max2, Int(2 * (max1 + max2)))
    G1_1, G2_1 = temporal_combination(params1, params2, g1s[1, 1, 1] + g2s[1, 1, 1], max1, max2)
    G1_2, G2_2 = temporal_combination(params2, params1, g1s[1, 1, 1] + g2s[1, 1, 1], max1, max2)
    p1 = ode_plotIt(params1, g1s[:, :, k1], g2s[:, :, k1], pop[:, :, k1], i, string(concl[i], " nM ", named1), false, 70.0, t_new)
    p2 = ode_plotIt(params2, g1s[:, :, k2], g2s[:, :, k2], pop[:, :, k2], j, string(concg[j], " nM ", named2), false, 70.0, t_new)
    p3 = helperPlotCombin(G1_1, G2_1, g1s[1, 1, 1] + g2s[1, 1, 1], string(concl[i], " nM ", named1, "+", concg[j], "nM ", named2), legend, 70.0) # first lapatinib, then gemcitabine
    p4 = helperPlotCombin(G1_2, G2_2, g1s[1, 1, 1] + g2s[1, 1, 1], string(concg[j], " nM ", named2, "+", concl[i], "nM ", named1), false, 70.0) # first gemcitabine then lapatinib
    plot(p1, p2, p3, p4, layout = (2, 2))
end

""" find the difference matrix betwrrn two orders of treatment. """
function find_combin_order(params1, params2, g1s, g2s)

    diff = zeros(17, 17)
    i = 1
    j = 1
    for max1 = 10.0:5:90.0 # 5-hour interval, from 10 hours to 90 hours
        for max2 = 10.0:5:90.0
            G1_1, G2_1 = temporal_combination(params1, params2, g1s[1, 1, 1] + g2s[1, 1, 1], max1, max2)
            G1_2, G2_2 = temporal_combination(params2, params1, g1s[1, 1, 1] + g2s[1, 1, 1], max1, max2)
            diff[i, j] = (G1_1[end] + G2_1[end]) - (G1_2[end] + G2_2[end])
            j += 1
        end
        j = 1
        i += 1
    end
    return diff
end

""" To find IC50 or IC90 for each drug, separately."""
function find_IC(population, which)
    lap = Array(population[189, :, 1])
    dox = Array(population[189, :, 2])
    gem = Array(population[189, :, 3])
    tax = Array(population[189, :, 4])
    pal = Array(population[189, :, 5])
    IC_lap = argmin(abs.(which * lap[1] .- lap)) #6
    IC_dox = argmin(abs.(which * dox[1] .- dox)) #3
    IC_gem = argmin(abs.(which * gem[1] .- gem)) #6
    IC_tax = argmin(abs.(which * tax[1] .- tax)) #4
    IC_pal = argmin(abs.(which * pal[1] .- pal)) #5
    return IC_lap, IC_dox, IC_gem, IC_tax, IC_pal # returns the argument
end

""" Plot the heatmap to describe the difference between the order of treatments. """
function plot_order_temporalCombin(params1, params2, g1s, g2s, named1, named2)
    diffs = DrugResponseModel.find_combin_order(params1, params2, g1s, g2s)
    heatmap(
        string.(10.0:5:90.0),
        string.(10.0:5:90.0),
        diffs,
        xlabel = string("time max2 [hr]"),
        ylabel = string("time max1 [hr]"),
        title = string(named1, " --> ", named2, " - ", named2, " --> ", named1),
    )
end

####--------------- Loewe additivity ------------------####
""" Find the inverse of a hill function (concentration), given the parameters and the effect. """
function inv_hill(p::Array{Float64, 1}, y)
    #p = [EC50, min, max, steepness], y:effect. it returns concentration
    conc = p[1] * (((y - p[2]) / (p[3] - y))^p[4])
    return conc
end

function costHill(ydata::Array{Float64, 1}, p::Array{Float64, 1}, conc::Array{Float64, 1})
    y = ydata[1] .+ (ydata[end] - ydata[1]) ./ (1 .+ (p[1] ./ conc) .^ p[2])
    return norm(y - ydata)
end

function optimizeHill(concs::Array{Float64, 2}, d1ind::Int, g1s::Array{Float64, 3}, g2s::Array{Float64, 3})
    nums1 = zeros(9)
    conc1 = zeros(9)
    conc1[1:8] = concs[:, d1ind]
    conc1[9] = 10000
    for i = 1:8
        nums1[i] = g1s[end, i, d1ind] + g2s[end, i, d1ind]
    end
    costs(p) = costHill(nums1, p, conc1)
    low = [conc1[2], 0.1]
    high = [conc1[7], 10.0]
    results_hill = bboptimize(
        costs;
        SearchRange = collect(zip(low, high)),
        NumDimensions = length(low),
        TraceMode = :silent,
        TraceInterval = 100,
        MaxSteps = 1E5,
    )
    par = best_candidate(results_hill)
    return [par[1], nums1[1], nums1[end], par[2]]
end

function low(d1, d2, p1, p2)
    f(x) = (d1 / inv_hill(p1, x)) + (d2 / inv_hill(p2, x)) - 1.0
    find_min = maximum([minimum([p2[2], p2[3]]), minimum([p1[2], p1[3]])])
    find_max = minimum([maximum([p2[2], p2[3]]), maximum([p1[2], p1[3]])])
    combined_effect = find_zero(f, [find_min, find_max])
    return combined_effect
end

function paramsAtEC50(p)
    ps = zeros(9, 5) # num_parameters x number of drugs.
    k = 1
    for i = 1:5
        ps[:, i] = [
            0.5 * (p[36] + p[k + 2]),
            0.5 * (p[37] + p[k + 3]),
            0.5 * p[k + 4],
            0.5 * p[k + 5],
            p[k + 6],
            floor(p[38]),
            floor(p[39]),
            floor(p[40]),
            floor(p[41]),
        ]
        k += 7
    end
    return ps
end

function loweCellNum(concs, d1ind, d2ind, g1s, g2s)
    pars1 = optimizeHill(concs, d1ind, g1s, g2s)
    pars2 = optimizeHill(concs, d2ind, g1s, g2s)
    combined_effs = zeros(9, 9)
    conc1 = zeros(9)
    conc2 = zeros(9)
    conc1[1:8] = concs[:, d1ind]
    conc1[9] = conc2[9] = 10000.0
    conc2[1:8] = concs[:, d2ind]
    for i = 1:9
        for j = 1:9
            combined_effs[i, j] = low(conc1[i], conc2[j], pars1, pars2)
        end
    end
    return combined_effs[1:8, 1:8]
end

function heatmap_combination(d1, d2, cellNum, i1, i2, d1name, d2name, effs, concs, g0)
    n = 8 # the number of concentrations we have
    combin = fullCombinationParam(d1, d2, effs, n)

    numscomb = zeros(n, n)
    for j = 1:n
        for m = 1:n
            numscomb[j, m] = numcells(combin[:, j, m], g0)
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
        clim = (0.0, 5.0),
    )
end
