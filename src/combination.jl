""" This file contains all the functions related to Bliss combinations. """

######---------------- Functions for Bliss combination ----------------########
function CombinationParam(p1, p2, n::Int)
    """ A function to calculate Bliss independence for drug combination assuming
    the two drugs hit different pathways and they effect independently. """

    pprime1 = copy(p1)
    pprime2 = copy(p2)
    param1 = Matrix{eltype(p1)}(undef, 4, 8)
    param2 = Matrix{eltype(p2)}(undef, 4, 8)
    # drug A
    param1[1:2, :] .= 1.0 .- (pprime1[1:2, :] ./ pprime1[1:2, 1]) # g1 and g2 prog. rates
    param1[3:4, :] .= pprime1[3:4, :]                             # g1 and g2 death rates
    # drug B
    param2[1:2, :] .= 1.0 .- (pprime2[1:2, :] ./ pprime2[1:2, 1])
    param2[3:4, :] .= pprime2[3:4, :]

    # make sure normalized correctly; the prog. rates in control condition must be zero in both drugs.
    @assert param1[1, 1] == param1[2, 1] == param2[1, 1] == param2[2, 1] == 0.0

    # For 8x8 combination of drug concentrations for G1 progression rate, G2 progression rate, and death rates in G1 and G2, respectively.
    combined = Matrix{eltype(p1)}(undef, n, n, 4)
    for j = 1:n
        for k = 1:n # param2 is changing, param1 is constant
            combined[j, k, 1:2] .=
                (1.0 .- (param1[1:2, j] .+ param2[1:2, k] .- param1[1:2, j] .* param2[1:2, k])) .* ((p1[1:2, 1] .+ p2[1:2, 1]) ./ 2)
            combined[j, k, 3:4] .= param1[3:4, j] .+ param2[3:4, k]
        end
    end
    @assert(all(combined .>= 0.0))
    @assert(all(combined .<= 5.0))

    return combined
end


""" To output the full ODE params for plotting the cell number. """
function fullCombinationParam(origP1, origP2, origFullParam, n::Int)
    """ Here we assume the base is origP1. """
    combined = CombinationParam(origP1, origP2, n)
    fullparam = zeros(9, n, n)
    fullparam[5:9, :, :] .= origFullParam[5:9, 1, 1]
    fullparam[1:4, :, :] .= permutedims(combined[:, :, 1:4], (3, 1, 2))
    return fullparam
end

""" This function calculates cell number for parameter sets that are the result of Bliss on prog. rates. """
function BlissModelComb(bliss_comb, g0)
    bliss_comb_cellnum = zeros(8, 8)
    for i = 1:8 # param1 is changing
        for j = 1:8 # param2 is changing
            g1, g2, _ = predict(bliss_comb[:, i, j], g0, 96.0)
            bliss_comb_cellnum[i, j] = g1 + g2
        end
    end
    return bliss_comb_cellnum
end

""" This function plots the heatmap of combined cell numbers given any two drugs. """
function Heatmap(concs, data, i1, i2, d1name, d2name, title; clim_min = 0.0, clim_max = 60.0)
    concs[1, :] .= 0.05
    heatmap(
        string.(round.(log.(concs[:, i2]), digits = 1)),
        string.(round.(log.(concs[:, i1]), digits = 1)),
        data,
        xlabel = string(d2name, " log[nM]"),
        ylabel = string(d1name, " log[nM]"),
        title = title,
        clim = (clim_min, clim_max),
    )
end


""" In this function, we apply the Bliss synergy to the cell numbers. 
This is to compare the traditional way of representing the combination effect, compare to the way we do in our model."""
function blissCellNum(g1s, g2s; n = 8)
    gs = g1s[end, :, :] + g2s[end, :, :]
    num = zeros(n, 5)

    for i = 1:5
        # num is a 8 x 5 matrix, holding scaled cell numbers for 5 drugs, in 8 concentrations, for a specific time point.
        num[:, i] = 1.0 .- (gs[:, i] ./ gs[1, i])
    end

    combined = zeros(n, n, 10)
    x = 1
    for i = 1:4
        for k = (i + 1):5
            # the base case for either of combinations is the average of drugA and drugB to scale the cell numbers back.
            # columns are drugA and rows are drug 2, and the third dimension shows which pair of drugs were combined.
            for j = 1:n
                combined[:, j, x] = -(num[:, i] .+ num[j, k] .- (num[:, i] .* num[j, k]) .- 1.0) .* (gs[1, i] + gs[1, k]) / 2
            end
            x += 1
        end
    end

    @assert all(combined .>= 0.0)
    return combined
end

""" This function calculates parameter values at their EC50 for all drugs given estimated parameters from AllDrugAtOnce fitting. """
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
