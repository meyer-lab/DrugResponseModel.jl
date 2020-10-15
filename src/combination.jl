""" This file contains all the functions related to Bliss combinations. """

######---------------- Functions for Bliss combination ----------------########
function CombinationParam(p1::Array{Float64, 2}, p2::Array{Float64, 2}, n::Int)
    """ A function to calculate Bliss independence for drug combination assuming
    the two drugs hit different pathways and they effect independently. """

    pprime1 = copy(p1)
    pprime2 = copy(p2)
    param1 = zeros(4, 8)
    param2 = zeros(4, 8)
    # drug A
    param1[1:2, :] .= 1.0 .- (pprime1[1:2, :] ./ pprime1[1:2, 1]) # g1 and g2 prog. rates
    param1[3:4, :] .= pprime1[3:4, :]                             # g1 and g2 death rates
    # drug B
    param2[1:2, :] .= 1.0 .- (pprime2[1:2, :] ./ pprime2[1:2, 1])
    param2[3:4, :] .= pprime2[3:4, :]

    # make sure normalized correctly; the prog. rates in control condition must be zero in both drugs.
    @assert param1[1,1] == param1[2,1] == param2[1,1] == param2[2,1] == 0.0

    """ For 8x8 combination of drug concentrations for G1 progression rate, G2 progression rate, and death rates in G1 and G2, respectively. """
    combined = zeros(n, n, 4)
    for j = 1:n
        for k = 1:n # param2 is changing, param1 is constant
            combined[j, k, 1:2] .= (1.0 .- (param1[1:2, j] .+ param2[1:2, k] .- param1[1:2, j] .* param2[1:2, k])) .* ((p1[1:2, 1] .+ p2[1:2, 1]) ./ 2)
            combined[j, k, 3:4] .= param1[3:4, j] .+ param2[3:4, k]
        end
    end
    @assert(all(combined .>= 0.0))
    @assert(all(combined .<= 5.0))

    return combined
end


""" To output the full ODE params for plotting the cell number. """
function fullCombinationParam(origP1::Array{Float64, 2}, origP2::Array{Float64, 2}, origFullParam::Array{Float64, 3}, n::Int)
    """ Here we assume the base is origP1, and we just want to get the params of EC50 from origP2. """
    combined = CombinationParam(origP1, origP2, n)
    fullparam = zeros(9, n, n)
    fullparam[5:9, :, :] .= origFullParam[5:9, 1, 1]
    for i = 1:4
        fullparam[i, :, :] .= combined[:, :, i]
    end
    return fullparam
end

""" This function calculates cell number for parameter sets that are the result of Bliss on prog. rates. """
function BlissModelComb(bliss_comb, g0)
      bliss_comb_cellnum = zeros(8,8)
      for i=1:8 # param2 is changing
            for j=1:8 # param1 is changing, param 2 is constant
                  g1, g2, _ = predict(bliss_comb[:, j, i], g0, 189)
                  bliss_comb_cellnum[j, i] = g1 + g2
            end
      end
      return bliss_comb_cellnum
end

""" This function plots the heatmap of combined cell numbers given any two drugs. """
function Heatmap(concs, data, i1, i2, d1name, d2name, title)
    concs[1, :] .= 0.05
    heatmap(
        string.(round.(log.(concs[:, i1]), digits = 1)),
        string.(round.(log.(concs[:, i2]), digits = 1)),
        data,
        xlabel = string(d1name, " log[nM]"),
        ylabel = string(d2name, " log [nM]"),
        title = title,
        clim = (0.0, 35.0),
    )
end


""" In this function, we apply the Bliss synergy to the cell numbers. 
This is to compare the traditional way of representing the combination effect, compare to the way we do in our model."""
function blissCellNum(g1s, g2s; T = 189, n = 8)
    num = zeros(n, 5)

    for i = 1:5
        # num is a 8 x 5 matrix, holding scaled cell numbers for 5 drugs, in 8 concenntration, for a specific time point.
        num[:, i] = 1.0 .- ((g1s[T, :, i] .+ g2s[T, :, i]) ./ (g1s[T, 1, i] + g2s[T, 1, i]))
    end
    combined = zeros(n, n, 10)
    for j = 1:n
        # the base case for either of combinations is the average of drugA and drugB to scale the cell numbers back.
        # columns are drugA and rows are drug 2, and the third dimension shows which pair of drugs were combined.
        combined[j, :, 1] = -(num[:, 1] .+ num[j, 2] .- (num[:, 1] .* num[j, 2]) .- 1.0) .* ((g1s[T, 1, 1] + g2s[T, 1, 1]) + (g1s[T, 1, 2] + g2s[T, 1, 2])) / 2 # lap w/ dox; meaning dox changes with rows and lap changes with columns
        combined[j, :, 2] = -(num[:, 1] .+ num[j, 3] .- (num[:, 1] .* num[j, 3]) .- 1.0) .* ((g1s[T, 1, 1] + g2s[T, 1, 1]) + (g1s[T, 1, 3] + g2s[T, 1, 3])) / 2 # lap w/ gem
        combined[j, :, 3] = -(num[:, 1] .+ num[j, 4] .- (num[:, 1] .* num[j, 4]) .- 1.0) .* ((g1s[T, 1, 1] + g2s[T, 1, 1]) + (g1s[T, 1, 4] + g2s[T, 1, 4])) / 2 # lap w/ pac
        combined[j, :, 4] = -(num[:, 1] .+ num[j, 5] .- (num[:, 1] .* num[j, 5]) .- 1.0) .* ((g1s[T, 1, 1] + g2s[T, 1, 1]) + (g1s[T, 1, 5] + g2s[T, 1, 5])) / 2# lap w/ palb
        combined[j, :, 5] = -(num[:, 2] .+ num[j, 3] .- (num[:, 2] .* num[j, 3]) .- 1.0) .* ((g1s[T, 1, 2] + g2s[T, 1, 2]) + (g1s[T, 1, 3] + g2s[T, 1, 3])) / 2 # dox w/ gem
        combined[j, :, 6] = -(num[:, 2] .+ num[j, 4] .- (num[:, 2] .* num[j, 4]) .- 1.0) .* ((g1s[T, 1, 2] + g2s[T, 1, 2]) + (g1s[T, 1, 4] + g2s[T, 1, 4])) / 2 # dox w/ pac
        combined[j, :, 7] = -(num[:, 2] .+ num[j, 5] .- (num[:, 2] .* num[j, 5]) .- 1.0) .* ((g1s[T, 1, 2] + g2s[T, 1, 2]) + (g1s[T, 1, 5] + g2s[T, 1, 5])) / 2 # dox w/ palb
        combined[j, :, 8] = -(num[:, 3] .+ num[j, 4] .- (num[:, 3] .* num[j, 4]) .- 1.0) .* ((g1s[T, 1, 3] + g2s[T, 1, 3]) + (g1s[T, 1, 4] + g2s[T, 1, 4])) / 2 # gem w/ pac
        combined[j, :, 9] = -(num[:, 3] .+ num[j, 5] .- (num[:, 3] .* num[j, 5]) .- 1.0) .* ((g1s[T, 1, 3] + g2s[T, 1, 3]) + (g1s[T, 1, 5] + g2s[T, 1, 5])) / 2 # gem w/ palb
        combined[j, :, 10] = -(num[:, 4] .+ num[j, 5] .- (num[:, 4] .* num[j, 5]) .- 1.0) .* ((g1s[T, 1, 4] + g2s[T, 1, 4]) + (g1s[T, 1, 5] + g2s[T, 1, 5])) / 2 # pac w/ palb
    end

    @assert(all(combined .>= 0.0))
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
