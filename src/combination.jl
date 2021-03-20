""" This file contains all the functions related to Bliss combinations. """

######---------------- Functions for Bliss combination ----------------########

""" Unit function to calculate the bliss for 2 drugs at one specific concentration. """
function Bliss_params_unit(pp1, pp2, control)
    # pp1 and pp2 are 1D arrays of size 8, including 8 parameters for a single concentration.
    p1 = copy(pp1)
    p2 = copy(pp2)
    # normalization
    p1[1:6] .= 1.0 .- (pp1[1:6] ./ control[1:6, 1]) # g1 and g2 prog. rates
    p1[7:12] .= pp1[7:12]                          # g1 and g2 death rates
    # drug B
    p2[1:6] .= 1.0 .- (pp2[1:6] ./ control[1:6, 2])
    p2[7:12] .= pp2[7:12]

    c = Array{eltype(pp1), 1}(undef, 12)
    c[1:6] .= (1.0 .- (p1[1:6] .+ p2[1:6] .- p1[1:6] .* p2[1:6])) .* ((control[1:6, 1] .+ control[1:6, 2]) ./ 2)
    c[7:12] .= p1[7:12] .+ p2[7:12]

    c
end

""" Using the unit function to find all combinations of parameters. """
function AllBliss_params(pp1, pp2)
    # pp1 and pp2 are 2D arrays [12 x 8] each includes the parameters fo all concentrations of a drug. 

    combined = Array{eltype(pp1), 3}(undef, 12, 8, 8)
    for i = 1:8
        for j = 1:8
            combined[:, i, j] .= Bliss_params_unit(pp1[:, i], pp2[:, j], hcat(pp1[:, 1], pp2[:, 1]))
        end
    end
    @assert all(combined[7:end, 1, 1] .== 0.0)
    # TODO: remember to uncomment this assertion after estimating correct set of parameters
    @assert(all(combined .>= 0.0))
    combined
end


""" This function calculates cell number for parameter sets that are the result of Bliss on prog. rates. """
function BlissModelComb(bliss_comb, pCtr)
    bliss_comb_cellnum = Matrix{eltype(bliss_comb)}(undef, 8, 8)

    for i = 1:8 # param1 is changing
        for j = 1:8 # param2 is changing
            g1, g2, _ = predict(bliss_comb[:, i, j], pCtr, 96.0)
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
function blissCellNum(g1s, g2s, tt; n = 8)
    gs = g1s[tt, :, :]
    gs2 = g2s[tt, :, :]
    num = zeros(n, 5)
    num2 = zeros(n, 5)

    for i = 1:5
        # num is a 8 x 5 matrix, holding scaled cell numbers for 5 drugs, in 8 concentrations, for a specific time point.
        num[:, i] = 1.0 .- (gs[:, i] ./ gs[1, i])
        num2[:, i] = 1.0 .- (gs2[:, i] ./ gs2[1, i])
    end

    combined = zeros(n, n, 10)
    combined2 = zeros(n, n, 10)
    x = 1
    for i = 1:4
        for k = (i + 1):5
            # the base case for either of combinations is the average of drugA and drugB to scale the cell numbers back.
            # columns are drugA and rows are drug 2, and the third dimension shows which pair of drugs were combined.
            for j = 1:n
                combined[:, j, x] = -(num[:, i] .+ num[j, k] .- (num[:, i] .* num[j, k]) .- 1.0) .* (gs[1, i] + gs[1, k]) / 2
                combined2[:, j, x] = -(num2[:, i] .+ num2[j, k] .- (num2[:, i] .* num2[j, k]) .- 1.0) .* (gs2[1, i] + gs2[1, k]) / 2
            end
            x += 1
        end
    end

    @assert all(combined .>= 0.0)
    return combined, combined2
end

""" only to find the combination of one concentration of drug 1, and once concentration of drug 2, over time. """
function pair_cellnum_Bliss(total1, total2)
    # note that each of the two inputs, should be vcat with control: total1: [control; condition_i]
    normedtotal1 = 1.0 .- (total1[:, 2] ./ total1[:, 1]) # normalize to control
    normedtotal2 = 1.0 .- (total2[:, 2] ./ total2[:, 1]) # normalize to control
    combined = -(normedtotal1 .+ normedtotal2 .- (normedtotal1 .* normedtotal2) .- 1.0) .* (total1[:, 1] .+ total2[:, 1])/2
    return combined
end