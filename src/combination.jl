""" This file contains all the functions related to Bliss combinations. """

######---------------- Functions for Bliss combination ----------------########


""" Unit function to calculate the bliss for 2 drugs at one specific concentration. """
function Bliss_params_unit(pp1, pp2, control)
    # pp1 and pp2 are 1D arrays of size 8, including 8 parameters for a single concentration.
    p1 = copy(pp1)
    p2 = copy(pp2)
    # normalization
    p1[1:8] .= 1.0 .- (pp1[1:8] ./ control[1:8, 1]) # g1 and g2 prog. rates

    # drug B
    p2[1:8] .= 1.0 .- (pp2[1:8] ./ control[1:8, 2])

    c = Array{eltype(pp1), 1}(undef, 16)
    c[1:8] .= (1.0 .- (p1[1:8] .+ p2[1:8] .- p1[1:8] .* p2[1:8])) .* ((control[1:8, 1] .+ control[1:8, 2]) ./ 2)
    c[9:16] .= pp1[9:16] .+ pp2[9:16]

    c
end

""" Using the unit function to find all combinations of parameters. """
function AllBliss_params(pp1, pp2)
    # pp1 and pp2 are 2D arrays [16 x 8] each includes the parameters fo all concentrations of a drug. 

    n = length(pp1[1, :])
    combined = Array{eltype(pp1), 3}(undef, 16, n, n)
    for i = 1:n
        for j = 1:n
            combined[:, i, j] .= Bliss_params_unit(pp1[:, i], pp2[:, j], hcat(pp1[:, 1], pp2[:, 1]))
        end
    end
    @assert all(combined[9:end, 1, 1] .== 0.0)
    @assert(all(combined .>= 0.0))
    combined
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
function blissCellNum(gs)
    nconc = length(gs[:, 1])
    ndrugs = length(gs[1, :])

    num = zeros(nconc, ndrugs)

    for i = 1:ndrugs
        # num is a 8 x 5 matrix, holding scaled cell numbers for 5 drugs, in 8 concentrations, for a specific time point.
        num[:, i] = 1.0 .- (gs[:, i] ./ gs[1, i])
    end

    ncombs = binomial(ndrugs, 2)
    combined = zeros(nconc, nconc, ncombs)
    x = 1
    for i = 1:ndrugs - 1
        for k = (i + 1):ndrugs
            # the base case for either of combinations is the average of drugA and drugB to scale the cell numbers back.
            # columns are drugA and rows are drug 2, and the third dimension shows which pair of drugs were combined.
            for j = 1:nconc
                combined[:, j, x] = -(num[:, i] .+ num[j, k] .- (num[:, i] .* num[j, k]) .- 1.0) .* (gs[1, i] + gs[1, k]) / 2
            end
            x += 1
        end
    end

    @assert all(combined .>= 0.0)
    return combined
end

""" only to find the combination of one concentration of drug 1, and one concentration of drug 2, over time. """
function pair_cellnum_Bliss(total1, total2)
    # note that each of the two inputs, should be hcat with control: total1: [control; condition_i]
    normedtotal1 = 1.0 .- (total1[:, 2] ./ total1[:, 1]) # normalize to control
    normedtotal2 = 1.0 .- (total2[:, 2] ./ total2[:, 1]) # normalize to control
    combined = -(normedtotal1 .+ normedtotal2 .- (normedtotal1 .* normedtotal2) .- 1.0) .* (total1[:, 1] .+ total2[:, 1]) / 2
    return combined
end

"""This function outputs the drug combination from cell counts over time."""
function output_Bliss_cellnum()
    # data import
    gt1, gt2 = DrugResponseModel.import_combination("AU01001")
    gt1_2, gt2_2 = DrugResponseModel.import_combination("AU01101")
    gt1_3, gt2_3 = DrugResponseModel.import_combination("AU00901")
    concs, _, _, _ = load(189, 1)
    GS1 = cat(gt1, gt1_2, gt1_3, dims = 4)
    GS2 = cat(gt2, gt2_2, gt2_3, dims = 4)
    meanGS1 = mean(GS1, dims = 4)
    meanGS2 = mean(GS2, dims = 4)
    meanGS2[:, :, 19] .= mean(cat(gt2[:, :, 19], gt2_2[:, :, 19], dims = 3), dims = 3)[:, :, 1]
    Total = zeros(193, 5, 5) # time x concentrations x 5 drugs
    Total[:, 1, :] .= meanGS1[3, :, 1] # controls
    Total[:, 2:5, 1] .= meanGS1[3, :, 2:5] # lapatinibs
    Total[:, 2, 2] .= meanGS1[3, :, 6] # dox 20 nM
    Total[:, 2:5, 3] .= meanGS1[3, :, 19:22] # gemcitabines
    Total[:, 2, 4] .= meanGS1[3, :, 13] # pax 2 nM
    Total[:, 2:5, 5] .= meanGS1[3, :, 7:10] # palbos
    bliss = zeros(189, 5, 5, 10) # the first one changes with rows, which is the drug that comes first (e.g., in lapatinib (rows-first) dox (columns-second))
    for i = 1:189
        bliss[i, :, :, :] .= blissCellNum(Total[i, :, :])
    end
    bliss
end