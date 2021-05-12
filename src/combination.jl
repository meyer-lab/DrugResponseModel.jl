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
function blissCellNum(gs; n = 8)
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

""" only to find the combination of one concentration of drug 1, and once concentration of drug 2, over time. """
function pair_cellnum_Bliss(total1, total2)
    # note that each of the two inputs, should be hcat with control: total1: [control; condition_i]
    normedtotal1 = 1.0 .- (total1[:, 2] ./ total1[:, 1]) # normalize to control
    normedtotal2 = 1.0 .- (total2[:, 2] ./ total2[:, 1]) # normalize to control
    combined = -(normedtotal1 .+ normedtotal2 .- (normedtotal1 .* normedtotal2) .- 1.0) .* (total1[:, 1] .+ total2[:, 1]) / 2
    return combined
end

function output_Bliss_cellnum()
    # data import
    gt1, gt2 = DrugResponseModel.import_combination("AU01001");
    gt1_2, gt2_2 = DrugResponseModel.import_combination("AU01101");
    gt1_3, gt2_3 = DrugResponseModel.import_combination("AU00901");
    concs, _, _, _ = load(189, 1);

    GS1 = cat(gt1, gt1_2, gt1_3, dims = 4);
    GS2 = cat(gt2, gt2_2, gt2_3, dims = 4);

    meanGS1 = mean(GS1, dims = 4);
    meanGS2 = mean(GS2, dims = 4);
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
        bliss[i, :, :, :] .= blissCellNum(Total[i, :, :]; n = 5)
    end
    bliss
end

# df1 = DataFrames.DataFrame(plb50_lpt25=bliss[:, 2, 3, 4], plb50_lpt50=bliss[:, 3, 3, 4], plb50_lpt100=bliss[:, 4, 3, 4], plb50_lpt250=bliss[:, 5, 3, 4])
# df2 = DataFrames.DataFrame(gem10_lpt25=bliss[:, 2, 3, 2], gem10_lpt50=bliss[:, 3, 3, 2], gem10_lpt100=bliss[:, 4, 3, 2], gem10_lpt250=bliss[:, 5, 3, 2])
# df3 = DataFrames.DataFrame(dox20_gem5=bliss[:, 1, 2, 5], dox20_gem10=bliss[:, 1, 3, 5], dox20_gem17=bliss[:, 1, 4, 5], dox20_gem30=bliss[:, 1, 5, 5])
# df4 = DataFrames.DataFrame(plb50_gem5=bliss[:, 2, 3, 9], plb50_gem10=bliss[:, 3, 3, 9], plb50_gem17=bliss[:, 4, 3, 9], plb50_gem30=bliss[:, 5, 3, 9])
# df5 = DataFrames.DataFrame(lpt100_palbo25=bliss[:, 4, 2, 4], lpt100_palbo50=bliss[:, 4, 3, 4], lpt100_palbo100=bliss[:, 4, 4, 4], lpt100_palbo250=bliss[:, 4, 5, 4])
# df6 = DataFrames.DataFrame(lpt100_gem5=bliss[:, 4, 2, 2], lpt100_gem10=bliss[:, 4, 3, 2], lpt100_gem17=bliss[:, 4, 4, 2], lpt100_gem30=bliss[:, 4, 5, 2])
# df7 = DataFrames.DataFrame(gem10_palbo25=bliss[:, 3, 2, 9], gem10_palbo50=bliss[:, 3, 3, 9], gem10_palbo100=bliss[:, 3, 4, 9], gem10_palbo250=bliss[:, 3, 5, 9])
# df8 = DataFrames.DataFrame(Pax2_lpt50 = bliss[:, 3, 1, 3], pax2_lpt100=bliss[:, 4, 1, 3], pax2_dox20=bliss[:, 1, 1, 6], pax2_gem10=bliss[:, 3, 1, 8], pax2_palbo50=bliss[:, 1, 3, 10], dox20_lpt100=bliss[:, 4, 1, 1], dox20_palbo50=bliss[:, 1, 3, 7])
# XLSX.writetable("Bliss_model.xlsx", df1, overwrite=true, sheetname="cell number")