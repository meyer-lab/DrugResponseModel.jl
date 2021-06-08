""" This file contains all the functions related to Bliss combinations. """

######---------------- Functions for Bliss combination ----------------########

""" Unit function to calculate the bliss for 2 drugs at one specific concentration. """
function Bliss_params_unit(pp1, pp2, control)
    # pp1 and pp2 are 1D arrays of size 8, including 8 parameters for a single concentration.
    p1 = copy(pp1)
    p2 = copy(pp2)
    # normalization
    p1[1:6] .= 1.0 .- (pp1[1:6] ./ control[1:6, 1]) # g1 and g2 prog. rates

    # drug B
    p2[1:6] .= 1.0 .- (pp2[1:6] ./ control[1:6, 2])

    c = Array{eltype(pp1), 1}(undef, 12)
    c[1:6] .= (1.0 .- (p1[1:6] .+ p2[1:6] .- p1[1:6] .* p2[1:6])) .* ((control[1:6, 1] .+ control[1:6, 2]) ./ 2)
    c[7:12] .= pp1[7:12] .+ pp2[7:12]

    c
end

""" Using the unit function to find all combinations of parameters. """
function AllBliss_params(pp1, pp2; n = 8)
    # pp1 and pp2 are 2D arrays [12 x 8] each includes the parameters fo all concentrations of a drug. 

    combined = Array{eltype(pp1), 3}(undef, 12, n, n)
    for i = 1:n
        for j = 1:n
            combined[:, i, j] .= Bliss_params_unit(pp1[:, i], pp2[:, j], hcat(pp1[:, 1], pp2[:, 1]))
        end
    end
    @assert all(combined[7:end, 1, 1] .== 0.0)
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
    concs, _, g1s1, g2s1 = load(189, 1);
    _, _, g1s2, g2s2 = load(189, 2);
    _, _, g1s3, g2s3 = load(189, 3);
    p = parameters()
    gem17 = DrugResponseModel.find_gem17(p)
    efc = getODEparams(p, concs)
    t = LinRange(0.0, 95.0, 189)

    g1m = (g1s1 .+ g1s2 .+ g1s3) ./ 3;
    g2m = (g2s1 .+ g2s2 .+ g2s3) ./ 3;
    g1m[:, 7, 3], g2m[:, 7, 3], _ = predict(gem17, efc[:, 1, 3], t)
    g1m[:, 8, 3] .= g1m[:, 7, 3]
    g2m[:, 8, 3] .= g2m[:, 7, 3]
    Total = g1m .+ g2m
    bliss = blissCellNum(Total[end, :, :]; n = 8)

    df1 = DataFrames.DataFrame(plb50_lpt25=bliss[4, 5, 4], plb50_lpt50=bliss[5, 5, 4], plb50_lpt100=bliss[6, 5, 4], plb50_lpt250=bliss[7, 5, 4])
    df2 = DataFrames.DataFrame(gem10_lpt25=bliss[4, 6, 2], gem10_lpt50=bliss[5, 6, 2], gem10_lpt100=bliss[6, 6, 2], gem10_lpt250=bliss[7, 6, 2])
    df3 = DataFrames.DataFrame(dox20_gem5=bliss[2, 5, 5], dox20_gem10=bliss[2, 6, 5], dox20_gem17=bliss[2, 7, 5], dox20_gem30=bliss[2, 8, 5])
    df4 = DataFrames.DataFrame(plb50_gem5=bliss[5, 3, 9], plb50_gem10=bliss[6, 3, 9], plb50_gem17=bliss[7, 3, 9], plb50_gem30=bliss[8, 3, 9])
    df5 = DataFrames.DataFrame(lpt100_palbo25=bliss[6, 4, 4], lpt100_palbo50=bliss[6, 5, 4], lpt100_palbo100=bliss[6, 6, 4], lpt100_palbo250=bliss[6, 7, 4])
    df6 = DataFrames.DataFrame(lpt100_gem5=bliss[6, 5, 2], lpt100_gem10=bliss[6, 6, 2], lpt100_gem17=bliss[6, 7, 2], lpt100_gem30=bliss[6, 8, 2])
    df7 = DataFrames.DataFrame(gem10_palbo25=bliss[6, 4, 9], gem10_palbo50=bliss[6, 5, 9], gem10_palbo100=bliss[6, 6, 9], gem10_palbo250=bliss[6, 7, 9])
    XLSX.writetable("Bliss_palbo50_lpt.xlsx", df1, overwrite=true, sheetname="cell number")
    XLSX.writetable("Bliss_gem10_lpt.xlsx", df2, overwrite=true, sheetname="cell number")
    XLSX.writetable("Bliss_dox20_gem.xlsx", df3, overwrite=true, sheetname="cell number")
    XLSX.writetable("Bliss_palbo50_gem.xlsx", df4, overwrite=true, sheetname="cell number")
    XLSX.writetable("Bliss_lpt100_palbo.xlsx", df5, overwrite=true, sheetname="cell number")
    XLSX.writetable("Bliss_lpt100_gem.xlsx", df6, overwrite=true, sheetname="cell number")
    XLSX.writetable("Bliss_gem10_palbo.xlsx", df7, overwrite=true, sheetname="cell number")

end


