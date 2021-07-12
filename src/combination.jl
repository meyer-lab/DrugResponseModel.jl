""" This file contains all the functions related to Bliss combinations. """

######---------------- Functions for Bliss combination ----------------########

function residHillAll3(hP, concentrations::Matrix, g1::Array, g2::Array)
    res = 0.0

    # Solve for all drugs
    t = 1
    for j = 1:3
        hill = hP[[t:(t + 17); 55:62]]
        for i = 3:10
            res += 20 * (maximum([0, (hill[i] - hill[i + 16])]))^2
        end
        res += residHill(hill, concentrations[:, j], g1[:, :, j], g2[:, :, j])
        t += 18
    end
    return res
end

function optim_all3(concs::Array{Float64, 2}, g1::Array{Float64, 3}, g2::Array{Float64, 3}; maxiter = 800000)
    f(x) = residHillAll3(x, concs, g1, g2)

    lP = [minimum(concs); 0.01; 5e-9 * ones(16)]
    low = vcat(lP, lP, lP, 5e-9, 5e-9, 5e-9, 5e-9, 5e-9, 5e-9, 5e-9, 5e-9)
    hP = [maximum(concs); 10.0; 2.0 * ones(16)]
    high = vcat(hP, hP, hP, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0)

    return optimize_helper(f, low, high, maxiter)
end

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
function AllBliss_params(pp1, pp2; n = 8)
    # pp1 and pp2 are 2D arrays [12 x 8] each includes the parameters fo all concentrations of a drug. 

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

    # data import (OLD EXPERIMENT)
    conc, _, g1s1, g2s1 = load(189, 1)
    _, _, g1s2, g2s2 = load(189, 2)
    _, _, g1s3, g2s3 = load(189, 3)
    p = parameters()
    gem17 = find_gem17(p)
    efc = getODEparams(p, conc)
    t = LinRange(0.0, 95.0, 189)

    g1m = (g1s1 .+ g1s2 .+ g1s3) ./ 3
    g2m = (g2s1 .+ g2s2 .+ g2s3) ./ 3
    g1m[:, 8, 3] .= g1m[:, 7, 3]
    g2m[:, 8, 3] .= g2m[:, 7, 3]
    g1m[:, 7, 3], g2m[:, 7, 3], _ = predict(gem17, efc[:, 1, 3], t)
    Total1 = g1m .+ g2m

    bliss = blissCellNum(Total1[end, :, :]; n = 8)
    df1 =
        DataFrames.DataFrame(plb50_lpt25 = bliss[4, 5, 4], plb50_lpt50 = bliss[5, 5, 4], plb50_lpt100 = bliss[6, 5, 4], plb50_lpt250 = bliss[7, 5, 4])
    df2 =
        DataFrames.DataFrame(gem10_lpt25 = bliss[4, 6, 2], gem10_lpt50 = bliss[5, 6, 2], gem10_lpt100 = bliss[6, 6, 2], gem10_lpt250 = bliss[7, 6, 2])
    df3 = DataFrames.DataFrame(dox20_gem5 = bliss[4, 5, 5], dox20_gem10 = bliss[4, 6, 5], dox20_gem17 = bliss[4, 7, 5], dox20_gem30 = bliss[4, 8, 5])
    df4 = DataFrames.DataFrame(plb50_gem5 = bliss[5, 5, 9], plb50_gem10 = bliss[6, 5, 9], plb50_gem17 = bliss[7, 5, 9], plb50_gem30 = bliss[8, 5, 9])
    df5 = DataFrames.DataFrame(
        lpt100_palbo25 = bliss[6, 4, 4],
        lpt100_palbo50 = bliss[6, 5, 4],
        lpt100_palbo100 = bliss[6, 6, 4],
        lpt100_palbo250 = bliss[6, 7, 4],
    )
    df6 = DataFrames.DataFrame(
        lpt100_gem5 = bliss[6, 5, 2],
        lpt100_gem10 = bliss[6, 6, 2],
        lpt100_gem17 = bliss[6, 7, 2],
        lpt100_gem30 = bliss[6, 8, 2],
    )
    df7 = DataFrames.DataFrame(
        gem10_palbo25 = bliss[6, 4, 9],
        gem10_palbo50 = bliss[6, 5, 9],
        gem10_palbo100 = bliss[6, 6, 9],
        gem10_palbo250 = bliss[6, 7, 9],
    )
    XLSX.writetable("Bliss_OLDexpCellNums.xlsx", df1, overwrite = true, sheetname = "cell number")

    # # data import (NEW EXPERIMENT)
    gt1, gt2 = DrugResponseModel.import_combination("AU01001")
    gt1_2, gt2_2 = DrugResponseModel.import_combination("AU01101")
    gt1_3, gt2_3 = DrugResponseModel.import_combination("AU00901")

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
    bliss = blissCellNum(Total[end, :, :]; n = 5)
    df1 = DataFrames.DataFrame(
        plb50_lpt25 = bliss[2, 3, 4],
        plb50_lpt50 = bliss[3, 3, 4],
        plb50_lpt100 = bliss[4, 3, 4],
        plb50_lpt250 = bliss[5, 3, 4],
        gem10_lpt25 = bliss[2, 3, 2],
        gem10_lpt50 = bliss[3, 3, 2],
        gem10_lpt100 = bliss[4, 3, 2],
        gem10_lpt250 = bliss[5, 3, 2],
        dox20_gem5 = bliss[2, 2, 5],
        dox20_gem10 = bliss[2, 3, 5],
        dox20_gem17 = bliss[2, 4, 5],
        dox20_gem30 = bliss[2, 5, 5],
        plb50_gem5 = bliss[2, 3, 9],
        plb50_gem10 = bliss[3, 3, 9],
        plb50_gem17 = bliss[4, 3, 9],
        plb50_gem30 = bliss[5, 3, 9],
        lpt100_palbo25 = bliss[4, 2, 4],
        lpt100_palbo50 = bliss[4, 3, 4],
        lpt100_palbo100 = bliss[4, 4, 4],
        lpt100_palbo250 = bliss[4, 5, 4],
        lpt100_gem5 = bliss[4, 2, 2],
        lpt100_gem10 = bliss[4, 3, 2],
        lpt100_gem17 = bliss[4, 4, 2],
        lpt100_gem30 = bliss[4, 5, 2],
        gem10_palbo25 = bliss[3, 2, 9],
        gem10_palbo50 = bliss[3, 3, 9],
        gem10_palbo100 = bliss[3, 4, 9],
        gem10_palbo250 = bliss[3, 5, 9],
    )
    XLSX.writetable("Bliss_NEWexpCellNums.xlsx", df1, overwrite = true, sheetname = "cell number")

    concs = zeros(5, 5)
    concs[1, :] .= 0.0
    concs[2:5, [1, 5]] .= [25.0, 50.0, 100.0, 250]
    concs[2:5, 2] .= 20.0
    concs[2:5, 3] .= [5.0, 10.0, 17.0, 30.0]
    concs[2:5, 4] .= 2.0

    p1 = plot(concs[:, 1], Total1[end, [1, 4, 5, 6, 7], 1], label = "old", title = "lapatinib", xlabel = "concentration", ylabel = "cell number")
    plot!(concs[:, 1], Total[end, :, 1], label = "new")
    p2 = plot(concs[:, 2], Total1[end, [1, 4, 4, 4, 4], 2], label = "old", title = "doxorubicin", xlabel = "concentration", ylabel = "cell number")
    plot!(concs[:, 2], Total[end, :, 2], label = "new")
    p3 = plot(concs[:, 3], Total1[end, [1, 5, 6, 7, 8], 3], label = "old", title = "gemcitabine", xlabel = "concentration", ylabel = "cell number")
    plot!(concs[:, 3], Total[end, :, 3], label = "new")
    p5 = plot(concs[:, 5], Total1[end, [1, 4, 5, 6, 7], 5], label = "old", title = "palbociclib", xlabel = "concentration", ylabel = "cell number")
    plot!(concs[:, 5], Total[end, :, 5], label = "new")
    p6 = plot(legend = false, grid = false, foreground_color_subplot = :white, top_margin = 1.5cm)
    p = plot(p1, p2, p3, p5, layout = (2, 2), size = (400, 350))
    ylims!((0.0, 4.0))
    savefig(p, "bliss.svg")
end
