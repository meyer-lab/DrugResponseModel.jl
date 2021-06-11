""" This file calculates loewe additivity """

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

function optimizeHill(concs::Array{Float64, 2}, d1ind::Int, Total)
    nums1 = zeros(8)
    conc1 = zeros(8)
    conc1[1:8] = concs[:, d1ind]
    for i = 1:8
        nums1[i] = Total[end, i, d1ind]
    end
    costs(p) = costHill(nums1, p, conc1)
    low = [conc1[2], 0.1]
    high = [conc1[8], 10.0]
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

TheHill(p, c) = p[2] + (p[3] - p[2]) / (1 + ((p[1] / c)^p[4]))

function loweCellNum(concs, d1ind, d2ind, Total)
    pars1 = optimizeHill(concs, d1ind, Total)
    pars2 = optimizeHill(concs, d2ind, Total)
    combined_effs = zeros(8, 8)
    conc1 = concs[:, d1ind]
    conc2 = concs[:, d2ind]
    for i = 1:8
        for j = 1:8
            combined_effs[i, j] = low(conc1[i], conc2[j], pars1, pars2)
        end
    end
    return combined_effs
end

function output_loewe()
    concs, _, g1s1, g2s1 = load(189, 1)
    _, _, g1s2, g2s2 = load(189, 2)
    _, _, g1s3, g2s3 = load(189, 3)
    p = parameters()
    gem17 = DrugResponseModel.find_gem17(p)
    efc = getODEparams(p, concs)
    t = LinRange(0.0, 95.0, 189)

    g1m = (g1s1 .+ g1s2 .+ g1s3) ./ 3
    g2m = (g2s1 .+ g2s2 .+ g2s3) ./ 3
    g1m[:, 7, 3], g2m[:, 7, 3], _ = predict(gem17, efc[:, 1, 3], t)
    g1m[:, 8, 3] .= g1m[:, 7, 3]
    g2m[:, 8, 3] .= g2m[:, 7, 3]
    Total = g1m .+ g2m
    lapatinib_palbo = loweCellNum(concs, 1, 5, Total)
    lapatinib_gemc = loweCellNum(concs, 1, 3, Total)
    gemc_palbo = loweCellNum(concs, 3, 5, Total)
    dox_gem = loweCellNum(concs, 2, 3, Total)

    cell = zeros(8, 5)
    for j = 1:5 # drug index
        par = optimizeHill(concs, j, Total)
        for i = 1:8 # conc index
            cell[i, j] = TheHill(par, concs[i, j])
        end
    end
    p1 = plot(concs[:, 1], cell[:, 1], label = "model", title = "lapatinib")
    plot!(concs[:, 1], Total[end, :, 1], label = "data")
    p2 = plot(concs[:, 2], cell[:, 2], label = "model", title = "doxorubicin")
    plot!(concs[:, 2], Total[end, :, 2], label = "data")
    p3 = plot(concs[:, 3], cell[:, 3], label = "model", title = "gemcitabine")
    plot!(concs[:, 3], Total[end, :, 3], label = "data")
    p4 = plot(concs[:, 4], cell[:, 4], label = "model", title = "paclitaxel")
    plot!(concs[:, 4], Total[end, :, 4], label = "data")
    p5 = plot(concs[:, 5], cell[:, 5], label = "model", title = "palbociclib")
    plot!(concs[:, 5], Total[end, :, 5], label = "data")
    p6 = plot(legend = false, grid = false, foreground_color_subplot = :white, top_margin = 1.5cm)
    p = plot(p1, p2, p3, p4, p5, p6, layout = (2, 3), size = (900, 350))
    ylims!((0.0, 4.0))
    savefig(p, "loewe.svg")

    df1 = DataFrames.DataFrame(
        plb50_lpt25 = lapatinib_palbo[4, 5],
        plb50_lpt50 = lapatinib_palbo[5, 5],
        plb50_lpt100 = lapatinib_palbo[6, 5],
        plb50_lpt250 = lapatinib_palbo[7, 5],
    )
    df2 = DataFrames.DataFrame(
        gem10_lpt25 = lapatinib_gemc[4, 6],
        gem10_lpt50 = lapatinib_gemc[5, 6],
        gem10_lpt100 = lapatinib_gemc[6, 6],
        gem10_lpt250 = lapatinib_gemc[7, 6],
    )
    df3 = DataFrames.DataFrame(dox20_gem5 = dox_gem[2, 5], dox20_gem10 = dox_gem[2, 6], dox20_gem17 = dox_gem[2, 7], dox20_gem30 = dox_gem[2, 8])
    df4 = DataFrames.DataFrame(
        plb50_gem5 = gemc_palbo[5, 3],
        plb50_gem10 = gemc_palbo[6, 3],
        plb50_gem17 = gemc_palbo[7, 3],
        plb50_gem30 = gemc_palbo[8, 3],
    )
    df5 = DataFrames.DataFrame(
        lpt100_palbo25 = lapatinib_palbo[6, 4],
        lpt100_palbo50 = lapatinib_palbo[6, 5],
        lpt100_palbo100 = lapatinib_palbo[6, 6],
        lpt100_palbo250 = lapatinib_palbo[6, 7],
    )
    df6 = DataFrames.DataFrame(
        lpt100_gem5 = lapatinib_gemc[6, 5],
        lpt100_gem10 = lapatinib_gemc[6, 6],
        lpt100_gem17 = lapatinib_gemc[6, 7],
        lpt100_gem30 = lapatinib_gemc[6, 8],
    )
    df7 = DataFrames.DataFrame(
        gem10_palbo25 = gemc_palbo[6, 4],
        gem10_palbo50 = gemc_palbo[6, 5],
        gem10_palbo100 = gemc_palbo[6, 6],
        gem10_palbo250 = gemc_palbo[6, 7],
    )
    XLSX.writetable("Loewe_palbo50_lpt.xlsx", df1, overwrite = true, sheetname = "cell number")
    XLSX.writetable("Loewe_gem10_lpt.xlsx", df2, overwrite = true, sheetname = "cell number")
    XLSX.writetable("Loewe_dox20_gem.xlsx", df3, overwrite = true, sheetname = "cell number")
    XLSX.writetable("Loewe_palbo50_gem.xlsx", df4, overwrite = true, sheetname = "cell number")
    XLSX.writetable("Loewe_lapt100_palbo.xlsx", df5, overwrite = true, sheetname = "cell number")
    XLSX.writetable("Loewe_lpt100_gem.xlsx", df6, overwrite = true, sheetname = "cell number")
    XLSX.writetable("Loewe_gem10_palb.xlsx", df7, overwrite = true, sheetname = "cell number")
end
