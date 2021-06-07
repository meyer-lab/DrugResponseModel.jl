""" This file calculates loewe additivity """

""" Find the inverse of a hill function (concentration), given the parameters and the effect. """
function inv_hill(p::Array{Float64, 1}, y)
    #p = [EC50, min, max, steepness], y:effect. it returns concentration
    conc = p[1] * (((y - p[2]) / (p[3] - y)) ^ p[4])
    return conc
end

function costHill(ydata::Array{Float64, 1}, p::Array{Float64, 1}, conc::Array{Float64, 1})
    y = ydata[1] .+ (ydata[end] - ydata[1]) ./ (1 .+ (p[1] ./ conc) .^ p[2])
    return norm(y - ydata)
end

function optimizeHill(concs::Array{Float64, 2}, d1ind::Int, Total)
    nums1 = zeros(6)
    conc1 = zeros(6)
    conc1[1:5] = concs[:, d1ind]
    conc1[6] = 10000
    for i=1:5
        nums1[i] = Total[end, i, d1ind]
    end
    costs(p) = costHill(nums1, p, conc1)
    low = [conc1[2], 0.1]
    high = [conc1[5], 10.0]
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

function loweCellNum(concs, d1ind, d2ind, Total)
    pars1 = optimizeHill(concs, d1ind, Total)
    pars2 = optimizeHill(concs, d2ind, Total)
    combined_effs = zeros(6, 6)
    conc1 = zeros(6)
    conc2 = zeros(6)
    conc1[1:5] = concs[:, d1ind]
    conc1[6] = conc2[6] = 10000.0
    conc2[1:5] = concs[:, d2ind]
    for i=1:6
        for j=1:6
            combined_effs[i,j] = low(conc1[i], conc2[j], pars1, pars2)
        end
    end
    return combined_effs[1:5, 1:5]
end

function output_loewe()
    concs = zeros(5, 5)
    concs[1, :] .= 0.0
    concs[2:5, [1, 5]] .= [25.0, 50.0, 100.0, 250]
    concs[2:5, 2] .= 20.0
    concs[2:5, 3] .= [5.0, 10.0, 17.0, 30.0]
    concs[2:5, 4] .= 2.0

    # data import
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
    lapatinib_palbo = loweCellNum(concs, 1, 5, Total)
    lapatinib_gemc = loweCellNum(concs, 1, 3, Total)
    gemc_palbo = loweCellNum(concs, 3, 5, Total)
    dox_gem = loweCellNum(concs, 2, 3, Total)

    df1 = DataFrames.DataFrame(plb50_lpt25=lapatinib_palbo[2, 2], plb50_lpt50=lapatinib_palbo[3, 2], plb50_lpt100=lapatinib_palbo[4, 2], plb50_lpt250=lapatinib_palbo[5, 2])
    df2 = DataFrames.DataFrame(gem10_lpt25=lapatinib_gemc[2, 3], gem10_lpt50=lapatinib_gemc[3, 3], gem10_lpt100=lapatinib_gemc[4, 3], gem10_lpt250=lapatinib_gemc[5, 3])
    df3 = DataFrames.DataFrame(dox20_gem5=dox_gem[2, 2], dox20_gem10=dox_gem[2, 3], dox20_gem17=dox_gem[2, 4], dox20_gem30=dox_gem[2, 5])
    df4 = DataFrames.DataFrame(plb50_gem5=gemc_palbo[2, 3], plb50_gem10=gemc_palbo[3, 3], plb50_gem17=gemc_palbo[4, 3], plb50_gem30=gemc_palbo[5, 3])
    df5 = DataFrames.DataFrame(lpt100_palbo25=lapatinib_palbo[4, 2], lpt100_palbo50=lapatinib_palbo[4, 3], lpt100_palbo100=lapatinib_palbo[4, 4], lpt100_palbo250=lapatinib_palbo[4, 5])
    df6 = DataFrames.DataFrame(lpt100_gem5=lapatinib_gemc[4, 2], lpt100_gem10=lapatinib_gemc[4, 3], lpt100_gem17=lapatinib_gemc[4, 4], lpt100_gem30=lapatinib_gemc[4, 5])
    df7 = DataFrames.DataFrame(gem10_palbo25=gemc_palbo[3, 2], gem10_palbo50=gemc_palbo[3, 3], gem10_palbo100=gemc_palbo[3, 4], gem10_palbo250=gemc_palbo[3, 5])
    XLSX.writetable("Loewe_palbo_lpt.xlsx", df1, overwrite=true, sheetname="cell number")
    XLSX.writetable("Loewe_gem_lpt.xlsx", df2, overwrite=true, sheetname="cell number")
    XLSX.writetable("Loewe_dox_gem.xlsx", df3, overwrite=true, sheetname="cell number")
    XLSX.writetable("Loewe_palbo_gem.xlsx", df4, overwrite=true, sheetname="cell number")
    XLSX.writetable("Loewe_lapt_palbo.xlsx", df5, overwrite=true, sheetname="cell number")
    XLSX.writetable("Loewe_lpt_gem.xlsx", df6, overwrite=true, sheetname="cell number")
    XLSX.writetable("Loewe_gem_palb.xlsx", df7, overwrite=true, sheetname="cell number")
end
