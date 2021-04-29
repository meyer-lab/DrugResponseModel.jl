""" Figure 2: drug combination. """

# data import
gt1, gt2 = DrugResponseModel.import_combination("AU01001");
gt1_2, gt2_2 = DrugResponseModel.import_combination("AU01101");
gt1_3, gt2_3 = DrugResponseModel.import_combination("AU00901");
concs, _, g1s1, g2s1 = load(189, 1);
_, _, g1s2, g2s2 = load(189, 2);
_, _, g1s3, g2s3 = load(189, 3);

g1m = (g1s1 .+ g1s2 .+ g1s3) ./ 3; # pure single drug data
g2m = (g2s1 .+ g2s2 .+ g2s3) ./ 3; # pure dingle drug data
totalm = g1m .+ g2m;

GS1 = cat(gt1, gt1_2, gt1_3, dims = 4);
GS2 = cat(gt2, gt2_2, gt2_3, dims = 4);

meanGS1 = mean(GS1, dims = 4);
meanGS2 = mean(GS2, dims = 4);
meanGS2[:, :, 19] .= mean(cat(gt2[:, :, 19], gt2_2[:, :, 19], dims = 3), dims = 3)[:, :, 1]

GC = zeros(2, length(meanGS1[1, :, 1, 1]), 6, 6)

GC[1:2, :, 1, :] .= meanGS1[1:2, :, 1] # control

# to find the Bliss_on_cellnumber for gem17 and existing combinations
function single_cellnum_combo(total1, total2, control1, control2)
    total1 = 1.0 .- (total1 ./ control1)
    total2 = 1.0 .- (total2 ./ control2)
    cmb = zeros(189)
    for i = 1:189
        cmb[i] = -(total1[i] + total2[i] .- (total1[i] * total2[i]) - 1.0) * (control1[i] + control2[i]) / 2
    end
    cmb
end

gem17 = meanGS1[3, 1:189, 21]
palbo50_gem17_cellnum = single_cellnum_combo(gem17, totalm[:, 5, 5], totalm[:, 1, 3], totalm[:, 1, 5])
lap100_gem17_cellnum = single_cellnum_combo(gem17, totalm[:, 6, 1], totalm[:, 1, 3], totalm[:, 1, 1])
dox20_gem17_cellnum = single_cellnum_combo(gem17, totalm[:, 4, 2], totalm[:, 1, 3], totalm[:, 1, 2])

t = LinRange(0.0, 95.0, 189)

# params from fitting all 5 drugs at once
p = [
    44.184,
    1.24076,
    0.0692788,
    0.0460918,
    0.3822,
    0.854034,
    0.605391,
    0.771326,
    0.0138293,
    0.00183699,
    0.000293753,
    0.0127534,
    0.00011816,
    0.0142541,
    60.6069,
    0.899573,
    1.99993,
    0.0748216,
    1.99326,
    0.468332,
    1.99864,
    1.22536,
    0.000141615,
    0.0318616,
    0.000216899,
    8.80158e-7,
    0.598489,
    0.00110572,
    6.68492,
    2.05974,
    1.99936,
    0.167588,
    0.507586,
    0.316074,
    0.248084,
    0.826596,
    1.6164e-5,
    3.10987e-6,
    3.55996e-5,
    7.73526e-6,
    0.0774056,
    8.26708e-5,
    3.34656,
    2.83739,
    0.0907361,
    0.108245,
    1.9758,
    1.96985,
    1.9993,
    0.210137,
    0.0690636,
    1.30442e-5,
    0.0767181,
    0.00991078,
    6.87891e-5,
    1.45086e-5,
    18.2253,
    1.1841,
    1.00505,
    0.0735852,
    1.97326,
    0.783828,
    0.45769,
    1.99355,
    0.0519941,
    0.000533671,
    0.00204743,
    9.52975e-5,
    5.23806e-5,
    0.0677505,
    0.339953,
    0.403341,
    0.802518,
    0.470576,
    1.298,
    0.423103,
];
efcs = getODEparams(p, concs);

# Bliss on cell numbers over time
Bliss_cellnum1 = zeros(189, 8, 8, 10)
Bliss_cellnum2 = zeros(189, 8, 8, 10)
for i = 1:189
    Bliss_cellnum1[i, :, :, :] .= blissCellNum(g1m, g2m, i)[1] # G1
    Bliss_cellnum2[i, :, :, :] .= blissCellNum(g1m, g2m, i)[2] # G2
end
Bliss_cellnum = Bliss_cellnum1 .+ Bliss_cellnum2 # total

# Bliss on Model
LPT_PLB = DrugResponseModel.AllBliss_params(efcs[:, :, 1], efcs[:, :, 5])
LPT_GEM = DrugResponseModel.AllBliss_params(efcs[:, :, 1], efcs[:, :, 3])
LPT_TAX = DrugResponseModel.AllBliss_params(efcs[:, :, 1], efcs[:, :, 4])
GEM_PLB = DrugResponseModel.AllBliss_params(efcs[:, :, 3], efcs[:, :, 5])
DOX_GEM = DrugResponseModel.AllBliss_params(efcs[:, :, 2], efcs[:, :, 3])

########### Palbociclib 50 nM + lapatinibs [25 nM, 50 nM, 100 nM, 250 nM]
palbo50Lap = zeros(3, 189, 4)
for i = 1:4
    palbo50Lap[1, :, i], palbo50Lap[2, :, i], _ = predict(LPT_PLB[:, i + 3, 5], LPT_PLB[:, 1, 1], t)
end
palbo50Lap[3, :, :] .= palbo50Lap[1, :, :] .+ palbo50Lap[2, :, :]
# well 2: 3,4,5,6

GC[1:2, :, 2, 1] .= meanGS1[1:2, :, 8]
GC[1:2, :, 3:6, 1] .= meanGS2[1:2, :, 3:6]

########### Palbociclib 50 nM + gemcitabines [5 nM, 10 nM, 17 nM, 30 nM]
palbo50Gem = zeros(3, 189, 4)
for i = 1:2
    palbo50Gem[1, :, i], palbo50Gem[2, :, i], _ = predict(GEM_PLB[:, i + 4, 5], GEM_PLB[:, 1, 1], t)
end
# Interpolation to find the parameters for 17 nM.
hill(p, c) = p[2] + (p[3] - p[2]) / (1 + ((p[1] / c)^p[4]))
gemc_hillParams = zeros(12, 4) # [a1, a2, b1, b2, b3, b4, d1, d2, d3, d4, d5, d6] x [EC50, min, max, k]
gemc_hillParams[:, 1] .= p[15] # ec50
gemc_hillParams[:, 4] .= p[16] # k
gemc_hillParams[1:6, 2] = p[71:76]
gemc_hillParams[7:12, 2] .= 0.0
gemc_hillParams[:, 3] .= p[17:28]
GEM17 = zeros(12)
for i = 1:length(GEM17)
    GEM17[i] = hill(gemc_hillParams[i, :], 17.0)
end
GEM17_PLB = DrugResponseModel.Bliss_params_unit(GEM17, efcs[:, 5, 5], hcat(efcs[:, 1, 3], efcs[:, 1, 5]))
palbo50Gem[1, :, 4], palbo50Gem[2, :, 4], _ = predict(GEM_PLB[:, 7, 5], GEM_PLB[:, 1, 1], t)
palbo50Gem[1, :, 3], palbo50Gem[2, :, 3], _ = predict(GEM17_PLB, GEM_PLB[:, 1, 1], t)
palbo50Gem[3, :, :] .= palbo50Gem[1, :, :] .+ palbo50Gem[2, :, :]
# well 2: 21, 22, 23, 24
GC[1:2, :, 2, 2] .= meanGS1[1:2, :, 8]
GC[1:2, :, 3:6, 2] .= meanGS2[1:2, :, 21:24]

########### Gemcitabine 10 nM + palbociclibs [25 nM, 50 nM, 100 nM, 250 nM]
Gem10palbo = zeros(3, 189, 3)
for i = 2:3
    Gem10palbo[1, :, i], Gem10palbo[2, :, i], _ = predict(GEM_PLB[:, 6, i + 4], GEM_PLB[:, 1, 1], t)
end
Gem10palbo[1, :, 1], Gem10palbo[2, :, 1], _ = predict(GEM_PLB[:, 6, 4], GEM_PLB[:, 1, 1], t)
Gem10palbo[3, :, :] .= Gem10palbo[1, :, :] .+ Gem10palbo[2, :, :]
# well 1: 18, 22, 23, 24
GC[1:2, :, 2, 3] .= meanGS1[1:2, :, 20]
GC[1:2, :, 3, 3] .= meanGS2[1:2, :, 18]
GC[1:2, :, 4:6, 3] .= meanGS2[1:2, :, 22:24]

########### Gemcitabine 10 nM + lapatinibs [25 nM, 50 nM, 100 nM, 250 nM]
Gem10Lap = zeros(3, 189, 4)
for i = 1:4
    Gem10Lap[1, :, i], Gem10Lap[2, :, i], _ = predict(LPT_GEM[:, i + 3, 6], LPT_GEM[:, 1, 1], t)
end
Gem10Lap[3, :, :] .= Gem10Lap[1, :, :] .+ Gem10Lap[2, :, :]
# well 2: 9, 10, 11, 12
GC[1:2, :, 2, 4] .= meanGS1[1:2, :, 20]
GC[1:2, :, 3:6, 4] .= meanGS2[1:2, :, 9:12]

########## Dox 20 nM + gemcitabines [5 nM, 10 nM, 17 nM, 30 nM]
dox20gem = zeros(3, 189, 4)
for i = 1:2
    dox20gem[1, :, i], dox20gem[2, :, i], _ = predict(DOX_GEM[:, 4, i + 4], DOX_GEM[:, 1, 1], t)
end
dox20gem[1, :, 4], dox20gem[2, :, 4], _ = predict(DOX_GEM[:, 4, 7], DOX_GEM[:, 1, 1], t)
DOX10_GEM17 = DrugResponseModel.Bliss_params_unit(efcs[:, 4, 2], GEM17, hcat(efcs[:, 1, 2], efcs[:, 1, 3]))
dox20gem[1, :, 3], dox20gem[2, :, 3], _ = predict(DOX10_GEM17, DOX_GEM[:, 1, 1], t)
dox20gem[3, :, :] .= dox20gem[1, :, :] .+ dox20gem[2, :, :]
# well 2: 15, 16, 17, 18
GC[1:2, :, 2, 5] .= meanGS1[1:2, :, 6]
GC[1:2, :, 3:6, 5] .= meanGS2[1:2, :, 15:18]

########### Lap 100 nM + gemcitabines [17 nM, 30 nM]
lap100gem = zeros(3, 189, 2)
LAP100_GEM17 = DrugResponseModel.Bliss_params_unit(efcs[:, 6, 1], GEM17, hcat(efcs[:, 1, 1], efcs[:, 1, 3]))
lap100gem[1, :, 1], lap100gem[2, :, 1], _ = predict(LAP100_GEM17, LPT_GEM[:, 1, 1], t)
lap100gem[1, :, 2], lap100gem[2, :, 2], _ = predict(LPT_GEM[:, 6, 7], LPT_GEM[:, 1, 1], t)
lap100gem[3, :, :] .= lap100gem[1, :, :] .+ lap100gem[2, :, :]
# well 2: 13, 19

########### Lap 100 nM + palbociclibs [25 nM, 50 nM, 100 nM, 250 nM]
Lap100palbo = zeros(3, 189, 3)
for i = 2:3
    Lap100palbo[1, :, i], Lap100palbo[2, :, i], _ = predict(LPT_PLB[:, 6, i + 4], LPT_PLB[:, 1, 1], t)
end
Lap100palbo[1, :, 1], Lap100palbo[2, :, 1], _ = predict(LPT_PLB[:, 6, 4], LPT_PLB[:, 1, 1], t)
Lap100palbo[3, :, :] .= Lap100palbo[1, :, :] .+ Lap100palbo[2, :, :]
# well 2: 8, 5, 14, 20
GC[1:2, :, 2, 6] .= meanGS1[1:2, :, 8]
GC[1:2, :, 3, 6] .= meanGS2[1:2, :, 4]
GC[1:2, :, 4, 6] .= meanGS2[1:2, :, 5]
GC[1:2, :, 5, 6] .= meanGS2[1:2, :, 14]
GC[1:2, :, 6, 6] .= meanGS2[1:2, :, 20]

########### Pax 2 nM + Lapatinib [50 nM, 100 nM]
Pax2_lap = zeros(3, 189, 2)
Pax2_lap[1, :, 1], Pax2_lap[2, :, 1], _ = predict(LPT_TAX[:, 5, 4], LPT_TAX[:, 1, 1], t)
Pax2_lap[1, :, 2], Pax2_lap[2, :, 2], _ = predict(LPT_TAX[:, 6, 4], LPT_TAX[:, 1, 1], t)
Pax2_lap[3, :, :] .= Pax2_lap[1, :, :] .+ Pax2_lap[2, :, :]
# well 2: 2 (50 nM) well 1: 16 (100 nM)

### The following 4 lines are for saving the Bliss on model and cell numbers in excel files.
### Remember to add XLSX and DataFrames packages before running this.
# df1 = DataFrames.DataFrame(pax2_lap50 = Pax2_lap[3, :, 1], pax2_lap100 = Pax2_lap[3, :, 2], lap100_palb25 = Lap100palbo[3, :, 1], lap100_palb100 = Lap100palbo[3, :, 2], lap100_palb250 = Lap100palbo[3, :, 3], lap100_gem17 = lap100gem[3, :, 1], lap100_gem30 = lap100gem[3, :, 2], dox20_gem5 = dox20gem[3, :, 1], dox20_gem10 = dox20gem[3, :, 2], dox20_gem17 = dox20gem[3, :, 3], dox20_gem30 = dox20gem[3, :, 4], gem10_lap25 = Gem10Lap[3, :, 1], gem10_lap50 = Gem10Lap[3, :, 2], gem10_lap100 = Gem10Lap[3, :, 3], gem10_lap250 = Gem10Lap[3, :, 4], gem10_palbo25 = Gem10palbo[3, :, 1], gem10_palbo100 = Gem10palbo[3, :, 2], gem10_palbo250 = Gem10palbo[3, :, 3], palbo50_gem5 = palbo50Gem[3, :, 1], palbo50_gem10 = palbo50Gem[3, :, 2], palbo50_gem17 = palbo50Gem[3, :, 3], palbo50_gem30 = palbo50Gem[3, :, 4], plb50_lpt25 = palbo50Lap[3, :, 1], plb50_lpt50 = palbo50Lap[3, :, 2], plb50_lpt100 = palbo50Lap[3, :, 3], plb50_lpt250 = palbo50Lap[3, :, 4], )
# df2 = DataFrames.DataFrame(pax2_lap50 = Bliss_cellnum[:, 5, 4, 3], pax2_lap100 = Bliss_cellnum[:, 6, 4, 3], lap100_palb25 = Bliss_cellnum[:, 6, 4, 4], lap100_palb100 = Bliss_cellnum[:, 6, 6, 4], lap100_palb250 = Bliss_cellnum[:, 6, 7, 4], lap100_gem17=lap100_gem17_cellnum, lap100_gem30 = Bliss_cellnum[:, 6, 7, 2], dox20_gem5 = Bliss_cellnum[:, 4, 5, 5], dox20_gem10 = Bliss_cellnum[:, 4, 6, 5], dox20_gem17 = dox20_gem17_cellnum, dox20_gem30 = Bliss_cellnum[:, 4, 7, 5], gem10_lap25 = Bliss_cellnum[:, 4, 6, 2], gem10_lap50 = Bliss_cellnum[:, 5, 6, 2], gem10_lap100 = Bliss_cellnum[:, 6, 6, 2], gem10_lap250 = Bliss_cellnum[:, 7, 6, 2], gem10_palbo25 = Bliss_cellnum[:, 6, 4, 7], gem10_palbo100 = Bliss_cellnum[:, 6, 6, 7], gem10_palbo250 = Bliss_cellnum[:, 6, 7, 7], palbo50_gem5 = Bliss_cellnum[:, 5, 5, 7], palbo50_gem10 = Bliss_cellnum[:, 6, 5, 7], palbo50_gem17 = palbo50_gem17_cellnum, palbo50_gem30 = Bliss_cellnum[:, 7, 5, 7], plb50_lpt25 = Bliss_cellnum[:, 4, 5, 4], plb50_lpt50 = Bliss_cellnum[:, 5, 5, 4], plb50_lpt100 = Bliss_cellnum[:, 6, 5, 4], plb50_lpt250 = Bliss_cellnum[:, 7, 5, 4])

# XLSX.writetable("Bliss_model.xlsx", df1, overwrite=true, sheetname="cell number")
# XLSX.writetable("Bliss_cellNumber.xlsx", df2, overwrite=true, sheetname="Bliss")

function SSEs_combination()
    SSEs = zeros(2, 7) # dim1: exp - Bliss on cell number, dim2: exp - Bliss on model

    SSEs[1, 1] = norm(Bliss_cellnum[:, 4:7, 5, 4] - meanGS2[3, 1:189, 3:6]) / 4
    SSEs[2, 1] = norm(palbo50Lap[3, :, :] - meanGS2[3, 1:189, 3:6]) / 4

    # find Bliss on cell number for palbo50 and gem 17 nM:
    palbo50gem17_cellnum = DrugResponseModel.pair_cellnum_Bliss(hcat(totalm[:, 1, 5], totalm[:, 5, 5]), hcat(totalm[:, 1, 3], meanGS1[3, 1:189, 21]))
    SSEs[1, 2] =
        (
            norm(Bliss_cellnum[:, 4:5, 5, 9] - meanGS2[3, 1:189, 21:22]) +
            norm(Bliss_cellnum[:, 7, 5, 9] - meanGS2[3, 1:189, 24]) +
            norm(palbo50gem17_cellnum - meanGS2[3, 1:189, 23])
        ) / 4
    SSEs[2, 2] = norm(palbo50Gem[3, :, :] - meanGS2[3, 1:189, 21:24]) / 4

    SSEs[1, 3] = (norm(Bliss_cellnum[:, 6, 6:7, 9] - meanGS1[3, 1:189, 23:24]) + norm(Bliss_cellnum[:, 6, 4, 9] - meanGS1[3, 1:189, 18])) / 3
    SSEs[2, 3] = (norm(Gem10palbo[3, :, 2:3] - meanGS1[3, 1:189, 23:24]) + norm(Gem10palbo[3, :, 1] - meanGS1[3, 1:189, 18])) / 3

    SSEs[1, 4] = norm(Bliss_cellnum[:, 4:7, 6, 2] - meanGS2[3, 1:189, 9:12]) / 4
    SSEs[2, 4] = norm(Gem10Lap[3, :, :] - meanGS2[3, 1:189, 9:12]) / 4

    dox20gem17_cellnum = DrugResponseModel.pair_cellnum_Bliss(hcat(totalm[:, 1, 2], totalm[:, 4, 2]), hcat(totalm[:, 1, 3], meanGS1[3, 1:189, 21]))
    SSEs[1, 5] =
        (
            norm(Bliss_cellnum[:, 4, 5:6, 5] - meanGS2[3, 1:189, 15:16]) +
            norm(dox20gem17_cellnum - meanGS2[3, 1:189, 17]) +
            norm(Bliss_cellnum[:, 4, 7, 5] - meanGS2[3, 1:189, 18])
        ) / 4
    SSEs[2, 5] = (norm(dox20gem[3, :, :] - meanGS2[3, 1:189, 15:18])) / 4

    lap100gem17_cellnum = DrugResponseModel.pair_cellnum_Bliss(hcat(totalm[:, 1, 1], totalm[:, 6, 1]), hcat(totalm[:, 1, 3], meanGS1[3, 1:189, 19]))
    SSEs[1, 6] = (norm(lap100gem17_cellnum - meanGS2[3, 1:189, 13]) + norm(Bliss_cellnum[:, 4, 7, 2] - meanGS2[3, 1:189, 19])) / 2
    SSEs[2, 6] = (norm(lap100gem[3, :, 2] - meanGS2[3, 1:189, 13]) + norm(lap100gem[3, :, 2] - meanGS2[3, 1:189, 19])) / 2

    SSEs[1, 7] =
        (
            norm(Bliss_cellnum[:, 6, 4, 4] - meanGS2[3, 1:189, 8]) +
            norm(Bliss_cellnum[:, 6, 6, 4] - meanGS2[3, 1:189, 14]) +
            norm(Bliss_cellnum[:, 6, 7, 4] - meanGS2[3, 1:189, 20])
        ) / 3
    SSEs[2, 7] =
        (
            norm(Lap100palbo[3, :, 1] - meanGS2[3, 1:189, 8]) +
            norm(Lap100palbo[3, :, 2] - meanGS2[3, 1:189, 14]) +
            norm(Lap100palbo[3, :, 3] - meanGS2[3, 1:189, 20])
        ) / 3

    return SSEs
end

function helper(Bliss_model1, Bliss_model2, i, title, subPlabel, label, GS2)
    t = LinRange(0.0, 95.0, 189)
    p = plot(
        t,
        Bliss_model1 .+ Bliss_model2,
        label = "total $label",
        xlabel = "time [hr]",
        lw = 3,
        ylabel = "cell number",
        color = "black",
        legend = :topright,
        title = title,
        titlefont = Plots.font("Helvetica", 12),
        guidefont = Plots.font("Helvetica", 12),
        xtickfont = Plots.font("Helvetica", 12),
        ytickfont = Plots.font("Helvetica", 12),
        bottom_margin = 1.5cm,
        fg_legend = :transparent,
        top_margin = 1.25cm,
        left_margin = 1.25cm,
        right_margin = 1.25cm,
    )
    plot!(
        t,
        mean(GS2, dims = 4)[3, 1:189, i],
        ribbon = std(GS2, dims = 4)[3, 1:189, i],
        alpha = 0.1,
        label = "",
        legend = :topright,
        titlefont = Plots.font("Helvetica", 12),
        lw = 3,
        color = "gray",
        guidefont = Plots.font("Helvetica", 12),
        xtickfont = Plots.font("Helvetica", 12),
        ytickfont = Plots.font("Helvetica", 12),
        bottom_margin = 1.5cm,
        fg_legend = :transparent,
        top_margin = 1.25cm,
        left_margin = 1.25cm,
        right_margin = 1.25cm,
    )
    plot!(
        t,
        Bliss_model1,
        label = "G1 $label",
        xlabel = "time [hr]",
        lw = 3,
        ylabel = "cell number",
        color = "mediumseagreen",
        titlefont = Plots.font("Helvetica", 12),
        guidefont = Plots.font("Helvetica", 12),
        xtickfont = Plots.font("Helvetica", 12),
        ytickfont = Plots.font("Helvetica", 12),
        bottom_margin = 1.5cm,
        fg_legend = :transparent,
        top_margin = 1.25cm,
        left_margin = 1.25cm,
        right_margin = 1.25cm,
    )
    plot!(
        t,
        Bliss_model2,
        label = "G2 $label",
        titlefont = Plots.font("Helvetica", 12),
        lw = 3,
        color = "mediumpurple3",
        guidefont = Plots.font("Helvetica", 12),
        xtickfont = Plots.font("Helvetica", 12),
        ytickfont = Plots.font("Helvetica", 12),
        bottom_margin = 1.5cm,
        fg_legend = :transparent,
        top_margin = 1.25cm,
        left_margin = 1.25cm,
        right_margin = 1.25cm,
    )
    plot!(
        t,
        mean(GS2, dims = 4)[1, 1:189, i],
        ribbon = std(GS2, dims = 4)[1, 1:189, i],
        alpha = 0.1,
        label = "",
        titlefont = Plots.font("Helvetica", 12),
        lw = 3,
        color = "darkseagreen3",
        guidefont = Plots.font("Helvetica", 12),
        xtickfont = Plots.font("Helvetica", 12),
        ytickfont = Plots.font("Helvetica", 12),
        bottom_margin = 1.5cm,
        fg_legend = :transparent,
        top_margin = 1.25cm,
        left_margin = 1.25cm,
        right_margin = 1.25cm,
    )
    plot!(
        t,
        mean(GS2, dims = 4)[2, 1:189, i],
        ribbon = std(GS2, dims = 4)[2, 1:189, i],
        alpha = 0.1,
        label = "",
        titlefont = Plots.font("Helvetica", 12),
        lw = 3,
        color = "plum",
        guidefont = Plots.font("Helvetica", 12),
        xtickfont = Plots.font("Helvetica", 12),
        ytickfont = Plots.font("Helvetica", 12),
        bottom_margin = 1.5cm,
        fg_legend = :transparent,
        top_margin = 1.25cm,
        left_margin = 1.25cm,
        right_margin = 1.25cm,
    )
    annotate!(-20.0, 47.0, text(subPlabel, :black, :left, Plots.font("Helvetica Bold", 15)))
    ylims!((0.0, 40.0))
    p
end

function figureS2()

    p1 = helper(Bliss_cellnum1[:, 4, 6, 2], Bliss_cellnum2[:, 5, 6, 2], 10, "gem 10 nM & lapt 50 nM", "a", "cell num", GS2)
    p2 = helper(Bliss_cellnum1[:, 5, 6, 2], Bliss_cellnum2[:, 7, 6, 2], 12, "gem 10 nM & lapt 250 nM", "b", "cell num", GS2)
    p3 = helper(Bliss_cellnum1[:, 6, 5, 9], Bliss_cellnum2[:, 6, 5, 9], 22, "palbo 50 nM & gem 10 nM", "c", "cell num", GS2)
    p4 = helper(Bliss_cellnum1[:, 7, 5, 9], Bliss_cellnum2[:, 7, 5, 9], 24, "palbo 50 nM & gem 30 nM", "d", "cell num", GS2)
    p5 = helper(Bliss_cellnum1[:, 5, 4, 3], Bliss_cellnum2[:, 5, 4, 3], 2, "pax 2 nM & lapt 50 nM", "e", "cell num", GS2)
    p6 = helper(Bliss_cellnum1[:, 6, 4, 3], Bliss_cellnum2[:, 6, 4, 3], 16, "pax 2 nM & lapt 100 nM", "f", "cell num", GS1)
    p7 = helper(Gem10Lap[1, :, 2], Gem10Lap[2, :, 2], 10, "gem 10 nM & lapt 50 nM", "g", "model", GS2)
    p8 = helper(Gem10Lap[1, :, 4], Gem10Lap[2, :, 4], 12, "gem 10 nM & lapt 250 nM", "h", "model", GS2)
    p9 = helper(palbo50Gem[1, :, 2], palbo50Gem[2, :, 2], 22, "palbo 50 nM & gem 10 nM", "i", "model", GS2)
    p10 = helper(palbo50Gem[1, :, 4], palbo50Gem[2, :, 4], 24, "palbo 50 nM & gem 30 nM", "j", "model", GS2)
    p11 = helper(Pax2_lap[1, :, 1], Pax2_lap[2, :, 1], 2, "pax 2 nM & lapt 50 nM", "k", "model", GS2)
    p12 = helper(Pax2_lap[1, :, 2], Pax2_lap[2, :, 2], 16, "pax 2 nM & lapt 100 nM", "l", "model", GS1)
    fig = plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, size = (2000, 700), layout = (2, 6))

    savefig(fig, "figureS2.svg")
    save("GC.jld", "GC", GC)
end
