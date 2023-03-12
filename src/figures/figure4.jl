"""Function to plot the combination figure4 for paper."""


function output_combination()
    concs, popul1, g1s1, g2s1 = load(189, 1)
    t = LinRange(0.0, 95.0, 189)
    p = parameters()
    efcs = getODEparams(p, concs)

    # Bliss on Model
    LPT_DOX = DrugResponseModel.AllBliss_params(efcs[:, :, 1], efcs[:, :, 2])
    LPT_GEM = DrugResponseModel.AllBliss_params(efcs[:, :, 1], efcs[:, :, 3])
    LPT_TAX = DrugResponseModel.AllBliss_params(efcs[:, :, 1], efcs[:, :, 4])
    LPT_PLB = DrugResponseModel.AllBliss_params(efcs[:, :, 1], efcs[:, :, 5])
    DOX_GEM = DrugResponseModel.AllBliss_params(efcs[:, :, 2], efcs[:, :, 3])
    DOX_TAX = DrugResponseModel.AllBliss_params(efcs[:, :, 2], efcs[:, :, 4])
    DOX_PLB = DrugResponseModel.AllBliss_params(efcs[:, :, 2], efcs[:, :, 5])
    GEM_TAX = DrugResponseModel.AllBliss_params(efcs[:, :, 3], efcs[:, :, 4])
    GEM_PLB = DrugResponseModel.AllBliss_params(efcs[:, :, 3], efcs[:, :, 5])
    TAX_PLB = DrugResponseModel.AllBliss_params(efcs[:, :, 4], efcs[:, :, 5])


    ########### lapatinibs + Doxorubicins
    LapDox = zeros(3, 189, 8, 8)
    for i = 1:8
        for j = 1:8
            LapDox[1, :, i, j], LapDox[2, :, i, j], _ = predict(LPT_DOX[:, i, j], LPT_DOX[:, 1, 1], t)
        end
    end
    LapDox[3, :, :, :] .= LapDox[1, :, :, :] .+ LapDox[2, :, :, :]

    ########### lapatinibs + Gemcitabines
    LapGem = zeros(3, 189, 8, 8)
    for i = 1:8
        for j = 1:8
            LapGem[1, :, i, j], LapGem[2, :, i, j], _ = predict(LPT_GEM[:, i, j], LPT_GEM[:, 1, 1], t)
        end
    end
    LapGem[3, :, :, :] .= LapGem[1, :, :, :] .+ LapGem[2, :, :, :]

    ########### lapatinibs + Taxol
    LapTax = zeros(3, 189, 8, 8)
    for i = 1:8
        for j = 1:8
            LapTax[1, :, i, j], LapTax[2, :, i, j], _ = predict(LPT_TAX[:, i, j], LPT_TAX[:, 1, 1], t)
        end
    end
    LapTax[3, :, :, :] .= LapTax[1, :, :, :] .+ LapTax[2, :, :, :]


    ########### lapatinibs + Palbociclibs
    LapPalbo = zeros(3, 189, 8, 8)
    for i = 1:8
        for j = 1:8
            LapPalbo[1, :, i, j], LapPalbo[2, :, i, j], _ = predict(LPT_PLB[:, i, j], LPT_PLB[:, 1, 1], t)
        end
    end
    LapPalbo[3, :, :, :] .= LapPalbo[1, :, :, :] .+ LapPalbo[2, :, :, :]

    ########### doxorubicins + gemcitabines 
    DoxGem = zeros(3, 189, 8, 8)
    for i = 1:8
        for j = 1:8
            DoxGem[1, :, i, j], DoxGem[2, :, i, j], _ = predict(DOX_GEM[:, i, j], DOX_GEM[:, 1, 1], t)
        end
    end
    DoxGem[3, :, :, :] .= DoxGem[1, :, :, :] .+ DoxGem[2, :, :, :]

    ########### doxorubicins + taxols 
    DoxTax = zeros(3, 189, 8, 8)
    for i = 1:8
        for j = 1:8
            DoxTax[1, :, i, j], DoxTax[2, :, i, j], _ = predict(DOX_TAX[:, i, j], DOX_TAX[:, 1, 1], t)
        end
    end
    DoxTax[3, :, :, :] .= DoxTax[1, :, :, :] .+ DoxTax[2, :, :, :]


    ########### doxorubicins + palbos
    DoxPalbo = zeros(3, 189, 8, 8)
    for i = 1:8
        for j = 1:8
            DoxPalbo[1, :, i, j], DoxPalbo[2, :, i, j], _ = predict(DOX_PLB[:, i, j], DOX_PLB[:, 1, 1], t)
        end
    end
    DoxPalbo[3, :, :, :] .= DoxPalbo[1, :, :, :] .+ DoxPalbo[2, :, :, :]

    ########### gemcitabines + taxols
    GemTax = zeros(3, 189, 8, 8)
    for i = 1:8
        for j = 1:8
            GemTax[1, :, i, j], GemTax[2, :, i, j], _ = predict(GEM_TAX[:, i, j], GEM_TAX[:, 1, 1], t)
        end
    end
    GemTax[3, :, :, :] .= GemTax[1, :, :, :] .+ GemTax[2, :, :, :]

    ########### gemcitabines + palbos
    GemPalb = zeros(3, 189, 8, 8)
    for i = 1:8
        for j = 1:8
            GemPalb[1, :, i, j], GemPalb[2, :, i, j], _ = predict(GEM_PLB[:, i, j], GEM_PLB[:, 1, 1], t)
        end
    end
    GemPalb[3, :, :, :] .= GemPalb[1, :, :, :] .+ GemPalb[2, :, :, :]

    ########### taxol + palbos
    TaxPlb = zeros(3, 189, 8, 8)
    for i = 1:8
        for j = 1:8
            TaxPlb[1, :, i, j], TaxPlb[2, :, i, j], _ = predict(TAX_PLB[:, i, j], TAX_PLB[:, 1, 1], t)
        end
    end
    TaxPlb[3, :, :, :] .= TaxPlb[1, :, :, :] .+ TaxPlb[2, :, :, :]

    # total
    df1 = DataFrames.DataFrame(
        control_control = LapDox[3, :, 1, 1],
        dox1 = LapDox[3, :, 1, 2],
        dox10 = LapDox[3, :, 1, 3],
        dox20 = LapDox[3, :, 1, 4],
        dox50 = LapDox[3, :, 1, 5],
        dox125 = LapDox[3, :, 1, 6],
        dox250 = LapDox[3, :, 1, 7],
        dox500 = LapDox[3, :, 1, 8],
        lpt5 = LapDox[3, :, 2, 1],
        lpt5dox1 = LapDox[3, :, 2, 2],
        lpt5dox10 = LapDox[3, :, 2, 3],
        lpt5dox20 = LapDox[3, :, 2, 4],
        lpt5dox50 = LapDox[3, :, 2, 5],
        lpt5dox100 = LapDox[3, :, 2, 6],
        lpt5dox250 = LapDox[3, :, 2, 7],
        lpt5dox500 = LapDox[3, :, 2, 8],
        lpt10 = LapDox[3, :, 3, 1],
        lpt10dox1 = LapDox[3, :, 3, 2],
        lpt10dox10 = LapDox[3, :, 3, 3],
        lpt10dox20 = LapDox[3, :, 3, 4],
        lpt10dox50 = LapDox[3, :, 3, 5],
        lpt10dox100 = LapDox[3, :, 3, 6],
        lpt10dox250 = LapDox[3, :, 3, 7],
        lpt10dox500 = LapDox[3, :, 3, 8],
        lpt25 = LapDox[3, :, 4, 1],
        lpt25dox1 = LapDox[3, :, 4, 2],
        lpt25dox10 = LapDox[3, :, 4, 3],
        lpt25dox20 = LapDox[3, :, 4, 4],
        lpt25dox50 = LapDox[3, :, 4, 5],
        lpt25dox100 = LapDox[3, :, 4, 6],
        lpt25dox250 = LapDox[3, :, 4, 7],
        lpt25dox500 = LapDox[3, :, 4, 8],
        lpt50 = LapDox[3, :, 5, 1],
        lpt50dox1 = LapDox[3, :, 5, 2],
        lpt50dox10 = LapDox[3, :, 5, 3],
        lpt50dox20 = LapDox[3, :, 5, 4],
        lpt50dox50 = LapDox[3, :, 5, 5],
        lpt50dox125 = LapDox[3, :, 5, 6],
        lpt50dox250 = LapDox[3, :, 5, 7],
        lpt50dox500 = LapDox[3, :, 5, 8],
        lpt100 = LapDox[3, :, 6, 1],
        lpt100dox1 = LapDox[3, :, 6, 2],
        lpt100dox10 = LapDox[3, :, 6, 3],
        lpt100dox20 = LapDox[3, :, 6, 4],
        lpt100dox50 = LapDox[3, :, 6, 5],
        lpt100dox125 = LapDox[3, :, 6, 6],
        lpt100dox250 = LapDox[3, :, 6, 7],
        lpt100dox500 = LapDox[3, :, 6, 8],
        lpt250 = LapDox[3, :, 7, 1],
        lpt250dox1 = LapDox[3, :, 7, 2],
        lpt250dox10 = LapDox[3, :, 7, 3],
        lpt250dox20 = LapDox[3, :, 7, 4],
        lpt250dox50 = LapDox[3, :, 7, 5],
        lpt250dox125 = LapDox[3, :, 7, 6],
        lpt250dox250 = LapDox[3, :, 7, 7],
        lpt250dox500 = LapDox[3, :, 7, 8],
        lpt500 = LapDox[3, :, 8, 1],
        lpt500dox1 = LapDox[3, :, 8, 2],
        lpt500dox10 = LapDox[3, :, 8, 3],
        lpt500dox20 = LapDox[3, :, 8, 4],
        lpt500dox50 = LapDox[3, :, 8, 5],
        lpt500dox125 = LapDox[3, :, 8, 6],
        lpt500dox250 = LapDox[3, :, 8, 7],
        lpt500dox500 = LapDox[3, :, 8, 8],
    )

    XLSX.writetable("LPT_DOX.xlsx", df1)

    df2 = DataFrames.DataFrame(
        control_control = LapGem[3, :, 1, 1],
        gem0_25 = LapGem[3, :, 1, 2],
        gem1 = LapGem[3, :, 1, 3],
        gem2_5 = LapGem[3, :, 1, 4],
        gem5 = LapGem[3, :, 1, 5],
        gem10 = LapGem[3, :, 1, 6],
        gem17 = LapGem[3, :, 1, 7],
        gem30 = LapGem[3, :, 1, 8],
        lpt5gem0_25 = LapGem[3, :, 2, 2],
        lpt5gem1 = LapGem[3, :, 2, 3],
        lpt5gem2_5 = LapGem[3, :, 2, 4],
        lpt5gem5 = LapGem[3, :, 2, 5],
        lpt5gem10 = LapGem[3, :, 2, 6],
        lpt5gem17 = LapGem[3, :, 2, 7],
        lpt5gem30 = LapGem[3, :, 2, 8],
        lpt10gem0_25 = LapGem[3, :, 3, 2],
        lpt10gem1 = LapGem[3, :, 3, 3],
        lpt10gem2_5 = LapGem[3, :, 3, 4],
        lpt10gem5 = LapGem[3, :, 3, 5],
        lpt10gem10 = LapGem[3, :, 3, 6],
        lpt10gem17 = LapGem[3, :, 3, 7],
        lpt10gem30 = LapGem[3, :, 3, 8],
        lpt25gem0_25 = LapGem[3, :, 4, 2],
        lpt25gem1 = LapGem[3, :, 4, 3],
        lpt25gem2_5 = LapGem[3, :, 4, 4],
        lpt25gem5 = LapGem[3, :, 4, 5],
        lpt25gem10 = LapGem[3, :, 4, 6],
        lpt25gem17 = LapGem[3, :, 4, 7],
        lpt25gem30 = LapGem[3, :, 4, 8],
        lpt50gem0_25 = LapGem[3, :, 5, 2],
        lpt50gem1 = LapGem[3, :, 5, 3],
        lpt50gem2_5 = LapGem[3, :, 5, 4],
        lpt50gem5 = LapGem[3, :, 5, 5],
        lpt50gem10 = LapGem[3, :, 5, 6],
        lpt50gem17 = LapGem[3, :, 5, 7],
        lpt50gem30 = LapGem[3, :, 5, 8],
        lpt100gem0_25 = LapGem[3, :, 6, 2],
        lpt100gem1 = LapGem[3, :, 6, 3],
        lpt100gem2_5 = LapGem[3, :, 6, 4],
        lpt100gem5 = LapGem[3, :, 6, 5],
        lpt100gem10 = LapGem[3, :, 6, 6],
        lpt100gem17 = LapGem[3, :, 6, 7],
        lpt100gem30 = LapGem[3, :, 6, 8],
        lpt250gem0_25 = LapGem[3, :, 7, 2],
        lpt250gem1 = LapGem[3, :, 7, 3],
        lpt250gem2_5 = LapGem[3, :, 7, 4],
        lpt250gem5 = LapGem[3, :, 7, 5],
        lpt250gem10 = LapGem[3, :, 7, 6],
        lpt250gem17 = LapGem[3, :, 7, 7],
        lpt250gem30 = LapGem[3, :, 7, 8],
        lpt500gem0_25 = LapGem[3, :, 8, 2],
        lpt500gem1 = LapGem[3, :, 8, 3],
        lpt500gem2_5 = LapGem[3, :, 8, 4],
        lpt500gem5 = LapGem[3, :, 8, 5],
        lpt500gem10 = LapGem[3, :, 8, 6],
        lpt500gem17 = LapGem[3, :, 8, 7],
        lpt500gem30 = LapGem[3, :, 8, 8],
    )

    XLSX.writetable("LPT_GEM.xlsx", df2)

    df3 = DataFrames.DataFrame(
        control_control = LapTax[3, :, 1, 1],
        taxol0_1 = LapTax[3, :, 1, 2],
        taxol1 = LapTax[3, :, 1, 3],
        taxol2 = LapTax[3, :, 1, 4],
        taxol3 = LapTax[3, :, 1, 5],
        taxol5 = LapTax[3, :, 1, 6],
        taxol7_5 = LapTax[3, :, 1, 7],
        taxol15 = LapTax[3, :, 1, 8],
        lpt5taxol0_1 = LapTax[3, :, 2, 2],
        lpt5taxol1 = LapTax[3, :, 2, 3],
        lpt5taxol2 = LapTax[3, :, 2, 4],
        lpt5taxol3 = LapTax[3, :, 2, 5],
        lpt5taxol5 = LapTax[3, :, 2, 6],
        lpt5taxol7_5 = LapTax[3, :, 2, 7],
        lpt5taxol15 = LapTax[3, :, 2, 8],
        lpt10taxol0_1 = LapTax[3, :, 3, 2],
        lpt10taxol1 = LapTax[3, :, 3, 3],
        lpt10taxol2 = LapTax[3, :, 3, 4],
        lpt10taxol3 = LapTax[3, :, 3, 5],
        lpt10taxol5 = LapTax[3, :, 3, 6],
        lpt10taxol7_5 = LapTax[3, :, 3, 7],
        lpt10taxol15 = LapTax[3, :, 3, 8],
        lpt25taxol0_1 = LapTax[3, :, 4, 2],
        lpt25taxol1 = LapTax[3, :, 4, 3],
        lpt25taxol2 = LapTax[3, :, 4, 4],
        lpt25taxol3 = LapTax[3, :, 4, 5],
        lpt25taxol5 = LapTax[3, :, 4, 6],
        lpt25taxol7_5 = LapTax[3, :, 4, 7],
        lpt25taxol15 = LapTax[3, :, 4, 8],
        lpt50taxol0_1 = LapTax[3, :, 5, 2],
        lpt50taxol1 = LapTax[3, :, 5, 3],
        lpt50taxol2 = LapTax[3, :, 5, 4],
        lpt50taxol3 = LapTax[3, :, 5, 5],
        lpt50taxol5 = LapTax[3, :, 5, 6],
        lpt50taxol7_5 = LapTax[3, :, 5, 7],
        lpt50taxol15 = LapTax[3, :, 5, 8],
        lpt100taxol0_1 = LapTax[3, :, 6, 2],
        lpt100taxol1 = LapTax[3, :, 6, 3],
        lpt100taxol2 = LapTax[3, :, 6, 4],
        lpt100taxol3 = LapTax[3, :, 6, 5],
        lpt100taxol5 = LapTax[3, :, 6, 6],
        lpt100taxol7_5 = LapTax[3, :, 6, 7],
        lpt100taxol15 = LapTax[3, :, 6, 8],
        lpt250taxol0_1 = LapTax[3, :, 7, 2],
        lpt250taxol1 = LapTax[3, :, 7, 3],
        lpt250taxol2 = LapTax[3, :, 7, 4],
        lpt250taxol3 = LapTax[3, :, 7, 5],
        lpt250taxol5 = LapTax[3, :, 7, 6],
        lpt250taxol7_5 = LapTax[3, :, 7, 7],
        lpt250taxol15 = LapTax[3, :, 7, 8],
        lpt500taxol0_1 = LapTax[3, :, 8, 2],
        lpt500taxol1 = LapTax[3, :, 8, 3],
        lpt500taxol2 = LapTax[3, :, 8, 4],
        lpt500taxol3 = LapTax[3, :, 8, 5],
        lpt500taxol5 = LapTax[3, :, 8, 6],
        lpt500taxol7_5 = LapTax[3, :, 8, 7],
        lpt500taxol15 = LapTax[3, :, 8, 8],
    )

    XLSX.writetable("LPT_TAX.xlsx", df3)

    df4 = DataFrames.DataFrame(
        control_control = LapPalbo[3, :, 1, 1],
        palbo5 = LapPalbo[3, :, 1, 2],
        palbo10 = LapPalbo[3, :, 1, 3],
        palbo25 = LapPalbo[3, :, 1, 4],
        palbo50 = LapPalbo[3, :, 1, 5],
        palbo100 = LapPalbo[3, :, 1, 6],
        palbo250 = LapPalbo[3, :, 1, 7],
        palbo500 = LapPalbo[3, :, 1, 8],
        lpt5palbo5 = LapPalbo[3, :, 2, 2],
        lpt5palbo10 = LapPalbo[3, :, 2, 3],
        lpt5palbo25 = LapPalbo[3, :, 2, 4],
        lpt5palbo50 = LapPalbo[3, :, 2, 5],
        lpt5palbo100 = LapPalbo[3, :, 2, 6],
        lptpalbo250 = LapPalbo[3, :, 2, 7],
        lpt5palbo500 = LapPalbo[3, :, 2, 8],
        lpt10palbo5 = LapPalbo[3, :, 3, 2],
        lpt10palbo10 = LapPalbo[3, :, 3, 3],
        lpt10palbo25 = LapPalbo[3, :, 3, 4],
        lpt10palbo50 = LapPalbo[3, :, 3, 5],
        lpt10palbo100 = LapPalbo[3, :, 3, 6],
        lpt10palbo250 = LapPalbo[3, :, 3, 7],
        lpt10palbo500 = LapPalbo[3, :, 3, 8],
        lpt25palbo5 = LapPalbo[3, :, 4, 2],
        lpt25palbo10 = LapPalbo[3, :, 4, 3],
        lpt25palbo25 = LapPalbo[3, :, 4, 4],
        lpt25palbo50 = LapPalbo[3, :, 4, 5],
        lpt25palbo100 = LapPalbo[3, :, 4, 6],
        lpt25palbo250 = LapPalbo[3, :, 4, 7],
        lpt25palbo500 = LapPalbo[3, :, 4, 8],
        lpt50palbo5 = LapPalbo[3, :, 5, 2],
        lpt50palbo10 = LapPalbo[3, :, 5, 3],
        lpt50palbo25 = LapPalbo[3, :, 5, 4],
        lpt50palbo50 = LapPalbo[3, :, 5, 5],
        lpt50palbo100 = LapPalbo[3, :, 5, 6],
        lpt50palbo250 = LapPalbo[3, :, 5, 7],
        lpt50palbo500 = LapPalbo[3, :, 5, 8],
        lpt100palbo5 = LapPalbo[3, :, 6, 2],
        lpt100palbo10 = LapPalbo[3, :, 6, 3],
        lpt100palbo25 = LapPalbo[3, :, 6, 4],
        lpt100palbo50 = LapPalbo[3, :, 6, 5],
        lpt100palbo100 = LapPalbo[3, :, 6, 6],
        lpt100palbo250 = LapPalbo[3, :, 6, 7],
        lpt100palbo500 = LapPalbo[3, :, 6, 8],
        lpt250palbo5 = LapPalbo[3, :, 7, 2],
        lpt250palbo10 = LapPalbo[3, :, 7, 3],
        lpt250palbo25 = LapPalbo[3, :, 7, 4],
        lpt250palbo50 = LapPalbo[3, :, 7, 5],
        lpt250palbo100 = LapPalbo[3, :, 7, 6],
        lpt250palbo250 = LapPalbo[3, :, 7, 7],
        lpt250palbo500 = LapPalbo[3, :, 7, 8],
        lpt500palbo5 = LapPalbo[3, :, 8, 2],
        lpt500palbo10 = LapPalbo[3, :, 8, 3],
        lpt500palbo25 = LapPalbo[3, :, 8, 4],
        lpt500palbo50 = LapPalbo[3, :, 8, 5],
        lpt500palbo100 = LapPalbo[3, :, 8, 6],
        lpt500palbo250 = LapPalbo[3, :, 8, 7],
        lpt500palbo500 = LapPalbo[3, :, 8, 8],
    )

    XLSX.writetable("LPT_PLB.xlsx", df4)

    df5 = DataFrames.DataFrame(
        control_control = DoxGem[3, :, 1, 1],
        dox1gem0_25 = DoxGem[3, :, 2, 2],
        dox1gem1 = DoxGem[3, :, 2, 3],
        dox1gem2_5 = DoxGem[3, :, 2, 4],
        dox1gem5 = DoxGem[3, :, 2, 5],
        dox1gem10 = DoxGem[3, :, 2, 6],
        dox1gem17 = DoxGem[3, :, 2, 7],
        dox1gem30 = DoxGem[3, :, 2, 8],
        dox10gem0_25 = DoxGem[3, :, 3, 2],
        dox10gem1 = DoxGem[3, :, 3, 3],
        dox10gem2_5 = DoxGem[3, :, 3, 4],
        dox10gem5 = DoxGem[3, :, 3, 5],
        dox10gem10 = DoxGem[3, :, 3, 6],
        dox10gem17 = DoxGem[3, :, 3, 7],
        dox10gem30 = DoxGem[3, :, 3, 8],
        dox20gem0_25 = DoxGem[3, :, 4, 2],
        dox20gem1 = DoxGem[3, :, 4, 3],
        dox20gem2_5 = DoxGem[3, :, 4, 4],
        dox20gem5 = DoxGem[3, :, 4, 5],
        dox20gem10 = DoxGem[3, :, 4, 6],
        dox20gem17 = DoxGem[3, :, 4, 7],
        dox20gem30 = DoxGem[3, :, 4, 8],
        dox50gem0_25 = DoxGem[3, :, 5, 2],
        dox50gem1 = DoxGem[3, :, 5, 3],
        dox50gem2_5 = DoxGem[3, :, 5, 4],
        dox50gem5 = DoxGem[3, :, 5, 5],
        dox50gem10 = DoxGem[3, :, 5, 6],
        dox50gem17 = DoxGem[3, :, 5, 7],
        dox50gem30 = DoxGem[3, :, 5, 8],
        dox125em0_25 = DoxGem[3, :, 6, 2],
        dox125gem1 = DoxGem[3, :, 6, 3],
        dox125gem2_5 = DoxGem[3, :, 6, 4],
        dox125gem5 = DoxGem[3, :, 6, 5],
        dox125gem10 = DoxGem[3, :, 6, 6],
        dox125gem17 = DoxGem[3, :, 6, 7],
        dox125gem30 = DoxGem[3, :, 6, 8],
        dox250gem0_25 = DoxGem[3, :, 7, 2],
        dox250gem1 = DoxGem[3, :, 7, 3],
        dox250gem2_5 = DoxGem[3, :, 7, 4],
        dox250gem5 = DoxGem[3, :, 7, 5],
        dox250gem10 = DoxGem[3, :, 7, 6],
        dox250gem17 = DoxGem[3, :, 7, 7],
        dox250gem30 = DoxGem[3, :, 7, 8],
        dox500gem0_25 = DoxGem[3, :, 8, 2],
        dox500gem1 = DoxGem[3, :, 8, 3],
        dox500gem2_5 = DoxGem[3, :, 8, 4],
        dox500gem5 = DoxGem[3, :, 8, 5],
        dox500gem10 = DoxGem[3, :, 8, 6],
        dox500gem17 = DoxGem[3, :, 8, 7],
        dox500gem30 = DoxGem[3, :, 8, 8],
    )

    XLSX.writetable("DOX_GEM.xlsx", df5)

    df6 = DataFrames.DataFrame(
        control_control = DoxTax[3, :, 1, 1],
        dox1taxol0_1 = DoxTax[3, :, 2, 2],
        dox1taxol1 = DoxTax[3, :, 2, 3],
        dox1taxol2 = DoxTax[3, :, 2, 4],
        dox1taxol3 = DoxTax[3, :, 2, 5],
        dox1taxol5 = DoxTax[3, :, 2, 6],
        dox1taxol7_5 = DoxTax[3, :, 2, 7],
        dox1taxol15 = DoxTax[3, :, 2, 8],
        dox10taxol0_1 = DoxTax[3, :, 3, 2],
        dox10taxol1 = DoxTax[3, :, 3, 3],
        dox10taxol2 = DoxTax[3, :, 3, 4],
        dox10taxol3 = DoxTax[3, :, 3, 5],
        dox10taxol5 = DoxTax[3, :, 3, 6],
        dox10taxol7_5 = DoxTax[3, :, 3, 7],
        dox10taxol15 = DoxTax[3, :, 3, 8],
        dox20taxol0_1 = DoxTax[3, :, 4, 2],
        dox20taxol1 = DoxTax[3, :, 4, 3],
        dox20taxol2 = DoxTax[3, :, 4, 4],
        dox20taxol3 = DoxTax[3, :, 4, 5],
        dox20taxol5 = DoxTax[3, :, 4, 6],
        dox20taxol7_5 = DoxTax[3, :, 4, 7],
        dox20taxol15 = DoxTax[3, :, 4, 8],
        dox50taxol0_1 = DoxTax[3, :, 5, 2],
        dox50taxol1 = DoxTax[3, :, 5, 3],
        dox50taxol2 = DoxTax[3, :, 5, 4],
        dox50taxol3 = DoxTax[3, :, 5, 5],
        dox50taxol5 = DoxTax[3, :, 5, 6],
        dox50taxol7_5 = DoxTax[3, :, 5, 7],
        dox50taxol15 = DoxTax[3, :, 5, 8],
        dox125taxol0_1 = DoxTax[3, :, 6, 2],
        dox125taxol1 = DoxTax[3, :, 6, 3],
        dox125taxol2 = DoxTax[3, :, 6, 4],
        dox125taxol3 = DoxTax[3, :, 6, 5],
        dox125taxol5 = DoxTax[3, :, 6, 6],
        dox125taxol7_5 = DoxTax[3, :, 6, 7],
        dox125taxol15 = DoxTax[3, :, 6, 8],
        dox250taxol0_1 = DoxTax[3, :, 7, 2],
        dox250taxol1 = DoxTax[3, :, 7, 3],
        dox250taxol2 = DoxTax[3, :, 7, 4],
        dox250taxol3 = DoxTax[3, :, 7, 5],
        dox250taxol5 = DoxTax[3, :, 7, 6],
        dox250taxol7_5 = DoxTax[3, :, 7, 7],
        dox250taxol15 = DoxTax[3, :, 7, 8],
        dox500taxol0_1 = DoxTax[3, :, 8, 2],
        dox500taxol1 = DoxTax[3, :, 8, 3],
        dox500taxol2 = DoxTax[3, :, 8, 4],
        dox500taxol3 = DoxTax[3, :, 8, 5],
        dox500taxol5 = DoxTax[3, :, 8, 6],
        dox500taxol7_5 = DoxTax[3, :, 8, 7],
        dox500taxol15 = DoxTax[3, :, 8, 8],
    )
    XLSX.writetable("DOX_TAX.xlsx", df6)

    df7 = DataFrames.DataFrame(
        control_control = DoxPalbo[3, :, 1, 1],
        dox1palbo5 = DoxPalbo[3, :, 2, 2],
        dox1palbo10 = DoxPalbo[3, :, 2, 3],
        dox1palbo25 = DoxPalbo[3, :, 2, 4],
        dox1palbo50 = DoxPalbo[3, :, 2, 5],
        dox1palbo100 = DoxPalbo[3, :, 2, 6],
        dox1palbo250 = DoxPalbo[3, :, 2, 7],
        dox1palbo500 = DoxPalbo[3, :, 2, 8],
        dox10palbo5 = DoxPalbo[3, :, 3, 2],
        dox10palbo10 = DoxPalbo[3, :, 3, 3],
        dox10palbo25 = DoxPalbo[3, :, 3, 4],
        dox10palbo50 = DoxPalbo[3, :, 3, 5],
        dox10palbo100 = DoxPalbo[3, :, 3, 6],
        dox10palbo250 = DoxPalbo[3, :, 3, 7],
        dox10palbo500 = DoxPalbo[3, :, 3, 8],
        dox20palbo5 = DoxPalbo[3, :, 4, 2],
        dox20palbo10 = DoxPalbo[3, :, 4, 3],
        dox20palbo25 = DoxPalbo[3, :, 4, 4],
        dox20palbo50 = DoxPalbo[3, :, 4, 5],
        dox20palbo100 = DoxPalbo[3, :, 4, 6],
        dox20palbo250 = DoxPalbo[3, :, 4, 7],
        dox20palbo500 = DoxPalbo[3, :, 4, 8],
        dox50palbo5 = DoxPalbo[3, :, 5, 2],
        dox50palbo10 = DoxPalbo[3, :, 5, 3],
        dox50palbo25 = DoxPalbo[3, :, 5, 4],
        dox50palbo50 = DoxPalbo[3, :, 5, 5],
        dox50palbo100 = DoxPalbo[3, :, 5, 6],
        dox50palbo250 = DoxPalbo[3, :, 5, 7],
        dox50palbo500 = DoxPalbo[3, :, 5, 8],
        dox125palbo5 = DoxPalbo[3, :, 6, 2],
        dox125palbo10 = DoxPalbo[3, :, 6, 3],
        dox125palbo25 = DoxPalbo[3, :, 6, 4],
        dox125palbo50 = DoxPalbo[3, :, 6, 5],
        dox125palbo100 = DoxPalbo[3, :, 6, 6],
        dox125palbo250 = DoxPalbo[3, :, 6, 7],
        dox125palbo500 = DoxPalbo[3, :, 6, 8],
        dox250palbo5 = DoxPalbo[3, :, 7, 2],
        dox250palbo10 = DoxPalbo[3, :, 7, 3],
        dox250palbo25 = DoxPalbo[3, :, 7, 4],
        dox250palbo50 = DoxPalbo[3, :, 7, 5],
        dox250palbo100 = DoxPalbo[3, :, 7, 6],
        dox250palbo250 = DoxPalbo[3, :, 7, 7],
        dox250palbo500 = DoxPalbo[3, :, 7, 8],
        dox500palbo5 = DoxPalbo[3, :, 8, 2],
        dox500palbo10 = DoxPalbo[3, :, 8, 3],
        dox500palbo25 = DoxPalbo[3, :, 8, 4],
        dox500palbo50 = DoxPalbo[3, :, 8, 5],
        dox500palbo100 = DoxPalbo[3, :, 8, 6],
        dox500palbo250 = DoxPalbo[3, :, 8, 7],
        dox500palbo500 = DoxPalbo[3, :, 8, 8],
    )

    XLSX.writetable("DOX_PLB.xlsx", df7)

    df8 = DataFrames.DataFrame(
        control_control = GemTax[3, :, 1, 1],
        gem0_25taxol0_1 = GemTax[3, :, 2, 2],
        gem0_25taxol1 = GemTax[3, :, 2, 3],
        gem0_25taxol2 = GemTax[3, :, 2, 4],
        gem0_25taxol3 = GemTax[3, :, 2, 5],
        gem0_25taxol5 = GemTax[3, :, 2, 6],
        gem0_25taxol7_5 = GemTax[3, :, 2, 7],
        gem0_25taxol15 = GemTax[3, :, 2, 8],
        gem1taxol0_1 = GemTax[3, :, 3, 2],
        gem1taxol1 = GemTax[3, :, 3, 3],
        gem1taxol2 = GemTax[3, :, 3, 4],
        gem1taxol3 = GemTax[3, :, 3, 5],
        gem1taxol5 = GemTax[3, :, 3, 6],
        gem1taxol7_5 = GemTax[3, :, 3, 7],
        gem1taxol15 = GemTax[3, :, 3, 8],
        gem2_5taxol0_1 = GemTax[3, :, 4, 2],
        gem2_5taxol1 = GemTax[3, :, 4, 3],
        gem2_5taxol2 = GemTax[3, :, 4, 4],
        gem2_5taxol3 = GemTax[3, :, 4, 5],
        gem2_5taxol5 = GemTax[3, :, 4, 6],
        gem2_5taxol7_5 = GemTax[3, :, 4, 7],
        gem2_5taxol15 = GemTax[3, :, 4, 8],
        gem5taxol0_1 = GemTax[3, :, 5, 2],
        gem5taxol1 = GemTax[3, :, 5, 3],
        gem5taxol2 = GemTax[3, :, 5, 4],
        gem5taxol3 = GemTax[3, :, 5, 5],
        gem5taxol5 = GemTax[3, :, 5, 6],
        gem5taxol7_5 = GemTax[3, :, 5, 7],
        gem5taxol15 = GemTax[3, :, 5, 8],
        gem10taxol0_1 = GemTax[3, :, 6, 2],
        gem10taxol1 = GemTax[3, :, 6, 3],
        gem10taxol2 = GemTax[3, :, 6, 4],
        gem10taxol3 = GemTax[3, :, 6, 5],
        gem10taxol5 = GemTax[3, :, 6, 6],
        gem10taxol7_5 = GemTax[3, :, 6, 7],
        gem10taxol15 = GemTax[3, :, 6, 8],
        gem17taxol0_1 = GemTax[3, :, 7, 2],
        gem17taxol1 = GemTax[3, :, 7, 3],
        gem17taxol2 = GemTax[3, :, 7, 4],
        gem17taxol3 = GemTax[3, :, 7, 5],
        gem17taxol5 = GemTax[3, :, 7, 6],
        gem17taxol7_5 = GemTax[3, :, 7, 7],
        gem17taxol15 = GemTax[3, :, 7, 8],
        gem30taxol0_1 = GemTax[3, :, 8, 2],
        gem30taxol1 = GemTax[3, :, 8, 3],
        gem30taxol2 = GemTax[3, :, 8, 4],
        gem30taxol3 = GemTax[3, :, 8, 5],
        gem30taxol5 = GemTax[3, :, 8, 6],
        gem30taxol7_5 = GemTax[3, :, 8, 7],
        gem30taxol15 = GemTax[3, :, 8, 8],
    )
    XLSX.writetable("GEM_TAX.xlsx", df8)

    df9 = DataFrames.DataFrame(
        control_control = GemPalb[3, :, 1, 1],
        gem0_25palbo5 = GemPalb[3, :, 2, 2],
        gem0_25palbo10 = GemPalb[3, :, 2, 3],
        gem0_25palbo25 = GemPalb[3, :, 2, 4],
        gem0_25palbo50 = GemPalb[3, :, 2, 5],
        gem0_25palbo100 = GemPalb[3, :, 2, 6],
        gem0_25palbo250 = GemPalb[3, :, 2, 7],
        gem0_25palbo500 = GemPalb[3, :, 2, 8],
        gem1palbo5 = GemPalb[3, :, 3, 2],
        gem1palbo10 = GemPalb[3, :, 3, 3],
        gem1palbo25 = GemPalb[3, :, 3, 4],
        gem1palbo50 = GemPalb[3, :, 3, 5],
        gem1palbo100 = GemPalb[3, :, 3, 6],
        gem1palbo250 = GemPalb[3, :, 3, 7],
        gem1palbo500 = GemPalb[3, :, 3, 8],
        gem2_5palbo5 = GemPalb[3, :, 4, 2],
        gem2_5palbo10 = GemPalb[3, :, 4, 3],
        gem2_5palbo25 = GemPalb[3, :, 4, 4],
        gem2_5palbo50 = GemPalb[3, :, 4, 5],
        gem2_5palbo100 = GemPalb[3, :, 4, 6],
        gem2_5palbo250 = GemPalb[3, :, 4, 7],
        gem2_5palbo500 = GemPalb[3, :, 4, 8],
        gem5palbo5 = GemPalb[3, :, 5, 2],
        gem5palbo10 = GemPalb[3, :, 5, 3],
        gem5palbo25 = GemPalb[3, :, 5, 4],
        gem5palbo50 = GemPalb[3, :, 5, 5],
        gem5palbo100 = GemPalb[3, :, 5, 6],
        gem5palbo250 = GemPalb[3, :, 5, 7],
        gem5palbo500 = GemPalb[3, :, 5, 8],
        gem10palbo5 = GemPalb[3, :, 6, 2],
        gem10palbo10 = GemPalb[3, :, 6, 3],
        gem10palbo25 = GemPalb[3, :, 6, 4],
        gem10palbo50 = GemPalb[3, :, 6, 5],
        gem10palbo100 = GemPalb[3, :, 6, 6],
        gem10palbo250 = GemPalb[3, :, 6, 7],
        gem10palbo500 = GemPalb[3, :, 6, 8],
        gem17palbo5 = GemPalb[3, :, 7, 2],
        gem17palbo10 = GemPalb[3, :, 7, 3],
        gem17palbo25 = GemPalb[3, :, 7, 4],
        gem17palbo50 = GemPalb[3, :, 7, 5],
        gem17palbo100 = GemPalb[3, :, 7, 6],
        gem17palbo250 = GemPalb[3, :, 7, 7],
        gem17palbo500 = GemPalb[3, :, 7, 8],
        gem30palbo5 = GemPalb[3, :, 8, 2],
        gem30palbo10 = GemPalb[3, :, 8, 3],
        gem30palbo25 = GemPalb[3, :, 8, 4],
        gem30palbo50 = GemPalb[3, :, 8, 5],
        gem30palbo100 = GemPalb[3, :, 8, 6],
        gem30palbo250 = GemPalb[3, :, 8, 7],
        gem30palbo500 = GemPalb[3, :, 8, 8],
    )

    XLSX.writetable("GEM_PLB.xlsx", df9)

    df10 = DataFrames.DataFrame(
        control_control = TaxPlb[3, :, 1, 1],
        tax0_1palbo5 = TaxPlb[3, :, 2, 2],
        tax0_1palbo10 = TaxPlb[3, :, 2, 3],
        tax0_1palbo25 = TaxPlb[3, :, 2, 4],
        tax0_1palbo50 = TaxPlb[3, :, 2, 5],
        tax0_1palbo100 = TaxPlb[3, :, 2, 6],
        tax0_1palbo250 = TaxPlb[3, :, 2, 7],
        tax0_1palbo500 = TaxPlb[3, :, 2, 8],
        tax1palbo5 = TaxPlb[3, :, 3, 2],
        tax1palbo10 = TaxPlb[3, :, 3, 3],
        tax1palbo25 = TaxPlb[3, :, 3, 4],
        tax1palbo50 = TaxPlb[3, :, 3, 5],
        tax1palbo100 = TaxPlb[3, :, 3, 6],
        tax1palbo250 = TaxPlb[3, :, 3, 7],
        tax1palbo500 = TaxPlb[3, :, 3, 8],
        tax2palbo5 = TaxPlb[3, :, 4, 2],
        tax2palbo10 = TaxPlb[3, :, 4, 3],
        tax2palbo25 = TaxPlb[3, :, 4, 4],
        tax2palbo50 = TaxPlb[3, :, 4, 5],
        tax2palbo100 = TaxPlb[3, :, 4, 6],
        tax2palbo250 = TaxPlb[3, :, 4, 7],
        tax2palbo500 = TaxPlb[3, :, 4, 8],
        tax3palbo5 = TaxPlb[3, :, 5, 2],
        tax3palbo10 = TaxPlb[3, :, 5, 3],
        tax3palbo25 = TaxPlb[3, :, 5, 4],
        tax3palbo50 = TaxPlb[3, :, 5, 5],
        tax3palbo100 = TaxPlb[3, :, 5, 6],
        tax3palbo250 = TaxPlb[3, :, 5, 7],
        tax3palbo500 = TaxPlb[3, :, 5, 8],
        tax5palbo5 = TaxPlb[3, :, 6, 2],
        tax5palbo10 = TaxPlb[3, :, 6, 3],
        tax5palbo25 = TaxPlb[3, :, 6, 4],
        tax5palbo50 = TaxPlb[3, :, 6, 5],
        tax5palbo100 = TaxPlb[3, :, 6, 6],
        tax5palbo250 = TaxPlb[3, :, 6, 7],
        tax5palbo500 = TaxPlb[3, :, 6, 8],
        taxol7_5palbo5 = TaxPlb[3, :, 7, 2],
        taxol7_5palbo10 = TaxPlb[3, :, 7, 3],
        taxol7_5palbo25 = TaxPlb[3, :, 7, 4],
        taxol7_5palbo50 = TaxPlb[3, :, 7, 5],
        taxol7_5palbo100 = TaxPlb[3, :, 7, 6],
        taxol7_5palbo250 = TaxPlb[3, :, 7, 7],
        taxol7_5palbo500 = TaxPlb[3, :, 7, 8],
        taxol15palbo5 = TaxPlb[3, :, 8, 2],
        taxol15palbo10 = TaxPlb[3, :, 8, 3],
        taxol15palbo25 = TaxPlb[3, :, 8, 4],
        taxol15palbo50 = TaxPlb[3, :, 8, 5],
        taxol15palbo100 = TaxPlb[3, :, 8, 6],
        taxol15palbo250 = TaxPlb[3, :, 8, 7],
        taxol15palbo500 = TaxPlb[3, :, 8, 8],
    )

    XLSX.writetable("TAX_PLB.xlsx", df10)
end


"""This function outputs the average duration of time that cells spend in G1 or S/G2 under combination treatments to excel files."""
function combined_phaseDurations()
    concs, popul1, g1s1, g2s1 = load(189, 1)
    t = LinRange(0.0, 95.0, 189)
    p = parameters()
    efcs = getODEparams(p, concs)
    # gem17 = DrugResponseModel.find_gem17(p)
    # efcs[:, 8, 3] = efcs[:, 7, 3]
    # efcs[:, 7, 3] = gem17

    # Bliss on Model
    LPT_DOX = DrugResponseModel.AllBliss_params(efcs[:, :, 1], efcs[:, :, 2])
    LPT_GEM = DrugResponseModel.AllBliss_params(efcs[:, :, 1], efcs[:, :, 3])
    LPT_TAX = DrugResponseModel.AllBliss_params(efcs[:, :, 1], efcs[:, :, 4])
    LPT_PLB = DrugResponseModel.AllBliss_params(efcs[:, :, 1], efcs[:, :, 5])
    DOX_GEM = DrugResponseModel.AllBliss_params(efcs[:, :, 2], efcs[:, :, 3])
    DOX_TAX = DrugResponseModel.AllBliss_params(efcs[:, :, 2], efcs[:, :, 4])
    DOX_PLB = DrugResponseModel.AllBliss_params(efcs[:, :, 2], efcs[:, :, 5])
    GEM_TAX = DrugResponseModel.AllBliss_params(efcs[:, :, 3], efcs[:, :, 4])
    GEM_PLB = DrugResponseModel.AllBliss_params(efcs[:, :, 3], efcs[:, :, 5])
    TAX_PLB = DrugResponseModel.AllBliss_params(efcs[:, :, 4], efcs[:, :, 5])

    function easier_(efcs)
        gi = zeros(2, 8, 8)
        gi[1, :, :] .= (2 ./ efcs[1, :, :] .+ 2 ./ efcs[2, :, :] .+ 2 ./ efcs[3, :, :] .+ 2 ./ efcs[4, :, :])
        gi[2, :, :] .= (5 ./ efcs[5, :, :] .+ 5 ./ efcs[6, :, :] .+ 5 ./ efcs[7, :, :] .+ 5 ./ efcs[8, :, :])
        gi
    end
    gim = zeros(2, 8, 8, 10)
    gim[:, :, :, 1] = easier_(LPT_DOX)
    gim[:, :, :, 2] = easier_(LPT_GEM)
    gim[:, :, :, 3] = easier_(LPT_TAX)
    gim[:, :, :, 4] = easier_(LPT_PLB)
    gim[:, :, :, 5] = easier_(DOX_GEM)
    gim[:, :, :, 6] = easier_(DOX_TAX)
    gim[:, :, :, 7] = easier_(DOX_PLB)
    gim[:, :, :, 8] = easier_(GEM_TAX)
    gim[:, :, :, 9] = easier_(GEM_PLB)
    gim[:, :, :, 10] = easier_(TAX_PLB)
    dfG1 = DataFrames.DataFrame(
        palbo50_lpt25 = gim[1, 4, 5, 4],
        palbo50_lpt50 = gim[1, 5, 5, 4],
        palbo50_lpt100 = gim[1, 6, 5, 4],
        palbo50_lpt250 = gim[1, 7, 5, 4],
        palbo50_gem5 = gim[1, 5, 5, 9],
        palbo50_gem10 = gim[1, 6, 5, 9],
        palbo50_gem30 = gim[1, 7, 5, 9],
        palbo50_gem100 = gim[1, 8, 5, 9],
        gem10_lpt25 = gim[1, 4, 6, 2],
        gem10_lpt50 = gim[1, 5, 6, 2],
        gem10_lpt100 = gim[1, 6, 6, 2],
        gem10_lpt250 = gim[1, 7, 6, 2],
        gem10_palbo25 = gim[1, 6, 4, 9],
        gem10_palbo50 = gim[1, 6, 5, 9],
        gem10_palbo100 = gim[1, 6, 6, 9],
        gem10_palbo250 = gim[1, 6, 7, 9],
        lap100_palbo25 = gim[1, 6, 4, 4],
        lap100_palbo50 = gim[1, 6, 5, 4],
        lap100_palbo100 = gim[1, 6, 6, 4],
        lap100_palbo250 = gim[1, 6, 7, 4],
        lap100_gem5 = gim[1, 6, 5, 2],
        lap100_gem10 = gim[1, 6, 6, 2],
        lap100_gem30 = gim[1, 6, 7, 2],
        lap100_gem100 = gim[1, 6, 8, 2],
        pax2_lpt25 = gim[1, 4, 4, 3],
        pax2_lpt50 = gim[1, 5, 4, 3],
        pax2_lpt100 = gim[1, 6, 4, 3],
        pax2_lpt250 = gim[1, 7, 4, 3],
    )
    dfG2 = DataFrames.DataFrame(
        palbo50_lpt25 = gim[2, 4, 5, 4],
        palbo50_lpt50 = gim[2, 5, 5, 4],
        palbo50_lpt100 = gim[2, 6, 5, 4],
        palbo50_lpt250 = gim[2, 7, 5, 4],
        palbo50_gem5 = gim[2, 5, 5, 9],
        palbo50_gem10 = gim[2, 6, 5, 9],
        palbo50_gem30 = gim[2, 7, 5, 9],
        palbo50_gem100 = gim[2, 8, 5, 9],
        gem10_lpt25 = gim[2, 4, 6, 2],
        gem10_lpt50 = gim[2, 5, 6, 2],
        gem10_lpt100 = gim[2, 6, 6, 2],
        gem10_lpt250 = gim[2, 7, 6, 2],
        gem10_palbo25 = gim[2, 6, 4, 9],
        gem10_palbo50 = gim[2, 6, 5, 9],
        gem10_palbo100 = gim[2, 6, 6, 9],
        gem10_palbo250 = gim[2, 6, 7, 9],
        lap100_palbo25 = gim[2, 6, 4, 4],
        lap100_palbo50 = gim[2, 6, 5, 4],
        lap100_palbo100 = gim[2, 6, 6, 4],
        lap100_palbo250 = gim[2, 6, 7, 4],
        lap100_gem5 = gim[2, 6, 5, 2],
        lap100_gem10 = gim[2, 6, 6, 2],
        lap100_gem30 = gim[2, 6, 7, 2],
        lap100_gem100 = gim[2, 6, 8, 2],
        pax2_lpt25 = gim[2, 4, 4, 3],
        pax2_lpt50 = gim[2, 5, 4, 3],
        pax2_lpt100 = gim[2, 6, 4, 3],
        pax2_lpt250 = gim[2, 7, 4, 3],
    )

    XLSX.writetable("G1CombinationDurations.xlsx", dfG1)
    XLSX.writetable("G2CombinationDurations.xlsx", dfG2)
end
