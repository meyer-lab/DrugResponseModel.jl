@testset "Combination tests where g1 prog. rate increased for one drug and decreases for the other drug and g2 prog. rates where in both drugs the rate is decreasing compared to control." begin
    p1 = ones(16, 8)
    p2 = ones(16, 8)

    p1[9:16, 1] .= 0.0
    p2[9:16, 1] .= 0.0
    # drugA
    p1[1, 1] = 0.6          # g1 prog. rate in control
    p1[1, 2:end] .= 0.8     # g1 prog. rate in all other concentrations
    p1[2, 1] = 1.0          # g2 prog. rate in control
    p1[2, 2:end] .= 0.8     # g2 prog. rate in all other concentrations
    # drugB
    p2[2, 1] = 1.0          # g2 prog. rate in control
    p2[2, 2:end] .= 0.5     # g2 prog. rate in all other concentrations
    p2[1, 1] = 0.6          # g1 prog. rate in control
    p2[1, 2:end] .= 0.3     # g1 prog. rate in all other concentrations
    combination = DrugResponseModel.AllBliss_params(p1, p2)
    @assert(all(combination[1:2, 2:end, 2] .>= 0.3))
    @assert(all(combination[1:2, 2:end, 2] .<= 0.5))
end

@testset "Combination tests from estimated parameters to converting to ODE parameters where for both drugs, rates are decreasing and one reaches to zero." begin
    p1 = hcat(vcat(ones(8), zeros(8)), 0.5 * ones(16))
    p2 = hcat(vcat(ones(8), zeros(8)), 1.5 * ones(16))

    cmb = DrugResponseModel.Bliss_params_unit(p1[:, 2], p2[:, 2], hcat(p1[:, 1], p2[:, 1]))
    @assert(cmb[1:8] == 0.75 * ones(8))
    @assert(cmb[9:16] == 2.0 * ones(8))
end

@testset "Test if the function that calculates bliss for cell numbers, works right." begin
    gt1, _ = DrugResponseModel.import_combination("AU01001") # [3, 193, 24]
    control = gt1[3, 1:50, 1]
    lpt50 = gt1[3, 1:50, 3]
    combin = DrugResponseModel.pair_cellnum_Bliss(hcat(control, lpt50), hcat(control, control))
    @assert(all(combin ≈ lpt50))

    Total = zeros(193, 5, 5) # time x concentrations x 5 drugs
    Total[:, 1, :] .= gt1[3, :, 1] # controls
    Total[:, 2:5, 1] .= gt1[3, :, 2:5] # lapatinibs
    Total[:, 2, 2] .= gt1[3, :, 6] # dox 20 nM
    Total[:, 2:5, 3] .= gt1[3, :, 19:22] # gemcitabines
    Total[:, 2, 4] .= gt1[3, :, 13] # pax 2 nM
    Total[:, 2:5, 5] .= gt1[3, :, 7:10] # palbos
    cellnum = zeros(30, 5, 5, 10) # the first one changes with rows, which is the drug that comes first (e.g., in lapatinib (rows-first) dox (columns-second))
    for i = 1:30
        cellnum[i, :, :, :] .= DrugResponseModel.blissCellNum(Total[i, :, :])
    end
    @assert(all(cellnum[:, 3, 1, 1] ≈ Total[1:30, 3, 1])) # control + lap50 ≈ lap50
    @assert(all(cellnum[:, 1, 1, 2] ≈ Total[1:30, 1, 3])) # control + control ≈ control
    @assert(all(cellnum[:, 1, 1, 2] ≈ Total[1:30, 1, 3])) # control + control ≈ control
    @assert(all(cellnum[:, 1, 4, 2] ≈ Total[1:30, 4, 3])) # control + gem 17 ≈ gem 17
end
