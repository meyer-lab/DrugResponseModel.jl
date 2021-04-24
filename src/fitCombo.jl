""" In this file the data will be organized for combinations to be fit to the model. """

########### 1. control, Palbociclib 50 nM + lapatinibs [25 nM, 50 nM, 100 nM, 250 nM]
########### 2. control, Palbociclib 50 nM + gemcitabines [5 nM, 10 nM, 17 nM, 30 nM]
########### 3. control, Gemcitabine 10 nM + palbociclibs [25 nM, 100 nM, 250 nM, 250 nM]
########### 4. control, Gemcitabine 10 nM + lapatinibs [25 nM, 50 nM, 100 nM, 250 nM]
########### 5. control, Doxorubicin 20 nM + gemcitabines [5 nM, 10 nM, 17 nM, 30 nM]
########### 6. control, Lapatinib  100 nM + palbociclibs [25 nM, 100 nM, 250 nM, 250 nM]

GC = JLD.load("GC.jld")["GC"] # 2 x 193 x 6 x 6 which is G1/G2, timePoints, concentrations, trials

function getODEparamsComboFit(p, conc)
    effects = zeros(eltype(p), 12, 6, 6)
    j=85
    # Scaled drug effect
    for i = 1:6
        xx = 1.0 ./ (1.0 .+ (p[k] ./ conc[:, i]) .^ p[k + 1])

        # [EC50, left, right, steepness]
        effects[1, :, i] = p[j] .+ (p[k + 2] - p[j]) .* xx # a1
        effects[2, :, i] = p[j + 1] .+ (p[k + 3] - p[j + 1]) .* xx # a2
        effects[3, :, i] = p[j + 2] .+ (p[k + 4] - p[j + 2]) .* xx # b1
        effects[4, :, i] = p[j + 3] .+ (p[k + 5] - p[j + 3]) .* xx # b2
        effects[5, :, i] = p[j + 4] .+ (p[k + 6] - p[j + 4]) .* xx # b3
        effects[6, :, i] = p[j + 5] .+ (p[k + 7] - p[j + 5]) .* xx # b4
        effects[7, :, i] = p[k + 8] .* xx   # g11
        effects[8, :, i] = p[k + 9] .* xx   # g12
        effects[9, :, i] = p[k + 10] .* xx  # g21
        effects[10, :, i] = p[k + 11] .* xx # g22
        effects[11, :, i] = p[k + 12] .* xx # g23
        effects[12, :, i] = p[k + 13] .* xx # g24

        k += 14
    end
    return effects
end
