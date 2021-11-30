"""
        Imports data works for the ODE model
"""

const init_cells = 1.0

function get_data(path_g2::String, path_total::String; max = 189)
    # Import data all the trials for each drug
    perc = readdlm(path_g2, ','; skipstart = 1)
    total = readdlm(path_total, ','; skipstart = 1)

    # Clip to data of interest
    perc = convert(Array, perc[1:max, 2:9])
    total = convert(Array, total[1:max, 2:9])

    # rescaling the experimental data assuming we have 1 initial cells for each trial
    gs = zeros(2, size(perc, 1), 8)
    population = init_cells * total
    gs[2, :, :] = 0.01 * population .* perc
    gs[1, :, :] = population - gs[2, :, :]

    # removing the peaks
    for i = 1:8
        gs[1, :, i] = savitzky_golay_filter(gs[1, :, i], 41, 3)
        gs[2, :, i] = savitzky_golay_filter(gs[2, :, i], 41, 3)
    end
    return gs
end

function import_combination(filename::String)
    basePath = joinpath(dirname(pathof(DrugResponseModel)), "..", "data")
    path_g2 = string("/", basePath, "/", filename, "_CYCLE.csv")
    path_total = string("/", basePath, "/", filename, "_CELL.csv")
    perc = readdlm(path_g2, ','; skipstart = 1)
    total = readdlm(path_total, ','; skipstart = 1)

    # Clip to data of interest
    perc = convert(Array{Float64, 2}, perc)
    total = convert(Array{Float64, 2}, total)

    # removing the peaks
    for i = 1:size(perc, 2)
        perc[:, i] = savitzky_golay_filter(perc[:, i], 41, 3)
        total[:, i] = savitzky_golay_filter(total[:, i], 41, 3)
    end

    # rescaling the experimental data assuming we have 1 initial cells for each trial
    gs = zeros(3, size(perc, 1), size(perc, 2))
    population = init_cells * total
    gs[1, :, :] = 0.01 * population .* perc
    gs[2, :, :] = population - gs[1, :, :]
    gs[3, :, :] = gs[1, :, :] + gs[2, :, :]

    return gs[:, :, 2:25], gs[:, :, 27:50], perc, total
end

""" This function imports the new round of data and filters. """
function import_data()
    # 96 conditions(not separating replicates) with 191 data points,
    # and 1158 conditions(not separating conditions) with 192 data points.
    # 18336 / 191 = 96
    # 222336 / 193 = 1152
    # so we have a total of 1254 conditions
    basePath = joinpath(dirname(pathof(DrugResponseModel)), "..", "data")

    d = zeros(193, 1248, 2) # time points x conditions x G1% and total

    dataname = "AU565_round2and3_model_subset_level2.csv"
    data = readdlm(joinpath(basePath, dataname), ',')

    conditions = Vector{String}()

    for i = 1:96 # first set that have 191 data points
        d[1:191, i, 1] = data[(191 *(i-1)) + 2 : (191 * i) + 1, 9]
        d[1:191, i, 2] = data[(191 *(i-1)) + 2 : (191 * i) + 1, 8]
        push!(conditions, data[(191 *(i-1)) + 2, 5])
    end
    for i=96:1247 # the rest that have 193 data points
        d[:, i+1, 1] = data[(193 *(i-1)) + 3 : (193 * i) + 2, 9]
        d[:, i+1, 2] = data[(193 *(i-1)) + 3 : (193 * i) + 2, 8]
        push!(conditions, data[(193 *(i-1)) + 3, 5])
    end

    # filter
    for i = 1:size(d, 2)
        d[:, i, 1] = savitzky_golay_filter(d[:, i, 1], 41, 3)
        d[:, i, 2] = savitzky_golay_filter(d[:, i, 2], 41, 3)
    end
    # separate G1 and SG2 cell#
    gs = zeros(2, size(d, 1), size(d, 2)) # G1# , G2# x time x condition
    population = init_cells * d[:, :, 2]
    gs[1, :, :] = population .* d[:, :, 1] # G1
    gs[2, :, :] = population - gs[1, :, :] # G2

    return gs, conditions
end

function reorganize(gs, conditions)
    new_condition = []; c = 0;
    control = []; id_c = []
    FU = []; f = 0; id_f = []
    pano = []; p = 0; id_p = []
    AZD = []; a = 0; id_az = []
    BEZ = []; bz = 0; id_bz = []
    Bort = []; b = 0; id_b = []
    MG = []; m = 0; id_mg= []
    Ever = []; e = 0; id_e = []
    MK = []; mk = 0; id_mk = []
    JQ1 = []; j = 0; id_j = []
    Tram = []; t = 0; id_t = []
    Cabo = []; cb = 0; id_cb = []
    for (key, val) in enumerate(conditions)
        if occursin("vehicle", val)
            append!(control, gs[:, key])
            push!(id_c, key)
            c += 1
        elseif occursin("control", val)
            append!(control, gs[:, key])
            push!(id_c, key)
            c += 1
        elseif occursin("5FU", val)
            append!(FU, gs[:, key])
            push!(id_f, key)
            f += 1
        elseif occursin("Panobinostat", val)
            append!(pano, gs[:, key])
            push!(id_p, key)
            p += 1
        elseif occursin("AZD", val)
            append!(AZD, gs[:, key])
            push!(id_az, key)
            a += 1
        elseif occursin("BEZ", val)
            append!(BEZ, gs[:, key])
            push!(id_bz, key)
            bz += 1
        elseif occursin("Bort", val)
            append!(Bort, gs[:, key])
            push!(id_b, key)
            b += 1
        elseif occursin("MG", val)
            append!(MG, gs[:, key])
            push!(id_mg, key)
            m += 1
        elseif occursin("Ever", val)
            append!(Ever, gs[:, key])
            push!(id_e, key)
            e += 1
        elseif occursin("MK", val)
            append!(MK, gs[:, key])
            push!(id_mk, key)
            mk += 1
        elseif occursin("JQ", val)
            append!(JQ1, gs[:, key])
            push!(id_j, key)
            j += 1
        elseif occursin("Tram", val)
            append!(Tram, gs[:, key])
            push!(id_t, key)
            t += 1
        elseif occursin("Cabo", val)
            append!(Cabo, gs[:, key])
            push!(id_cb, key)
            cb += 1
        end
    end

    new_c = similar(conditions)
    new_cond = vcat(id_c, id_f, id_p, id_az, id_bz, id_b, id_mg, id_e, id_mk, id_j, id_t, id_cb)
    for i=1:length(conditions)
        new_c[i] = conditions[new_cond[i]]
    end
    control = reshape(control, 193, c); FU = reshape(FU, 193, f); pano = reshape(pano, 193, p); AZD = reshape(AZD, 193, a)
    BEZ = reshape(BEZ, 193, bz); Bort = reshape(Bort, 193, b); MG = reshape(MG, 193, m); Ever = reshape(Ever, 193, e)
    MK = reshape(MK, 193, mk); JQ1 = reshape(JQ1, 193, j); Tram = reshape(Tram, 193, t); Cabo = reshape(Cabo, 193, cb)
    new_gs = hcat(control, FU, pano, AZD, BEZ, Bort, MG, Ever, MK, JQ1, Tram, Cabo)
    @assert(size(gs) == size(new_gs))
    return new_gs, new_c
end

function reorganize_data()
    gs, conditions = import_data()
    new_g1, _ = reorganize(gs[1, :, :], conditions)
    new_g2, new_cond = reorganize(gs[2, :, :], conditions)
    gss = zeros(2, 193, 1248)
    gss[1, :, :] = new_g1
    gss[2, :, :] = new_g2

    return gss, new_cond
end

function setup_data(drug_name::String)
    basePath = joinpath(dirname(pathof(DrugResponseModel)), "..", "data")

    dfname = string(drug_name, ".csv")
    dfname2 = string(drug_name, "_pop.csv")

    if occursin("Lapatinib", drug_name)
        idx = 1
    elseif occursin("Doxorubicin", drug_name)
        idx = 2
    elseif occursin("Gemcitabine", drug_name)
        idx = 3
    elseif occursin("Paclitaxel", drug_name)
        idx = 4
    elseif occursin("Palbociclib", drug_name)
        idx = 1
    end

    #----------- import concentrations
    concentration = readdlm(joinpath(basePath, "concentrations.csv"), ','; skipstart = 1)
    conc_l = [Float64(concentration[idx, col]) for col = 2:9]
    conc_l[1] = 0.0

    #------------ import cell data
    gs = get_data(joinpath(basePath, dfname), joinpath(basePath, dfname2))

    return conc_l, gs[2, :, :], gs[1, :, :]
end

function load(max, repi)
    g1s = zeros(max, 8, 5)
    g2s = zeros(max, 8, 5)
    concentrations = zeros(8, 5)
    drugs = ["Lapatinib", "Doxorubicin", "Gemcitabine", "Paclitaxel", "Palbociclib"]

    for i = 1:5
        concentrations[:, i], g2s[:, :, i], g1s[:, :, i] = setup_data(string(drugs[i], repi))
    end

    return concentrations, g1s + g2s, g1s, g2s
end

function savitzky_golay_filter(y::AbstractVector, window_size::Integer, polynomial_order::Integer)
    # input validity checks
    @assert isodd(window_size) "Window size must be an odd integer, i.e. fitting 2m + 1 points around the current value."
    @assert polynomial_order < window_size "Polynomial order must be less than the window size."

    # window size is 2m + 1 points
    m = (window_size - 1) ÷ 2

    # build the Vandermonde design matrix A. Each row corresponds to a point in the fitting window -m:m
    # and each columns correspond to powers in the range 0:polynomial_order
    fitting_points = (-m):m
    A = Matrix{Float64}(undef, window_size, polynomial_order + 1)
    for i = 1:window_size, j = 1:(polynomial_order + 1)
        A[i, j] = fitting_points[i]^(j - 1)
    end

    # for interpolation we'll want the full pseudo-inverse so we can calculate all the fit values at the edges
    # Ap = y
    C = pinv(A)

    # the filter coefficients are the rows of `C`
    filter_coeffs = C[1, :] * factorial(0)

    # convolve with the filter coefficients with a couple extra steps:
    # 1. because of convolution will reverse coefficients we flip before
    # 2. c = conv(a,b) will return a vector of length(c) = length(a) + length(b) - 1 so we chop off the first and last m points
    smoothed = conv(reverse(filter_coeffs), y)[(m + 1):(end - m)]

    # for interpolation edge handling calculate the full fits
    # if we are just smoothing then we can use the design and coefficient matrix as is
    AC = A * C
    smoothed[1:m] = (AC * y[1:window_size])[1:m]
    smoothed[(end - m + 1):end] = (AC * y[(end - window_size + 1):end])[(end - m + 1):end]

    return smoothed
end
