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
