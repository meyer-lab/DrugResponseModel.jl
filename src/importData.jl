"""
        Imports data works for both ODE and DDE model
"""
function get_data(path_g2::String, path_total::String; max = 189)
    # Import data all the trials for each drug
    data = CSV.read(path_g2)
    total = CSV.read(path_total)

    # delete the extra index column
    data = data[:, 2:9]
    total = total[:, 2:9]

    # getting all the 8 trials
    drug = data[1:max, 1:8]
    pop = total[1:max, 1:8]

    # rescaling the experimental data assuming we have 20 initial cells for each trial
    g1 = zeros(size(drug, 1), 8)
    g2 = zeros(size(drug, 1), 8)
    g1_0 = zeros(8)
    g2_0 = zeros(8)
    population = zeros(size(drug, 1), 8)

    init_cells = 20.0

    # Unifying the dataset to be all in the unit of [# of cells] at each time point forall the trials for a drug
    for i = 1:8
        population[:, i] = init_cells * pop[:, i]
        g2[:, i] = 0.01 * population[:, i] .* drug[:, i]
        g1[:, i] = population[:, i] .- g2[:, i]
        g2_0[i] = init_cells * (drug[1, i] / 100.0)
        g1_0[i] = init_cells * (1 - drug[1, i] / 100.0)
    end
    # removing the peaks
    for i = 1:8
        population[:, i] = savitzky_golay_filter(population[:, i], 41, 3)
        g2[:, i] = savitzky_golay_filter(g2[:, i], 41, 3)
        g1[:, i] = savitzky_golay_filter(g1[:, i], 41, 3)
    end
    return population, g2, g1
end

""" This function takes in the drug name which is a string and must be among this list: ["lapatinib", "doxorubicin", "paclitaxel", "gemcitabine"]. It returns the cnocentrations, population, cell, and initial cell number for that drug."""
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
    concentration = CSV.read(joinpath(basePath, "concentrations.csv"))
    conc_l = [Float64(concentration[idx, col]) for col = 2:9]
    conc_l[1] = 0.05

    #------------ import cell data
    pop_l, g2_l, g1_l = get_data(joinpath(basePath, dfname), joinpath(basePath, dfname2))

    return conc_l, pop_l, g2_l, g1_l
end

function load(max, i)
    # i is the replicate number
    concl, popl, g2l, g1l = setup_data(string("Lapatinib", i))
    concd, popd, g2d, g1d = setup_data(string("Doxorubicin", i))
    concg, popg, g2g, g1g = setup_data(string("Gemcitabine", i))
    concp, popp, g2p, g1p = setup_data(string("Paclitaxel", i))
    concpal, poppal, g2pal, g1pal = setup_data(string("Palbociclib", i))
    concentrations = hcat(concl, concd, concg, concp, concpal)

    #     populations = [popl, popd, popg, popp, poppal]
    g1s = zeros(max, 8, 5)
    g2s = zeros(max, 8, 5)
    pops = zeros(max, 8, 5)
    g1s[:, :, 1] = g1l
    g1s[:, :, 2] = g1d
    g1s[:, :, 3] = g1g
    g1s[:, :, 4] = g1p
    g1s[:, :, 5] = g1pal
    g2s[:, :, 1] = g2l
    g2s[:, :, 2] = g2d
    g2s[:, :, 3] = g2g
    g2s[:, :, 4] = g2p
    g2s[:, :, 5] = g2pal
    pops[:, :, 1] = popl
    pops[:, :, 2] = popd
    pops[:, :, 3] = popg
    pops[:, :, 4] = popp
    pops[:, :, 5] = poppal

    return concentrations, pops, g1s, g2s
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
