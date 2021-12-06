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
    basePath = joinpath(dirname(pathof(DrugResponseModel)), "..", "data")
    data = CSV.read(joinpath(basePath, "AU565_round2and3_model_subset_level2.csv"), DataFrame)

    # add a new column of "_" to stick to the end of "treatment" column so we can reshape without loosing columns with the same name
    new_col = similar(data[:, 5])
    for i=1:length(new_col)
        new_col[i] = "_"
    end
    data_ = hcat(data, new_col)

    # transform so the columns are drug conditions
    tmp_norm = transform(data_, [:treatment, :x1, :field, :well] => ByRow(string) => :treat)
    norm_T0 = unstack(tmp_norm, :elapsed_minutes, :treat, :cell_count_norm_T0, allowduplicates=true)
    select!(norm_T0, Not(:elapsed_minutes))

    # remove the time 0:120 to make all conditions the same length and syncronize
    norm_T0 = norm_T0[Not(190:end), :]

    G1_prop = unstack(tmp_norm, :elapsed_minutes, :treat, :G1_proportion, allowduplicates=true)
    G1_prop = G1_prop[Not(190:end), :]
    select!(G1_prop, Not(:elapsed_minutes)) # delete the time column

    norm_t = Matrix{Float64}(norm_T0)
    g_prop = Matrix{Float64}(G1_prop)

    # filter
    for i = 1:size(norm_t, 2)
        norm_t[:, i] = savitzky_golay_filter(norm_t[:, i], 41, 3)
        g_prop[:, i] = savitzky_golay_filter(g_prop[:, i], 41, 3)
    end

    # separate G1 and SG2 cell#
    gs = zeros(2, size(norm_t, 1), size(norm_t, 2)) # G1# , G2# x time x condition
    population = init_cells * norm_t
    gs[1, :, :] = population .* g_prop # G1
    gs[2, :, :] = population - gs[1, :, :] # G2

    condition = similar(names(norm_T0))
    for (index, item) in enumerate(names(norm_T0))
        condition[index] = rsplit(item, "_", limit=2)[1]
    end
    return gs, condition
end

""" Extract 4 replicates of each condition. """
function trim_data(g, c)

    # find unique conditions
    uniq_c = unique(c)
    # remove "control_0" and "vehicle_0" from the conditions
    filter!(e -> e != "vehicle_0", uniq_c)
    filter!(e -> e != "control_0", uniq_c)
    # append the index of 4 replicates of the same condition
    inds = []
    for item in uniq_c
        tm = findall(x -> x == item, c)
        push!(inds, tm[1:4])
    end

    # create a new G with a bit more organized conditions based on unique 
    new_g = zeros(2, 4, 189, length(inds)) # G1/G2 x 4 replicates x 189 data points x 87 conditions
    for i in 1:length(uniq_c)
        new_g[1, :, :, i] = g[1, :, inds[i]]' # G1
        new_g[2, :, :, i] = g[2, :, inds[i]]' # SG2
    end
    @assert size(new_g, 4) == size(uniq_c)[1] == size(inds)[1]

    return new_g, uniq_c
end

drugs = ["5FU", "AZD5438", "Panobinostat", "MG132", "BEZ235", "Everolimus", "JQ1", "Bortezomib", "MK1775", "Trametinib", "Cabozantinib"]

""" Create a tensor form of the data """
function form_tensor(new_g, uniq_c)
    new_ind = []
    for (ind, drug) in enumerate(drugs)
        tm2 = findall( y -> occursin(drug, y), uniq_c)
        push!(new_ind, tm2[1:7])
    end

    tensor = zeros(2, 4, 189, 7, 11)
    conditions = []
    for i = 1:11
        tensor[1, :, :, :, i] = new_g[1, :, :, new_ind[i]]
        tensor[2, :, :, :, i] = new_g[2, :, :, new_ind[i]]
        push!(conditions, uniq_c[new_ind[i]])
    end
    return tensor, conditions
end

""" create one csv file for each drug. """
function output_drugs(g, c)
    uniq_c = unique(c)
    vehicle_index = findall(x -> occursin("vehicle", x), c)
    control_index = findall(x -> occursin("control", x), c)
    # the index of drugs
    new_ind = []
    for (ind, drug) in enumerate(drugs)
        tm2 = findall(y -> occursin(drug, y), uniq_c)
        push!(new_ind, tm2)
    end

    # the indexes 
    inds = []
    for item in uniq_c
        tm = findall(x -> x == item, c)
        push!(inds, tm)
    end

    for i in 1:length(uniq_c)
        new_g[1, :, :, i] = g[1, :, inds[i]]' # G1
        new_g[2, :, :, i] = g[2, :, inds[i]]' # SG2
    end

    for indx, drug in enumerate(drugs)
        df = [DataFrames.DataFrame(g[1, :, new_ind[indx]]) for 1:length(inds)]
        rename!(df, uniq_c[new_ind[indx]])
    end
end

""" Easy way of loading the data, either in the form of a tensor, or each drug, separately. """
function load_data(tensor=true, drug_name=false, average=true)
    g, c = import_data()
    newg, newc = trim_data(g, c)
    if tensor
        return form_tensor(new_g, newc)
        # TODO
    elseif drug_name
        return
    end
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
