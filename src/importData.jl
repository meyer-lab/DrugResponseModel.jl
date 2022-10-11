"""
Imports data for the ODE model
"""

const init_cells = 1.0
basePath = joinpath(dirname(pathof(DrugResponseModel)), "..", "data")

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

""" This function imports the _new round_ of data and filters. """
function import_data(hc=false)
    if hc == false
        basePath = joinpath(dirname(pathof(DrugResponseModel)), "..", "data")
        data = CSV.read(joinpath(basePath, "AU565_round2and3_model_subset_level2.csv"), DataFrame)
    else
        basePath = joinpath(dirname(pathof(DrugResponseModel)), "..", "data", "HC")
        data = CSV.read(joinpath(basePath, hc), DataFrame)
        data.treatment = replace.(data.treatment, Ref("Untreated" => "Untreated_0"))
    end

    # add a new column of "_" to stick to the end of "treatment" column so we can reshape without loosing columns with the same name
    new_col = similar(data.treatment)
    for i=1:length(new_col)
        new_col[i] = "_"
    end
    data_ = hcat(data, new_col)

    # transform so the columns are drug conditions
    tmp_norm = DataFrames.transform(data_, [:treatment, :x1, :field, :well] => ByRow(string) => :treat)
    norm_T0 = unstack(tmp_norm, :elapsed_minutes, :treat, :cell_count_norm_T0, allowduplicates=true)
    select!(norm_T0, Not(:elapsed_minutes))

    # remove the time 0:120 to make all conditions the same length and syncronize
    if hc == false
        norm_T0 = norm_T0[Not(190:end), :]
    else
        norm_T0 = norm_T0[Not(183:end), :]
    end

    @assert(sum(collect(any(ismissing, c) for c in eachcol(norm_T0))) == 0)
    G1_prop = unstack(tmp_norm, :elapsed_minutes, :treat, :G1_proportion, allowduplicates=true)
    if hc == false
        G1_prop = G1_prop[Not(190:end), :]
    else
        G1_prop = G1_prop[Not(183:end), :]
    end
    select!(G1_prop, Not(:elapsed_minutes)) # delete the time column

    @assert(sum(collect(any(ismissing, c) for c in eachcol(G1_prop))) == 0)
    norm_t = Matrix{Float64}(norm_T0)
    g_prop = Matrix{Float64}(G1_prop)

    # filter
    for i = 1:size(norm_t, 2)
        norm_t[:, i] = savitzky_golay_filter(norm_t[:, i], 61, 3)
        g_prop[:, i] = savitzky_golay_filter(g_prop[:, i], 61, 3)
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
    if hc == false
        return gs, condition
    else
        return gs, condition, unique(data.Drug1)
    end
end

""" Extract 4 replicates of each condition. """
function trim_data(g, c)

    # find unique conditions
    uniq_c = unique(c)
    # remove "control_0" and "vehicle_0" from the conditions
    filter!(e -> e != "vehicle_0", uniq_c)
    filter!(e -> e != "control_0", uniq_c)
    # filter!(e -> e != "Untreated_0", uniq_c)
    # append the index of 4 replicates of the same condition
    inds_lengths = []
    for item in uniq_c
        tm = findall(x -> x == item, c)
        push!(inds_lengths, length(tm))
    end

    l = minimum(inds_lengths) # the GCD of the number of replicates, since not all conditions have 4 replicates...
    inds = []
    for item in uniq_c
        tm = findall(x -> x == item, c)
        push!(inds, tm[1:l])
    end
    # create a new G with a bit more organized conditions based on unique 
    new_g = zeros(2, l, size(g)[2], length(inds)) # G1/G2 x 4 replicates x 189 data points x conditions
    for i in 1:length(uniq_c)
        new_g[1, :, :, i] = g[1, :, inds[i]]' # G1
        new_g[2, :, :, i] = g[2, :, inds[i]]' # SG2
    end
    @assert size(new_g, 4) == size(uniq_c)[1] == size(inds)[1]

    return new_g
end

drugs = ["BEZ235", "Trametinib", "5FU", "AZD5438", "Panobinostat", "MG132", "Everolimus", "JQ1", "Bortezomib", "MK1775", "Cabozantinib"]

""" Create a tensor form of the data """
function form_tensor(new_g, c)
    uniq_c = unique(c)
    filter!(e -> e != "vehicle_0", uniq_c)
    filter!(e -> e != "control_0", uniq_c)
    # filter Bortezomib high concentrations, because not the trend of interest
    new_ind = []
    for (ind, drug) in enumerate(drugs)
        tm2 = findall( y -> occursin(drug, y), uniq_c)
        push!(new_ind, tm2[1:7])
    end

    tensor = zeros(2, 4, 189, 8, 11)
    vehicle_indx = findall( y -> occursin("vehicle", y), c)
    for i=1:4
        tensor[:, :, :, 1, :] .= new_g[:, :, :, vehicle_indx[i]]
    end
    conditions = []
    for i = 1:11
        tensor[1, :, :, 2:8, i] = new_g[1, :, :, new_ind[i]]
        tensor[2, :, :, 2:8, i] = new_g[2, :, :, new_ind[i]]
        push!(conditions, pushfirst!(uniq_c[new_ind[i]], "Vehicle"))
    end
    return tensor, conditions
end

""" Form tensor for the HCC1143 dataset. """
function hcc_tensor(gs, new_g, cc, ndrugs)
    uniq_c = unique(cc)
    filter!(e -> e != "Untreated", ndrugs)
    new_ind = []
    for (ind, drug) in enumerate(ndrugs)
        tm2 = findall( y -> occursin(drug, y), uniq_c)
        push!(new_ind, tm2[1:7])
    end

    print(new_ind)
    tensor = zeros(2, size(new_g)[2], size(new_g)[3], 8, 3)
    vehicle_indx = findall( y -> occursin("Untreated", y), cc)

    control = mean(gs[:, :, vehicle_indx], dims=3)[:, :, 1]
    b = zeros(2, size(new_g)[2], size(new_g)[3])
    for i=1:size(new_g)[2]
        b[:, i, :] = control
    end

    tensor[:, :, :, 1, :] .= b
    conditions = []
    for i = 1:3
        tensor[1, :, :, 2:8, i] = new_g[1, :, :, new_ind[i]]
        tensor[2, :, :, 2:8, i] = new_g[2, :, :, new_ind[i]]
        push!(conditions, pushfirst!(uniq_c[new_ind[i]], "Untreated_0"))
    end
    return tensor, conditions, gs[:, :, vehicle_indx]
end

""" This function imports all 11 drugs of the new round of AU565 experiments, with average of replicates. """
function au565_all()
    a, b = DrugResponseModel.import_data()
    newg = DrugResponseModel.trim_data(a, b)
    tens, conds = DrugResponseModel.form_tensor(newg, b)

    names = []
    cons = []
    for cc in conds
        c = []
        for item in cc[2:end]
            push!(c, parse(Float64, rsplit(item, "_")[2]))
        end
        push!(names, rsplit(cc[2], "_")[1])
        push!(cons, c)
    end
    concs = []
    for item in cons
        item = vcat([0.0], item)
        item = convert(Array{Float64, 1}, item)
        push!(concs, item)
    end

    tensor = mean(tens, dims=2)[:, 1, :, :, :]

    return tensor, names, concs, conds
end

""" This function puts together the data from all drug treatments of HCC1143. """
function hcc_all()
    g1, c1, d1 = DrugResponseModel.import_data("HC01301_level_2.csv")
    newg1 = DrugResponseModel.trim_data(g1, c1)
    ten1, cond1, cl1 = DrugResponseModel.hcc_tensor(g1, newg1, c1, d1)
    t1 = mean(ten1, dims=2)

    g2, c2, d2 = DrugResponseModel.import_data("HC01401_level_2.csv")
    newg2 = DrugResponseModel.trim_data(g2, c2)
    ten2, cond2, cl2 = DrugResponseModel.hcc_tensor(g2, newg2, c2, d2)
    t2 = mean(ten2, dims=2)

    conds = append!(cond1, cond2)
    names = []
    concs = []
    for cc in conds
        c = []
        for item in cc
            push!(c, parse(Float64, rsplit(item, "_")[2]))
        end
        push!(names, rsplit(cc[2], "_")[1])
        push!(concs, c)
    end
    
    # g3, c3, d3 = import_data("HC00901_level_2.csv")
    # newg3 = trim_data(g3, c3)
    # ten3, cond3 = hcc_tensor(g3, newg3, c3, d3)

    # g4, c4, d4 = import_data("HC01001_level_2.csv")
    # newg4 = trim_data(g4, c4)
    # ten4, cond4 = hcc_tensor(g4, newg4, c4, d4)

    # return cat(ten3, ten4, dims=4), cond3, cond4
    ten = cat(t1, t2, dims=5)[:, 1, :, :, :]
    return ten[:, :, :, [4, 5, 2, 3]], names[[4, 5, 2, 3]], concs[[4, 5, 2, 3]], conds[[4, 5, 2, 3]]
end

""" create one csv file for each drug. """
function output_drugs(g, c, drugname)
    uniq_c = unique(c)
    filter!(x->x!="vehicle_0", uniq_c)
    filter!(x->x!="control_0", uniq_c)

    inds = []
    for item in uniq_c
        tm = findall(x -> x == item, c)
        push!(inds, tm)
    end

    # create a new G with a bit more organized conditions based on unique
    new_g = zeros(2, 189, 8, length(uniq_c)) # control
    for i in 1:length(uniq_c)
        new_g[:, :, 1:length(inds[i]), i] .= g[:, :, inds[i]] # G1
    end

    # conditions
    new_ind = []
    for (ind, drug) in enumerate(drugs)
        tm2 = findall( y -> occursin(drug, y), uniq_c)
        push!(new_ind, tm2)
    end

    i = findall(y -> y==drugname, drugs)[1]
    gg1 = new_g[1, :, :, new_ind[i]]
    g1 = zeros(size(gg1)[1], size(gg1)[2], size(gg1)[3]+1)
    g1[:, :, 2:1+size(gg1)[3]] .= gg1
    g1[:, :, 1] = g[1, :, findall(y -> occursin("vehicle", y), c)[1:8]]
    gg2 = new_g[2, :, :, new_ind[i]]
    g2 = zeros(size(gg2)[1], size(gg2)[2], size(gg2)[3]+1)
    g2[:, :, 2:1+size(gg2)[3]] .= gg2
    g2[:, :, 1] = g[2, :, findall(y -> occursin("vehicle", y), c)[1:8]]

    return g1, g2
end

function concentration(drugname)
    # import concentrations
    xlx = XLSX.readxlsx(string("/", basePath, "/conditions.xlsx"))
    conds = xlx["Sheet1"][:]
    concentrations = conds[2:end, 2:end]
    # find index of the drug
    drugs = conds[2:end, 1]
    i = findall(y -> y==drugname, drugs)[1]
    return concentrations[i, :]
end


function import_mt1(mt, dr)
    basePath = joinpath(dirname(pathof(DrugResponseModel)), "..", "data", dr)
    data = CSV.read(joinpath(basePath, mt), DataFrame)
    data.treatment = replace.(data.treatment, Ref("Untreated" => "Untreated_0"))

    # add a new column of "_" to stick to the end of "treatment" column so we can reshape without loosing columns with the same name
    new_col = similar(data.treatment)
    for i=1:length(new_col)
        new_col[i] = "_"
    end
    data_ = hcat(data, new_col)

    # transform so the columns are drug conditions
    tmp_norm = DataFrames.transform(data_, [:treatment, :x1, :field, :well] => ByRow(string) => :treat)
    norm_T0 = unstack(tmp_norm, :elapsed_minutes, :treat, :cell_count_norm_T0, allowduplicates=true)
    select!(norm_T0, Not(:elapsed_minutes))
    norm_T0 = DataFrames.transform(norm_T0, names(norm_T0) .=> Impute.locf, renamecols=false)
    # norm_T0 = dropmissing(norm_T0)

    if sum(collect(any(ismissing, c) for c in eachcol(norm_T0))) != 0
        print(mt, sum(collect(any(ismissing, c) for c in eachcol(norm_T0))))
    end
    @assert(sum(collect(any(ismissing, c) for c in eachcol(norm_T0))) == 0)
    G1_prop = unstack(tmp_norm, :elapsed_minutes, :treat, :G1_proportion, allowduplicates=true)
    select!(G1_prop, Not(:elapsed_minutes)) # delete the time column
    G1_prop = DataFrames.transform(G1_prop, names(G1_prop) .=> Impute.locf, renamecols=false)

    # G1_prop = dropmissing(G1_prop)
    @assert(sum(collect(any(ismissing, c) for c in eachcol(G1_prop))) == 0)
    norm_t = Matrix{Float64}(norm_T0)
    g_prop = Matrix{Float64}(G1_prop)

    # filter
    for i = 1:size(norm_t, 2)
        norm_t[:, i] = savitzky_golay_filter(norm_t[:, i], 61, 3)
        g_prop[:, i] = savitzky_golay_filter(g_prop[:, i], 61, 3)
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
    return gs, condition, unique(data.Drug1)
end


""" Form tensor for the 21MT1 dataset. """
function mt_tensor(gs, new_g, cc, ndrugs)
    uniq_c = unique(cc)
    filter!(e -> e != "Untreated", ndrugs)
    new_ind = []
    for (ind, drug) in enumerate(ndrugs)
        tm2 = findall( y -> occursin(drug, y), uniq_c)
        push!(new_ind, tm2[1:7])
    end

    tensor = zeros(2, size(new_g)[2], size(new_g)[3], 8, length(new_ind))
    vehicle_indx = findall( y -> occursin("Untreated", y), cc)

    control = mean(gs[:, :, vehicle_indx], dims=3)[:, :, 1]
    b = zeros(2, size(new_g)[2], size(new_g)[3])
    for i=1:size(new_g)[2]
        b[:, i, :] = control
    end

    tensor[:, :, :, 1, :] .= b
    conditions = []
    for i = 1:length(new_ind)
        tensor[1, :, :, 2:8, i] = new_g[1, :, :, new_ind[i]]
        tensor[2, :, :, 2:8, i] = new_g[2, :, :, new_ind[i]]
        push!(conditions, pushfirst!(uniq_c[new_ind[i]], "Untreated_0"))
    end
    return tensor, conditions, gs[:, :, vehicle_indx]
end


""" This function puts together the data from all drug treatments of 21MT1 data. """
function mt1_all()
    g1, c1, d1 = DrugResponseModel.import_mt1("2101001_CtcK_level_2.csv", "21MT1")
    newg1 = DrugResponseModel.trim_data(g1, c1)
    ten1, cond1, cl1 = DrugResponseModel.mt_tensor(g1, newg1, c1, d1)
    t1 = mean(ten1, dims=2)

    g2, c2, d2 = DrugResponseModel.import_mt1("2101201_CtcK_level_2.csv", "21MT1")
    newg2 = DrugResponseModel.trim_data(g2, c2)
    ten2, cond2, cl2 = DrugResponseModel.mt_tensor(g2, newg2, c2, d2)
    t2 = mean(ten2, dims=2)

    # g3, c3, d3 = DrugResponseModel.import_mt1("2101301_CtcK_level_2.csv", "21MT1")
    # newg3 = DrugResponseModel.trim_data(g3, c3)
    # ten3, cond3, cl3 = DrugResponseModel.mt_tensor(g3, newg3, c3, d3)
    # t3 = mean(ten3, dims=2)

    g4, c4, d4 = DrugResponseModel.import_mt1("2101401_CtcK_level_2.csv", "21MT1")
    newg4 = DrugResponseModel.trim_data(g4, c4)
    ten4, cond4, cl4 = DrugResponseModel.mt_tensor(g4, newg4, c4, d4)
    t4 = mean(ten4, dims=2)

    g5, c5, d5 = DrugResponseModel.import_mt1("2101501_CtcK_level_2.csv", "21MT1")
    newg5 = DrugResponseModel.trim_data(g5, c5)
    ten5, cond5, cl5 = DrugResponseModel.mt_tensor(g5, newg5, c5, d5)
    t5 = mean(ten5, dims=2)

    g6, c6, d6 = DrugResponseModel.import_mt1("2101601_CtcK_level_2.csv", "21MT1")
    newg6 = DrugResponseModel.trim_data(g6, c6)
    ten6, cond6, cl6 = DrugResponseModel.mt_tensor(g6, newg6, c6, d6)
    t6 = mean(ten6, dims=2)

    g7, c7, d7 = DrugResponseModel.import_mt1("2101701_CtcK_level_2.csv", "21MT1")
    newg7 = DrugResponseModel.trim_data(g7, c7)
    ten7, cond7, cl7 = DrugResponseModel.mt_tensor(g7, newg7, c7, d7)
    t7 = mean(ten7, dims=2)

    g8, c8, d8 = DrugResponseModel.import_mt1("2101801_CtcK_level_2.csv", "21MT1")
    newg8 = DrugResponseModel.trim_data(g8, c8)
    ten8, cond8, cl8 = DrugResponseModel.mt_tensor(g8, newg8, c8, d8)
    t8 = mean(ten8, dims=2)

    conds = append!(cond1, cond2, cond4, cond5, cond6, cond7, cond8)
    nams = []
    concs = []
    for cc in conds
        c = []
        for item in cc
            push!(c, parse(Float64, rsplit(item, "_")[2]))
        end
        push!(nams, rsplit(cc[2], "_")[1])
        push!(concs, c)
    end
    
    ten =  cat(t1, t2, t4, t5, t6, t7, t8, dims=5)[:, 1, :, :, :]
    tensor = zeros(2, 193, 8, 6)
    tensor[:, :, 1, :] .= mean(ten[:, :, 1, :], dims=4)[:, :, 1] # untreated
    tensor[:, :, :, 1] = mean(ten[:, :, :, [1, 3, 5, 11, 17]], dims=4) # taxol
    tensor[:, :, :, 2] = mean(ten[:, :, :, [2, 4, 6, 12, 18]], dims=4) # palbo
    tensor[:, :, :, 3] = mean(ten[:, :, :, [7, 13, 19]], dims=4) # trametinib
    tensor[:, :, :, 4] = mean(ten[:, :, :, [8, 14]], dims=4) # bez
    tensor[:, :, :, 5] = mean(ten[:, :, :, [9, 15]], dims=4) # dox
    tensor[:, :, :, 6] = mean(ten[:, :, :, [10, 16]], dims=4) # gem

    css = [concs[1], concs[2], concs[7], concs[8], concs[9], concs[10]]
    conditions = [conds[1], conds[2], conds[7], conds[8], conds[9], conds[10]]
    names = [nams[1], nams[2], nams[7], nams[8], nams[9], nams[10]]

    return tensor[:, :, :, [1, 2, 5, 6]], names[[1, 2, 5, 6]], css[[1, 2, 5, 6]], conditions[[1, 2, 5, 6]]
end

""" This function puts together the data from all drug treatments of MDAMB157 data. """
function mda_all()
    g1, c1, d1 = DrugResponseModel.import_mt1("MD00101_CtcK_level_2.csv", "MDAMB")
    newg1 = DrugResponseModel.trim_data(g1, c1)
    ten1, cond1, cl1 = DrugResponseModel.mt_tensor(g1, newg1, c1, d1)
    t1 = mean(ten1, dims=2)

    g2, c2, d2 = DrugResponseModel.import_mt1("MD00201_CtcK_level_2.csv", "MDAMB")
    newg2 = DrugResponseModel.trim_data(g2, c2)
    ten2, cond2, cl2 = DrugResponseModel.mt_tensor(g2, newg2, c2, d2)
    t2 = mean(ten2, dims=2)

    g3, c3, d3 = DrugResponseModel.import_mt1("MD00301_CtcK_level_2.csv", "MDAMB")
    newg3 = DrugResponseModel.trim_data(g3, c3)
    ten3, cond3, cl3 = DrugResponseModel.mt_tensor(g3, newg3, c3, d3)
    t3 = mean(ten3, dims=2)

    g4, c4, d4 = DrugResponseModel.import_mt1("MD00401_CtcK_level_2.csv", "MDAMB")
    newg4 = DrugResponseModel.trim_data(g4, c4)
    ten4, cond4, cl4 = DrugResponseModel.mt_tensor(g4, newg4, c4, d4)
    t4 = mean(ten4, dims=2)

    g5, c5, d5 = DrugResponseModel.import_mt1("MD00501_CtcK_level_2.csv", "MDAMB")
    newg5 = DrugResponseModel.trim_data(g5, c5)
    ten5, cond5, cl5 = DrugResponseModel.mt_tensor(g5, newg5, c5, d5)
    t5 = mean(ten5, dims=2)

    g6, c6, d6 = DrugResponseModel.import_mt1("MD00601_CtcK_level_2.csv", "MDAMB")
    newg6 = DrugResponseModel.trim_data(g6, c6)
    ten6, cond6, cl6 = DrugResponseModel.mt_tensor(g6, newg6, c6, d6)
    t6 = mean(ten6, dims=2)

    g7, c7, d7 = DrugResponseModel.import_mt1("MD00701_CtcK_level_2.csv", "MDAMB")
    newg7 = DrugResponseModel.trim_data(g7, c7)
    ten7, cond7, cl7 = DrugResponseModel.mt_tensor(g7, newg7, c7, d7)
    t7 = mean(ten7, dims=2)

    g8, c8, d8 = DrugResponseModel.import_mt1("MD00801_CtcK_level_2.csv", "MDAMB")
    newg8 = DrugResponseModel.trim_data(g8, c8)
    ten8, cond8, cl8 = DrugResponseModel.mt_tensor(g8, newg8, c8, d8)
    t8 = mean(ten8, dims=2)

    conds = append!(cond1, cond2, cond3, cond4, cond5, cond6, cond7, cond8)
    nams = []
    concs = []
    for cc in conds
        c = []
        for item in cc
            push!(c, parse(Float64, rsplit(item, "_")[2]))
        end
        push!(nams, rsplit(cc[2], "_")[1])
        push!(concs, c)
    end
    
    ten =  cat(t1, t2, t3, t4, t5, t6, t7, t8, dims=5)[:, 1, :, :, :]
    tensor = zeros(2, 193, 8, 6)
    tensor[:, :, 1, :] .= mean(ten[:, :, 1, :], dims=4)[:, :, 1] # untreated
    tensor[:, :, :, 1] = mean(ten[:, :, :, [7, 13, 19]], dims=4) # taxol
    tensor[:, :, :, 2] = mean(ten[:, :, :, [2, 8, 14, 20]], dims=4) # palbo
    tensor[:, :, :, 3] = mean(ten[:, :, :, [3, 9, 15, 21]], dims=4) # trametinib
    tensor[:, :, :, 4] = mean(ten[:, :, :, [4, 10, 16]], dims=4) # bez
    tensor[:, :, :, 5] = mean(ten[:, :, :, [5, 11, 17]], dims=4) # dox
    tensor[:, :, :, 6] = mean(ten[:, :, :, [1, 6, 12, 18]], dims=4) # gem

    css = [concs[7], concs[2], concs[3], concs[4], concs[5], concs[1]]
    conditions = [conds[7], conds[2], conds[3], conds[4], conds[5], conds[1]]
    names = [nams[7], nams[2], nams[3], nams[4], nams[5], nams[1]]

    return tensor, names, css, conditions
end
