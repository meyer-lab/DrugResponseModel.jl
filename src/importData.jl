"""
        Imports data works for both ODE and DDE model
"""
function get_data(path_g2::String, path_total::String)
    # Import data all the trials for each drug
    data = CSV.read(path_g2)
    total = CSV.read(path_total)

    # delete the extra index column
    select!(data, Not(1:2))
    select!(total, Not(1:2))


    # getting all the 8 trials
    drug = data[1:192, 1:8]
    pop = total[1:192, 1:8]

    # rescaling the experimental data assuming we have 20 initial cells for each trial
    g1 = zeros(size(drug, 1), 8)
    g2 = zeros(size(drug, 1), 8)
    g1_0 = zeros(8)
    g2_0 = zeros(8)

    init_cells = 20.0

    # Unifying the dataset to be all in the unit of [# of cells] at each time point forall the trials for a drug
    for i = 1:8
        pop[:, i] = init_cells * pop[:, i]
        g2[:, i] = 0.01 * pop[:, i] .* drug[:, i]
        g1[:, i] = pop[:, i] .- g2[:, i]
        g2_0[i] = init_cells * (drug[1, i] / 100.0)
        g1_0[i] = init_cells * (1 - drug[1, i] / 100.0)
    end
    # removing the peaks
    for i = 1:8
        pop[:, i] = remove_peaks(pop[:, i])
        g2[:, i] = remove_peaks(g2[:, i])
        g1[:, i] = remove_peaks(g1[:, i])
    end
    return pop, g2, g1, g2_0, g1_0
end

""" Removes noise peaks from the raw data. """
function remove_peaks(data)
    data = copy(data)

    for i = 1:(length(data) - 3)
        if (abs(data[i + 2] - data[i + 1]) > 20 * abs(data[i + 3] - data[i + 2])) &&
           (abs(data[i + 1] - data[i]) > 20 * abs(data[i + 3] - data[i + 2]))
            data[i + 1] = (data[i] + data[i + 2]) / 2
        end
    end

    return data
end

""" This function takes in the drug name which is a string and must be among this list: ["lapatinib", "doxorubicin", "paclitaxel", "gemcitabine"]. It returns the cnocentrations, population, cell, and initial cell number for that drug."""
function setup_data(drug_name::String)
    basePath = joinpath(dirname(pathof(DrugResponseModel)), "..", "data")

    if drug_name == "lapatinib"
        dfname, dfname2, idx = "lap.csv", "lap_pop.csv", 1
    elseif drug_name == "doxorubicin"
        dfname, dfname2, idx = "dox.csv", "dox_pop.csv", 2
    elseif drug_name == "gemcitabine"
        dfname, dfname2, idx = "gem.csv", "gem_pop.csv", 3
    elseif drug_name == "paclitaxel"
        dfname, dfname2, idx = "taxol1.csv", "taxol1_pop.csv", 4
    else
        error("The drug is not amongst the data, please check the drug_name.")
    end

    #----------- import concentrations
    concentration = CSV.read(joinpath(basePath, "concentrations.csv"))
    conc_l = [Float64(concentration[idx, col]) for col = 2:9]

    #------------ import cell data
    pop_l, g2_l, g1_l, g2_0_l, g1_0_l = get_data(joinpath(basePath, dfname), joinpath(basePath, dfname2))

    return conc_l, pop_l, g2_l, g1_l, g2_0_l, g1_0_l
end
