import CSV, DataFrames
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
    for i in 1:8
        pop[:, i] = init_cells*pop[:, i]
        g2[:, i] = 0.01*pop[:, i] .* drug[:, i]
        g1[:, i] = pop[:, i] .- g2[:, i]
        g2_0[i] = init_cells*(drug[1, i]/100.0)
        g1_0[i] = init_cells*(1 - drug[1, i]/100.0)
    end
    # removing the peaks
    for i in 1:8
        pop[:, i] = remove_peaks(pop[:, i])
        g2[:, i] = remove_peaks(g2[:, i])
        g1[:, i] = remove_peaks(g1[:, i])
    end
    return pop, g2, g1, g2_0, g1_0
end

# to remove the peaks from the raw data
function remove_peaks(data)
    data = copy(data)

    for i in 1:length(data) - 3
        if (abs(data[i+2] - data[i+1]) > 20*abs(data[i+3] - data[i+2])) && (abs(data[i+1] - data[i]) > 20*abs(data[i+3] - data[i+2]))
            data[i+1] = (data[i] + data[i+2])/2
        end
    end

    return data
end

function setup_data(drug_name::String)
    """ This function takes in the drug name which is a string and must be among this list: ["lapatinib", "doxorubicin", "paclitaxel", "gemcitabine"]. It returns the cnocentrations, population, cell, and initial cell number for that drug."""

    #----------- import concentrations
    concentration = CSV.read(joinpath("..", "data", "concentrations.csv"))
    # lapatinib
    conc_l = [Float64(concentration[1, col]) for col in 2:9]
    # doxorubicin
    conc_d = [Float64(concentration[2, col]) for col in 2:9]
    # gemcitabine
    conc_g = [Float64(concentration[3, col]) for col in 2:9]
    # paclitaxel
    conc_t = [Float64(concentration[4, col]) for col in 2:9]

    #------------ import cell data
    # lapatinib
    pop_l, g2_l, g1_l, g2_0_l, g1_0_l = get_data(joinpath("..", "data", "lap.csv"),
                                       joinpath("..", "data", "lap_pop.csv"));
    # doxorubicin
    pop_d, g2_d, g1_d, g2_0_d, g1_0_d = get_data(joinpath("..", "data", "dox.csv"),
                                       joinpath("..", "data", "dox_pop.csv"));
    # gemcitabine
    pop_g, g2_g, g1_g, g2_0_g, g1_0_g = get_data(joinpath("..", "data", "gem.csv"),
                                       joinpath("..", "data", "gem_pop.csv"));
    # paclitaxel
    pop_t, g2_t, g1_t, g2_0_t, g1_0_t = get_data(joinpath("..", "data", "taxol1.csv"),
                                       joinpath("..", "data", "taxol1_pop.csv"));
    if drug_name == "lapatinib"
        return conc_l, pop_l, g2_l, g1_l, g2_0_l, g1_0_l
    elseif drug_name == "doxorubicin"
        return conc_d, pop_d, g2_d, g1_d, g2_0_d, g1_0_d
    elseif drug_name == "gemcitabine"
        return conc_g, pop_g, g2_g, g1_g, g2_0_g, g1_0_g
    elseif drug_name == "paclitaxel"
        return conc_t, pop_t, g2_t, g1_t, g2_0_t, g1_0_t
    else
        error("The drug is not amongst the data, please check the drug_name.")
    end
end