import CSV, DataFrames, Distributed, Statistics
"""
        Imports data, works for both ODE and DDE model
"""

function get_data(path_g2::String, path_total::String)
    # Import data all the trials for each drug
    data = CSV.read(path_g2)
    total = CSV.read(path_total)

    # delete the extra index column

    select(data, 1)
    select(data, 2)
    select(total, 1)
    select(total, 2)


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
    return pop, g2, g1, g2_0, g1_0
end

function remove_peaks(data)
    data = copy(data)

    for i in 1:length(data) - 3
        if (abs(data[i+2] - data[i+1]) > 20*abs(data[i+3] - data[i+2])) && (abs(data[i+1] - data[i]) > 20*abs(data[i+3] - data[i+2]))
            data[i+1] = (data[i] + data[i+2])/2
        end
    end

    return data
end

