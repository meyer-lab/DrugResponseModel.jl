using LeastSquaresOptim, DifferentialEquations, DelayDiffEq, DiffEqBase, Optim, Plots, Statistics, DataFrames, CSV, Distributed

"""
        This file contains functions to:
    1. Import data, extracting G1 and G2 data in # of cells.  
    2. Fit the data to a Delay Differential Equation model, and find the parameters.
"""

function get_data(path_g2::String, path_total::String)
    # Import data all the trials for each drug
    data = CSV.read(path_g2)
    total = CSV.read(path_total)

    # delete the extra index column
    deletecols!(data, 1:2)
    deletecols!(total, 1:2)

    # getting all the 8 trials
    drug = data[1:192, 1:8]
    pop = total[1:192, 1:8]

    # rescaling the experimental data assuming we have 20 initial cells for each trial
    g1 = zeros(size(drug, 1), 8)
    g2 = zeros(size(drug, 1), 8)
    g1_0 = zeros(8)
    g2_0 = zeros(8)
    initial_val = fill([], 8)

    init_cells = 20.0

    # Unifying the dataset to be all in the unit of [# of cells] at each time point forall the trials for a drug
    for i in 1:8
        pop[:, i] = init_cells*pop[:, i]
        g2[:, i] = 0.01*pop[:, i] .* drug[:, i]
        g1[:, i] = pop[:, i] .- g2[:, i]
        g2_0[i] = init_cells*(drug[1, i]/100.0)
        g1_0[i] = init_cells*(1 - drug[1, i]/100.0)
        initial_val[i] = [g1_0[i], g2_0[i]]
    end

    return pop, g2, g1, initial_val
end

# Specifying which trial of the data is being used; could be an integer belonging to {1, 3, ..., 8}
# 1 refers to Control conditions and 8 refers to maximum drug concentration used for treatment.

function DDEmodel(du, u, h, p, t)
    du[1] = -p[1]*(h(p, t-p[3])[1]) + 2*p[2]*(h(p, t-p[4])[2]) - p[6]*u[1]
    du[2] = p[1]*(h(p, t-p[3])[1]) - p[2]*(h(p, t-p[4])[2]) - p[7]*u[2]
end

function DDEsolve(pp, i)
    lags = [pp[3], pp[4]]
    h(pp, t) = pp[5]*ones(2)
    t = LinRange(0.0, 95.5, 192)
    tspan = (0.0, 95.5)
    u0 = initial_val[i]
    prob = DDEProblem(DDEmodel, u0, h, tspan, pp; constant_lags = lags)
    sol = solve(prob, MethodOfSteps(Tsit5()))
    return sol
end

function resid(pp)
    t = LinRange(0.0, 95.5, 192)
    res = zeros(2, 192)
    sol = DDEsolve(pp, i)
    res[1, :] = sol(t, idxs=1).u - g1[:, i]
    res[2, :] = sol(t, idxs=2).u - g2[:, i]
    return res
end

function optimIt(initial_guess, lower_bound, upper_bound)
    results_dde = optimize(resid, initial_guess, LevenbergMarquardt(), lower = lower_bound, upper = upper_bound)
    return results_dde.minimizer
end

# ----------- an example for how to use the functions ----------- #
"""
gem_pop, gem2, gem1, gem2_0, gem1_0 = get_data("..//data//gem.csv", "..//data//gem_pop.csv");
i = 3
p  = [0.0571958, 0.0326492, 0.7504399, 4.501804, 19.3794, 0.00560244, 0.00236748]
# setting lowest delay for tau1 to be half an hour and for tau2 to be 3 hours.
low = [0.0001, 0.0001, 0.5, 3.0, 1.0, 0.0001, 0.0001]
upp = [0.5, 0.5, 30.0, 30.0, 200.0, 0.5, 0.5]
params = optimIt(p, low, upp)
plotIt(params, "Gemcitabine Trial 3")
"""
# --------------------------- Hill-function fitting ----------------------- #
