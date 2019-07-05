"""This file contains functions to:
1. Import data and do some pre-processing 
2. Fit the data to a Delay Differential Equation model, and find the parameters 
3. Fit the parameteres to a Hill-function and find Hill-curve parameters
"""
import CSV
using LeastSquaresOptim, DifferentialEquations, DelayDiffEq, DiffEqBase, Optim, Plots, Statistics, DataFrames


function get_data(path_g2, path_total)
    # Import data all the trials for each drug
    data = CSV.read(path_g2)
    total = CSV.read(path_total)

    # delete the extra index column
    deletecols!(data, 1:2)
    deletecols!(total, 1:2)

    # getting all the 8 trials
    drug = data[:, 1:8]
    pop = total[:, 1:8]

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

# Specifying which trial of the data is being used; could be an integer belonging to {1, 3, ..., 8}
# 1 refers to Control conditions and 8 refers to maximum drug concentration used for treatment.
i = 8;

function DDEmodel(du, u, h, p, t)
    du[1] = -p[1]*(h(p, t-p[3])[1]) + 2*p[2]*(h(p, t-p[4])[2]) - p[6]*u[1]
    du[2] = p[1]*(h(p, t-p[3])[1]) - p[2]*(h(p, t-p[4])[2]) - p[7]*u[2]
end

function DDEsolve(pp)
    lags = [pp[3], pp[4]]
    h(p, t) = pp[5]*ones(2)
    t = LinRange(0.0, 95.5, 192)
    tspan = (0.0, 95.5)
    u0 = [gem1_0[i], gem2_0[i]]
    prob = DDEProblem(DDEmodel, u0, h, tspan, pp; constant_lags = lags)
    solve(prob)
end

function resid(pp)
    t = LinRange(0.0, 95.5, 192)
    res = zeros(2, 192)
    sol = DDEsolve(pp, i)
    res[1, :] = sol(t, idxs=1).u - gem1[:, i]
    res[2, :] = sol(t, idxs=2).u - gem2[:, i]
    return res
end

function optimIt(p, low, upp)
    results_dde = optimize(resid, p, LevenbergMarquardt(), lower = low, upper = upp)
    return results_dde.minimizer
end


# --------------------------- Hill-function fitting ----------------------- #
