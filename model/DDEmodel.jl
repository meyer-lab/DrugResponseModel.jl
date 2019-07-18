using LeastSquaresOptim, DifferentialEquations, DelayDiffEq, DiffEqBase, Optim, Plots, Statistics, DataFrames, CSV, Distributed

"""
        This file contains functions to fit the data to a Delay Differential Equation model, and find the parameters.
"""

function DDEmodel(du, u, h, p, t)
    du[1] = -p[1]*(h(p, t-p[3])[1]) + 2*p[2]*(h(p, t-p[4])[2]) - p[5]*u[1]
    du[2] = p[1]*(h(p, t-p[3])[1]) - p[2]*(h(p, t-p[4])[2]) - p[6]*u[2]
end

function DDEsolve(pp::Array, i::Int, g1_0::Array, g2_0::Array)
    lags = [pp[3], pp[4]]
    t = LinRange(0.0, 95.5, 192)
    h(pp, t) = 10*ones(2)
    tspan = (0.0, 95.5)
    u0 = [g1_0[i], g2_0[i]]
    prob = DDEProblem(DDEmodel, u0, h, tspan, pp; constant_lags = lags)
    solve(prob, MethodOfSteps(Tsit5()))
end

function resid(pp::Array, i::Int, g1::Matrix, g2::Matrix)
    t = LinRange(0.0, 95.5, 192)
    res = zeros(2, 192)
    sol = DDEsolve(pp, i, g1_0, g2_0)
    res[1, :] = sol(t, idxs=1).u - g1[:, i]
    res[2, :] = sol(t, idxs=2).u - g2[:, i]
    return res
end

function optimIt(initial_guess::Array, lower_bound::Array, upper_bound::Array, i::Int, g1::Matrix, g2::Matrix)
    residuals(pp) = resid(pp, i, g1, g2)
    results_dde = optimize(residuals, initial_guess, LevenbergMarquardt(), lower = lower_bound, upper = upper_bound)
    return results_dde.minimizer
end
