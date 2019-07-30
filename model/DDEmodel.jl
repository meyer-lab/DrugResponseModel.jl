using LeastSquaresOptim, DifferentialEquations, DelayDiffEq, DiffEqBase, Optim, Plots, Statistics, DataFrames, CSV, Distributed, LsqFit

"""
        This file contains functions to fit the data to a Delay Differential Equation model, and find the parameters.
"""

function DDEmodel(du, u, h, p, t)
    du[1] = -p[1]*(h(p, t-p[3])[1]) + 2*p[2]*(h(p, t-p[4])[2]) - p[5]*u[1]
    du[2] = p[1]*(h(p, t-p[3])[1]) - p[2]*(h(p, t-p[4])[2]) - p[6]*u[2]
end

exp_model(t, p) = @. p[1]*exp(t*p[2]) # exponential model

function find_history(g1::Matrix, g2::Matrix)

    control_g1 = g1[:, 1] # y1
    control_g2= g2[:, 1] # y2
    p0 = [1.0, 0.5]
    times = range(0.0; stop = 95.5, length = 192) # x
    fit_g1 = curve_fit(exp_model, times, control_g1, p0) 
    fit_g2 = curve_fit(exp_model, times, control_g2, p0) 

    return coef(fit_g1), coef(fit_g2)
end
    
function DDEsolve(pp::Array, i::Int, g1_0::Array, g2_0::Array)
    lags = [pp[3], pp[4]]
    fit1, fit2 = find_history(g1, g2)
    h(pp, t) = [exp_model(t, coef(fit1)); exp_model(t, coef(fit2))]
    times = range(0.0; stop = 95.5, length = 192) # x
    u0 = [g1_0[i], g2_0[i]]
    alg = MethodOfSteps(AutoTsit5(Rosenbrock23()))
    prob = DDEProblem(DDEmodel, u0, h, extrema(times), pp; constant_lags = lags)
end

function prob_generator(prob, p)
    exp_p = exp.(p)
    remake(prob; p = exp_p, constant_lags = [exp_p[3], exp_p[4]])
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
