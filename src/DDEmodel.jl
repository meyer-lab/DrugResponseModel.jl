using DelayDiffEq, DiffEqParamEstim, BlackBoxOptim, Plots, LsqFit, Optim, Printf

"""
        This file contains functions to fit the data to a Delay Differential Equation model, and find the parameters
"""

# DDE model 
function DDEmodel(du, u, h, p, t)
    du[1] = -p[1]*(h(p, t-p[3])[1]) + 2*p[2]*(h(p, t-p[4])[2]) - p[5]*u[1]
    du[2] = p[1]*(h(p, t-p[3])[1]) - p[2]*(h(p, t-p[4])[2]) - p[6]*u[2]
end

# estimate history function
exp_model(t, p) = @. p[1]*exp(t*p[2]) # exponential model

function find_history(g1::Matrix{Float64}, g2::Matrix{Float64})
    control_g1 = g1[:, 1] # y1
    control_g2= g2[:, 1] # y2
    p0 = [1.0, 0.5]
    times = range(0.0; stop = 95.5, length = 192) # x
    fit_g1 = curve_fit(exp_model, times, control_g1, p0) 
    fit_g2 = curve_fit(exp_model, times, control_g2, p0) 

    return coef(fit_g1), coef(fit_g2)
end

# define problem generator (optimization in log space)
function prob_generator(prob, p)
    exp_p = exp.(p)
    remake(prob; p = exp_p, constant_lags = [exp_p[3], exp_p[4]])
end


function ddesolve(times::Vector{Float64}, g1::Matrix{Float64}, g2::Matrix{Float64}, g1_0::Vector{Float64}, g2_0::Vector{Float64}, params::Vector{Float64}, j::Int)
    
    # history function
    fit1, fit2 = find_history(g1, g2)
    h(p, t) = [exp_model(t, fit1); exp_model(t, fit2)]
    data = vcat(g1[:, j]', g2[:, j]')

    # problem
    prob = DDEProblem(DDEmodel, [g1_0[j], g2_0[j]], h, extrema(times), params;
                      constant_lags = [params[3], params[4]])
    # algorithm to solve
    alg = MethodOfSteps(AutoTsit5(Rosenbrock23()); constrained=true)

    # returning estimated parameteres and the objective function
    return alg, prob, data
end


function optimization(g1::Matrix{Float64}, g2::Matrix{Float64}, g1_0::Vector{Float64}, g2_0::Vector{Float64}, initial_guess::Vector{Float64}, j::Int, lower::Vector{Float64}, upper::Vector{Float64}, num_steps::Int)
    times = range(0.0; stop = 95.5, length = 192)

    alg, prob, data = ddesolve(collect(times), g1, g2, g1_0, g2_0, initial_guess, j)

    # objective function
    obj = build_loss_objective(prob, alg, L2Loss(times, data);
                               prob_generator=prob_generator,
                               verbose_opt = false)
    # optimizing
    results_dde = bboptimize(obj; SearchRange=bound,

                                    NumDimensions=6, TraceMode=:silent, MasSteps=num_steps, Method=:adaptive_de_rand_1_bin_radiuslimited)
    new_guess = best_candidate(results_dde)
    return best_fitness(results_dde), exp.(new_guess)
end