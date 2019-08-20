using DelayDiffEq, DiffEqParamEstim, BlackBoxOptim, Plots, LsqFit
gr()
"""
        This file contains functions to fit the data to a Delay Differential Equation model, and find the parameters.
"""

# DDE model 
function DDEmodel(du, u, h, p, t)
    du[1] = -p[1]*(h(p, t- 1/p[1])[1]) + 2*p[2]*(h(p, t- 1/p[2])[2]) - p[3]*u[1]
    du[2] = p[1]*(h(p, t- 1/p[1])[1]) - p[2]*(h(p, t- 1/p[2])[2]) - p[4]*u[2]
end

# estimate history function
exp_model(t, p) = @. p[1]*exp(t*p[2]) # exponential model

function find_history(g1, g2)

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
    remake(prob; p = exp_p, constant_lags = [1/exp_p[1], 1/exp_p[2]])
end


function ddesolve(g1, g2, g1_0, g2_0, params, j)
    times = range(0.0; stop = 95.5, length = 192)
    data = vcat(g1[:, j]', g2[:, j]')
    
    # history function
    fit1, fit2 = find_history(g1, g2)
    h(p, t) = [exp_model(t, fit1); exp_model(t, fit2)]

    # problem
    prob = DDEProblem(DDEmodel, [g1_0[j], g2_0[j]], h, extrema(times), params;
                      constant_lags = [1/params[1], 1/params[2]])
    # algorithm to solve
    alg = MethodOfSteps(AutoTsit5(Rosenbrock23()); constrained=true)

    # objective function
    obj = build_loss_objective(prob, alg, L2Loss(times, data);
                               prob_generator=prob_generator,
                               verbose_opt=false)

    # returning estimated parameteres and the objective function
    return obj(params)
end

function optimization(g1, g2, g1_0, g2_0, initial_guess, j, bound, num_steps)
    times = range(0.0; stop = 95.5, length = 192)
    data = vcat(g1[:, j]', g2[:, j]')
    
    # history function
    fit1, fit2 = find_history(g1, g2)
    h(p, t) = [exp_model(t, fit1); exp_model(t, fit2)]

    # problem
    prob = DDEProblem(DDEmodel, [g1_0[j], g2_0[j]], h, extrema(times), initial_guess;
                      constant_lags = [1/initial_guess[1], 1/initial_guess[2]])
    # algorithm to solve
    alg = MethodOfSteps(AutoTsit5(Rosenbrock23()); constrained=true)

    # objective function
    obj = build_loss_objective(prob, alg, L2Loss(times, data);
                               prob_generator=prob_generator,
                               verbose_opt = false)
    # optimizing
    results_dde = bboptimize(obj; SearchRange=bound,
                                    NumDimensions=6, TraceInterval=100, MasSteps=num_steps, Method=:adaptive_de_rand_1_bin_radiuslimited)
    # returning estimated parameteres and the objective function
    return exp.(best_candidate(results_dde))
end
