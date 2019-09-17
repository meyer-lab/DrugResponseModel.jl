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
    remake(prob; p = exp_p, constant_lags = [exp_p[3], exp_p[4]])
end

# for hill fitting 
function ddesolve(g1, g2, g1_0, g2_0, params, j)
    times = range(0.0; stop = 95.5, length = 192)
    data = vcat(g1[:, j]', g2[:, j]')
    
    # history function
    fit1, fit2 = find_history(g1, g2)
    h(p, t) = [exp_model(t, fit1); exp_model(t, fit2)]

    # problem
    prob = DDEProblem(DDEmodel, [g1_0[j], g2_0[j]], h, extrema(times), params;
                      constant_lags = [params[3], params[4]])
    # algorithm to solve
    alg = MethodOfSteps(AutoTsit5(Rosenbrock23()); constrained=true)

    # objective function
    obj = build_loss_objective(prob, alg, L2Loss(times, data);
                               prob_generator=prob_generator,
                               verbose_opt=false)
    # returning estimated parameteres and the objective function
    return obj(params)
end


function optimization(g1, g2, g1_0, g2_0, initial_guess, j, lower, upper, num_steps)
    times = range(0.0; stop = 95.5, length = 192)
    data = vcat(g1[:, j]', g2[:, j]')
    
    # history function
    fit1, fit2 = find_history(g1, g2)
    h(p, t) = [exp_model(t, fit1); exp_model(t, fit2)]
    bound = collect(zip(lower, upper))
    # problem
    prob = DDEProblem(DDEmodel, [g1_0[j], g2_0[j]], h, extrema(times), initial_guess;
                      constant_lags = [initial_guess[3], initial_guess[4]])
    # algorithm to solve
    alg = MethodOfSteps(AutoTsit5(Rosenbrock23()); constrained=true)

    # objective function
    obj = build_loss_objective(prob, alg, L2Loss(times, data);
                               prob_generator=prob_generator,
                               verbose_opt = false)
    # optimizing
    println("blackbox optim begins ...")
    results_dde = bboptimize(obj; SearchRange=bound,
                                    NumDimensions=6, TraceInterval=100, MasSteps=num_steps, Method=:adaptive_de_rand_1_bin_radiuslimited)
    println("fitness before local optimization:  ")
    println(best_fitness(results_dde))
    new_guess = best_candidate(results_dde)
    return best_fitness(results_dde), exp.(new_guess)
#     println("local optimization begins")
#     optim_res = optimize(obj, lower, upper, new_guess, Fminbox(), Optim.Options(show_trace=true))
#     println("the fitness after local optimization : ")
#     println(Optim.minimum(optim_res))
#     return Optim.minimum(optim_res), exp.(Optim.minimizer(optim_res))
end
