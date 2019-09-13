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

# this function is used for Hill fitting
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

# this function is used for DDE fitting
function DDEsolve(initial_guess, g1, g2, g1_0, g2_0, j)

    times = range(0.0; stop = 95.5, length = 192)
    # history function
    fit1, fit2 = find_history(g1, g2)
    h(p, t) = [exp_model(t, fit1); exp_model(t, fit2)]
    pars = append!(initial_guess[4*j - 3 : 4*j], [initial_guess[33], initial_guess[34]])
    # problem
    prob = DDEProblem(DDEmodel, [g1_0[j], g2_0[j]], h, extrema(times), pars;
                          constant_lags = [initial_guess[33], initial_guess[34]])
    alg = MethodOfSteps(AutoTsit5(Rosenbrock23()); constrained=true)

    println("hello")
    solve(prob)
    println("hello again")
end

function cost_function(initial_guess, g1, g2, g1_0, g2_0)
    times = range(0.0; stop = 95.5, length = 192)

    total_cost = zeros(8)
    for j in 1:8
        sol = DDEsolve(initial_guess, g1, g2, g1_0, g2_0, j)
        res1 = sum((sol(t, idxs=1).u - g1[:, j]).^2)
        res2 = sum((sol(t, idxs=2).u - g2[:, j]).^2)
        total_cost[j] = res1 + res2
    end
    return sum(total_cost)
end

function optimization(g1, g2, g1_0, g2_0, initial_guess, lower, upper, num_steps)

    # the cost function must just be a function of parameters
    cost_function(initial_guess) = cost_function(initial_guess, g1, g2, g1_0, g2_0)
    bound = collect(zip(lower, upper))

    # optimizing
    results_dde = bboptimize(cost_function; SearchRange=bound,
                                        NumDimensions=34,
                                        TraceInterval=100,
                                        MasSteps=num_steps,
                                        Method=:adaptive_de_rand_1_bin_radiuslimited)

    new_guess = best_candidate(results_dde)
    return best_fitness(results_dde), exp.(new_guess)
end
