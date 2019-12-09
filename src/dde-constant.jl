using DelayDiffEq, DiffEqParamEstim, BlackBoxOptim, LsqFit

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


# this function is used for DDE fitting with constant delay
function DDEsolve(g1, g2, g1_0, g2_0, initial_guess, j)

    # initial guess has 34 values, there are 8 * 4 parameters for transition and death rates for each conentration, and the last two are delays for G1 and G2.
    times = range(0.0; stop = 95.5, length = 192)

    # history function
    fit1, fit2 = find_history(g1, g2)
    h(p, t) = [exp_model(t, fit1); exp_model(t, fit2)]

    # the parameters for a specific concentration, j determines that
    pars = append!([initial_guess[4*j - 3], initial_guess[4*j - 2]], [initial_guess[33], initial_guess[34], initial_guess[4*j - 1], initial_guess[4*j]])
    # problem
    prob = DDEProblem(DDEmodel, [g1_0[j], g2_0[j]], h, extrema(times), pars;
                          constant_lags = [pars[5], pars[6]])

    # the algorithm to solve the dde problem
    alg = MethodOfSteps(AutoTsit5(Rosenbrock23()); constrained=true)

    solve(prob, alg)
end

function cost_function(initial_guess, g1, g2, g1_0, g2_0)
    times = range(0.0; stop = 95.5, length = 192)

    total_cost = zeros(8)
    for j in 1:8
        sol = DDEsolve(g1, g2, g1_0, g2_0, initial_guess, j)
        res1 = sum((sol(times, idxs=1).u - g1[:, j]).^2)
        res2 = sum((sol(times, idxs=2).u - g2[:, j]).^2)
        total_cost[j] = res1 + res2
    end
    # the output of the cost function must be an scalar, so we sum all the squared errors.
    return sum(total_cost)
end

function optimization_constantDelay(g1, g2, g1_0, g2_0, initial_guess, lower, upper, num_steps)

    # the cost function must just be a function of parameters
    cost_fcn(initial_guess) = cost_function(initial_guess, g1, g2, g1_0, g2_0)
    bound = collect(zip(lower, upper))
    # optimizing
    results_dde = bboptimize(cost_fcn; SearchRange=bound, MaxSteps=num_steps, Method=:adaptive_de_rand_1_bin_radiuslimited)

    new_guess = best_candidate(results_dde)
    return best_fitness(results_dde), new_guess
end