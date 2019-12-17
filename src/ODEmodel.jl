"""
        In this file we want to estimate parameters of an ODE model describing the number of cells in G1 or G2 phase of the cell cycle .
"""


""" Time- and state-invariant model Jacobian. """
function update_coef(A, u, p, t)
    A[:, :] .= [-p[1]-p[3] 2*p[2]; p[1] -p[2]-p[4]]
end

""" Predicts the model given a set of parametrs. """
function predict(p, g1_0, g2_0, i, t)
    u0 = [g1_0[i], g2_0[i]]
    A = zeros(eltype(p), 2, 2)
    update_coef(A, nothing, p, nothing)
    Op = DiffEqArrayOperator(A, update_func=update_coef)
    prob = ODEProblem(Op, u0, extrema(t), p)
    solution = solve(prob, AutoTsit5(Rosenbrock23()))
    return prob, solution
end

""" Calculates the cost function for a given set of parameters. """
function cost(p, g1_0, g2_0, g1, g2, i)
    t = range(0.0; stop = 95.5, length = 192)
    _, solution = predict(p, g1_0, g2_0, i, t)
    res = zeros(2, 192)
    G1 = solution(t, idxs=1).u
    G2 = solution(t, idxs=2).u
    res[1, :] = (G1 - g1[:, i]).^2
    res[2, :] = (G2 - g2[:, i]).^2
    summ = sum(res[1,:]) + sum(res[2,:])
    return summ
end
""" Fit the ODE model to data. """
function ODEoptimizer(par::Array, i::Int, g1::Matrix, g2::Matrix, g1_0::Array, g2_0::Array)
    residuals(par) = cost(par, g1_0, g2_0, g1, g2, i)
    # lower and upper bounds for the parameters
    lower_bound = 0.0001*ones(4)
    upper_bound = 2*ones(4)
    bound = collect(zip(lower_bound, upper_bound))
    # global optimization with black box optimization
    results_ode = bboptimize(residuals; SearchRange=bound, NumDimensions=4, TraceMode=:silent, MaxSteps=10000)

    return best_fitness(results_ode), best_candidate(results_ode)
end


function ode_plotIt(params::Vector{Float64}, g1::Matrix, g2::Matrix, g1_0::Array, g2_0::Array, pop, i::Int, title::String, legend::Any)
    """ Given estimated parameters for each trial, solve the DDE model plot the predicted curve 
    for # of cells in G1, G2, or total, along with their corresponding real data,
    for a longer time which is 2 times of the original time (~195 hours)
    """
    t = LinRange(0.0, 95.5, 192)
    t_new = LinRange(0.0, 195.5, 292)
    _, solution = predict(params, g1_0, g2_0, i, t_new)

    plot(t_new, solution(t_new, idxs=1).u, label = "G1 est", dpi = 150, xlabel = "time [hours]", ylabel = "# of cells", lw=2.0, alpha = 0.6, color=:green)
    plot!(t, g1[:, i], label = "G1", dpi = 150, markersize = 1.0, color=:darkgreen)
    plot!(t_new, solution(t_new, idxs=2).u, label = "G2 est", legend=legend, legendfontsize=6, fg_legend = :transparent, lw=2.0, alpha = 0.6, color=:sienna)
    plot!(t, g2[:, i], label = "G2", dpi = 150, markersize = 1.0, color=:darkorange)
    plot!(t_new, (solution(t_new, idxs=2).u + solution(t_new, idxs=1).u), label = "total est", dpi = 150, lw=2.0, alpha = 0.6, color=:hotpink)
    plot!(t, pop[!, i], label = "total", dpi = 150, markersize = 1.0, color=:indigo)
    plot!( annotation=[ (75,90, text(title, 12)) ])
end

""" Plot the data and curves for all concentrations. """
function ODEplot_all(params_ode, g1_l::Matrix, g2_l::Matrix, g1_0_l::Array, g2_0_l::Array, pop_l)
    # plotting the fitted curves
    rl = [ode_plotIt(params_ode[:, i], g1_l, g2_l, g1_0_l, g2_0_l, pop_l, i, "", false) for i in 1:7]
    r8 = ode_plotIt(params_ode[:, 8], g1_l, g2_l, g1_0_l, g2_0_l, pop_l, 8, "", :topleft)
    plot(rl..., r8, layout = (2,4))
    plot!(size=(1200, 600), layout = (4,2), dpi=200)
    ylims!((0.0, 120.0))
end
