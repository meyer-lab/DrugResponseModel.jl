using OrdinaryDiffEq, DiffEqParamEstim, Plots, CSV, Optim, DiffEqBase, BlackBoxOptim
gr()

"""
        In this file we want to estimate parameters of an ODE model describing the number of cells in G1 or G2 phase of the cell cycle 
"""


##---------------------------- Building the function and residuals -------------------------##
function ODEmodel(du, u, p, t)
    # p = [alpha, beta, gamma1, gamma2, initg1, initg2]
    du[1] = -p[1]*u[1] + 2*p[2]*u[2] - p[3]*u[1]
    du[2] = p[1]*u[1] - p[2]*u[2] - p[4]*u[2]
end

function ODEoptimizer(lower_bound::Array, upper_bound::Array, par::Array, i::Int, g1::Matrix, g2::Matrix, g1_0::Array, g2_0::Array)
    t = LinRange(0.0, 95.5, 192)
    times = range(0.0; stop = 95.5, length = 192)
    data = vcat(g1[:, i]', g2[:, i]')

    u0 = [g1_0[i], g2_0[i]]
    # generating the ODEproblem
    prob = ODEProblem(ODEmodel, u0, extrema(times), par)
    # solver algorithm
    alg = Tsit5()
    # lower and upper bounds for the parameters
    bound = collect(zip(lower_bound, upper_bound))
    # objective function
    obj = build_loss_objective(prob, alg, L2Loss(times, data);
                               verbose_opt = false)
    # global optimization with black box optimization
    results_ode = bboptimize(obj; SearchRange=bound, NumDimensions=4, TraceInterval=100)

    return best_candidate(results_ode)
end

function ode_plotIt(params::Array, g1::Matrix, g2::Matrix, g1_0::Array, g2_0::Array, pop, i::Int, title::String, legend::Any)
    """ Given estimated parameters for each trial, 
    solve the DDE model plot the predicted curve 
    for # of cells in G1, G2, or total, 
    along with their corresponding real data,
    for a longer time which is 2 times of the 
    original time (~195 hours)
    """
    t = LinRange(0.0, 95.5, 192)
    t_new = LinRange(0.0, 195.5, 292)
    tspan_new = (0.0, 195.5)
    u0_new = [g1_0[i], g2_0[i]]
    prob_new = ODEProblem(ODEmodel, u0_new, tspan_new, params)
    solution = solve(prob_new, Tsit5())

    plot(t_new, solution(t_new, idxs=1).u, label = "G1 est", dpi = 150, xlabel = "time [hours]", ylabel = "# of cells", lw=2.0, alpha = 0.6, color=:green)
    plot!(t, g1[:, i], label = "G1", dpi = 150, markersize = 1.0, color=:darkgreen)
    plot!(t_new, solution(t_new, idxs=2).u, label = "G2 est", legend=legend, legendfontsize=6, fg_legend = :transparent, lw=2.0, alpha = 0.6, color=:sienna)
    plot!(t, g2[:, i], label = "G2", dpi = 150, markersize = 1.0, color=:darkorange)
    plot!(t_new, (solution(t_new, idxs=2).u + solution(t_new, idxs=1).u), label = "total est", dpi = 150, lw=2.0, alpha = 0.6, color=:hotpink)
    plot!(t, pop[i], label = "total", dpi = 150, markersize = 1.0, color=:indigo)
    plot!( annotation=[ (75,90, text(title, 12)) ])
end