"""
        In this file we want to estimate parameters of an ODE model describing the number of cells in G1 or G2 phase of the cell cycle .
"""


""" Time- and state-invariant model Jacobian. """
function update_coef(A, u, p, t)
    A[:, :] .= [-p[1]-p[3] 2*p[2]; p[1] -p[2]-p[4]]
end

function ODEmodel(p)
    # p = [alpha, beta, gamma1, gamma2]
    @assert all(p .>= 0.0)
    return DiffEqArrayOperator([-p[1]-p[3] 2*p[2]; p[1] -p[2]-p[4]], update_func=update_coef)
end

""" Fit the ODE model to data. """
function ODEoptimizer(lower_bound::Array, upper_bound::Array, par::Array, i::Int, g1::Matrix, g2::Matrix, g1_0::Array, g2_0::Array)
    times = range(0.0; stop = 95.5, length = 192)
    data = vcat(g1[:, i]', g2[:, i]')

    u0 = [g1_0[i], g2_0[i]]
    # generating the ODEproblem
    prob = ODEProblem(ODEmodel(par), u0, extrema(times))
    # lower and upper bounds for the parameters
    bound = collect(zip(lower_bound, upper_bound))
    # objective function
    obj = build_loss_objective(prob, LinearExponential(), L2Loss(times, data); verbose_opt=false)
    # global optimization with black box optimization
    results_ode = bboptimize(obj; SearchRange=bound, NumDimensions=4, TraceMode=:silent, MaxSteps=50000)

    return best_candidate(results_ode)
end


function ode_plotIt(params::Vector{Float64}, g1::Matrix, g2::Matrix, g1_0::Array, g2_0::Array, pop, i::Int, title::String, legend::Any)
    """ Given estimated parameters for each trial, solve the DDE model plot the predicted curve 
    for # of cells in G1, G2, or total, along with their corresponding real data,
    for a longer time which is 2 times of the original time (~195 hours)
    """
    t = LinRange(0.0, 95.5, 192)
    t_new = LinRange(0.0, 195.5, 292)
    tspan_new = (0.0, 195.5)
    u0_new = [g1_0[i], g2_0[i]]
    prob_new = ODEProblem(ODEmodel(params), u0_new, tspan_new, params)
    solution = solve(prob_new, LinearExponential())

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
    r1 = ode_plotIt(params_ode[:, 1], g1_l, g2_l, g1_0_l, g2_0_l, pop_l, 1, "", false)
    r2 = ode_plotIt(params_ode[:, 2], g1_l, g2_l, g1_0_l, g2_0_l, pop_l, 2, "", false)
    r3 = ode_plotIt(params_ode[:, 3], g1_l, g2_l, g1_0_l, g2_0_l, pop_l, 3, "", false)
    r4 = ode_plotIt(params_ode[:, 4], g1_l, g2_l, g1_0_l, g2_0_l, pop_l, 4, "", false)
    r5 = ode_plotIt(params_ode[:, 5], g1_l, g2_l, g1_0_l, g2_0_l, pop_l, 5, "", false)
    r6 = ode_plotIt(params_ode[:, 6], g1_l, g2_l, g1_0_l, g2_0_l, pop_l, 6, "", false)
    r7 = ode_plotIt(params_ode[:, 7], g1_l, g2_l, g1_0_l, g2_0_l, pop_l, 7, "", false)
    r8 = ode_plotIt(params_ode[:, 8], g1_l, g2_l, g1_0_l, g2_0_l, pop_l, 8, "", :topleft)
    plot(r1, r2, r3, r4, r5, r6, r7, r8, layout = (2,4))
    plot!(size=(1200, 600), layout = (4,2), dpi=200)
    ylims!((0.0, 120.0))
end
