"""
        In this file we want to estimate parameters of an ODE model describing the number of cells in G1 or G2 phase of the cell cycle 
"""


""" Actual differential equation. """
function ODEmodelFlex(du, u, p, t, nG1)
    # p = [alpha, beta, gamma1, gamma2, initg1, initg2]

    # G1
    du[2:nG1] .= p[1] .* u[1:(nG1-1)]
    du[1:nG1] .= -p[1] .* u[1:nG1] .- p[3] .* u[1:nG1]
    du[1] += 2*p[2]*u[end]

    # G2
    du[(nG1 + 2):end] .= p[2] .* u[(nG1 + 1):(length(u)-1)]
    du[(nG1 + 1):end] .= -p[1] .* u[(nG1 + 1):end] .- p[3] .* u[(nG1 + 1):end]
end


""" Predicts the model given a set of parametrs. """
function predict(p, g1_0, g2_0, i, t)
    u0 = [p[7]*g1_0[i], (1-p[7])*g1_0[i], p[8]*g2_0[i], (1-p[8])*g2_0[i]]
    prob = ODEProblem((a, b, c, d) -> ODEmodelFlex(a, b, c, d, 2), u0, extrema(t), p[1:6])
    solution = solve(prob, AutoTsit5(Rosenbrock23()))
    return prob, solution
end


""" Calculates the cost function for a given set of parameters. """
function cost(p, g1_0, g2_0, g1, g2, i)
    t = range(0.0; stop = 95.5, length = 192)
    _, solution = predict(p, g1_0, g2_0, i, t)
    res = zeros(2, 192)
    G1 = solution(t, idxs=1).u + solution(t, idxs=2).u
    G2 = solution(t, idxs=3).u + solution(t, idxs=4).u
    res[1, :] = (G1 - g1[:, i]).^2
    res[2, :] = (G2 - g2[:, i]).^2
    summ = sum(res[1,:]) + sum(res[2,:])
    return summ
end

""" Fit the ODE model to data. """
function ODEoptimizer4(p::Array, i::Int, g1::Matrix, g2::Matrix, g1_0::Array, g2_0::Array)
    
    residuals(p) = cost(p, g1_0, g2_0, g1, g2, i)
    # lower and upper bounds for the parameters
    lower_bound = zeros(8)
    upper_bound = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0]
    bound = collect(zip(lower_bound, upper_bound))

    # global optimization with black box optimization
    results_ode = bboptimize(residuals; SearchRange=bound, NumDimensions=8, TraceMode=:silent, MaxSteps=50000)

    return best_fitness(results_ode), best_candidate(results_ode)
end


function ode_plotIt4(params::Vector{Float64}, g1::Matrix, g2::Matrix, g1_0::Array, g2_0::Array, pop, i::Int, title::String, legend::Any)
    """ Given estimated parameters for each trial, solve the DDE model plot the predicted curve 
    for # of cells in G1, G2, or total, along with their corresponding real data,
    for a longer time which is 2 times of the original time (~195 hours)
    """
    t = range(0.0; stop = 95.5, length = 192)
    t_new = LinRange(0.0, 195.5, 292)
    _, solution = predict(params, g1_0, g2_0, i, t_new)

    plot(t_new, solution(t_new, idxs=1).u + solution(t_new, idxs=2).u, label = "G1 est", dpi = 150, xlabel = "time [hours]", ylabel = "# of cells", lw=2.0, alpha = 0.6, color=:green)
    plot!(t, g1[:, i], label = "G1", dpi = 150, markersize = 1.0, color=:darkgreen)
    plot!(t_new, solution(t_new, idxs=3).u + solution(t_new, idxs=4).u, label = "G2 est", legend=legend, legendfontsize=6, fg_legend = :transparent, lw=2.0, alpha = 0.6, color=:sienna)
    plot!(t, g2[:, i], label = "G2", dpi = 150, markersize = 1.0, color=:darkorange)
    plot!(t_new, (solution(t_new, idxs=1).u + solution(t_new, idxs=2).u+ solution(t_new, idxs=3).u + solution(t_new, idxs=4).u), label = "total est", dpi = 150, lw=2.0, alpha = 0.6, color=:hotpink)
    plot!(t, pop[!, i], label = "total", dpi = 150, markersize = 1.0, color=:indigo)
    plot!( annotation=[ (75,90, text(title, 12)) ])
end

""" Plot all the concentrations. """
function ODEplot_all4(params_ode, g1_l::Matrix, g2_l::Matrix, g1_0_l::Array, g2_0_l::Array, pop_l)
    # plotting the fitted curves
    r1 = ode_plotIt4(params_ode[:, 1], g1_l, g2_l, g1_0_l, g2_0_l, pop_l, 1, "", false)
    r2 = ode_plotIt4(params_ode[:, 2], g1_l, g2_l, g1_0_l, g2_0_l, pop_l, 2, "", false)
    r3 = ode_plotIt4(params_ode[:, 3], g1_l, g2_l, g1_0_l, g2_0_l, pop_l, 3, "", false)
    r4 = ode_plotIt4(params_ode[:, 4], g1_l, g2_l, g1_0_l, g2_0_l, pop_l, 4, "", false)
    r5 = ode_plotIt4(params_ode[:, 5], g1_l, g2_l, g1_0_l, g2_0_l, pop_l, 5, "", false)
    r6 = ode_plotIt4(params_ode[:, 6], g1_l, g2_l, g1_0_l, g2_0_l, pop_l, 6, "", false)
    r7 = ode_plotIt4(params_ode[:, 7], g1_l, g2_l, g1_0_l, g2_0_l, pop_l, 7, "", false)
    r8 = ode_plotIt4(params_ode[:, 8], g1_l, g2_l, g1_0_l, g2_0_l, pop_l, 8, "", :topleft)
    plot(r1, r2, r3, r4, r5, r6, r7, r8, layout = (2,4))
    plot!(size=(800, 400), layout = (4,2), dpi=200)
    ylims!((0.0, 120.0))
end

