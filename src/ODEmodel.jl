using OrdinaryDiffEq, DiffEqParamEstim, Plots, CSV, Optim, DiffEqBase, LeastSquaresOptim


"""
        In this file we want to estimate parameters of an ODE model describing the number of cells in G1 or G2 phase of the cell cycle 
"""


##---------------------------- Building the function and residuals -------------------------##
function ODEmodel(du, u, p, t)
    # p = [alpha, beta, gamma1, gamma2, initg1, initg2]
    du[1] = -p[1]*u[1] + 2*p[2]*u[2] - p[3]*u[1]
    du[2] = p[1]*u[1] - p[2]*u[2] - p[4]*u[2]
end

function ODEsolve(par::Array, i::Int, g1_0::Array, g2_0::Array)
    t = LinRange(0.0, 95.5, 192)
    tspan = (0.0, 95.5)
    u0 = [g1_0[i], g2_0[i]]
    prob = ODEProblem(ODEmodel, u0, tspan, par)
    solve(prob, Tsit5())
end

function residual(par, i, g1_0, g2_0, g1, g2)
    t = LinRange(0.0, 95.5, 192)
    res = zeros(2, 192)
    sol = ODEsolve(par, i, g1_0, g2_0)
    res[1, :] = sol(t, idxs=1).u - g1[:, i]
    res[2, :] = sol(t, idxs=2).u - g2[:, i]
    return res
end

function ode_optimIt(initial_guess::Array, lower_bound::Array, upper_bound::Array, i::Int, g1::Matrix, g2::Matrix, g1_0, g2_0)
    residuals(pp) = residual(pp, i, g1_0, g2_0, g1, g2)
    results_dde = optimize(residuals, initial_guess, Dogleg(), lower = lower_bound, upper = upper_bound)
    return results_dde.minimizer
end

function ode_plotIt(params::Array, g1::Matrix, g2::Matrix, g1_0::Array, g2_0::Array, pop, i::Int, title::String)
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

    plot(t_new, solution(t_new, idxs=1).u, label = "G1 est", dpi = 150, title = title, xlabel = "time [hours]", ylabel = "# of cells")
    plot!(t, g1[:, i], label = "G1", dpi = 150)
    plot!(t_new, solution(t_new, idxs=2).u, label = "G2 est", legend=:topleft, dpi = 150)
    plot!(t, g2[:, i], label = "G2", dpi = 150)
    plot!(t_new, (solution(t_new, idxs=2).u + solution(t_new, idxs=1).u), label = "total est", dpi = 150)
    plot!(t, pop[i], label = "total", dpi = 150)
end