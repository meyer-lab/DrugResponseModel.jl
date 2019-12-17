"""
        In this file we want to estimate parameters of an ODE model describing the number of cells in G1 or G2 phase of the cell cycle 
"""


""" Actual differential equation. """
function ODEmodelFlex(du, u, p, t, nG1)
    # p = [alpha, beta, gamma1, gamma2]

    # G1
    du[1] = 2*p[2]*u[end]
    du[2:(nG1+1)] .= p[1] .* u[1:nG1]
    du[1:nG1] .+= -p[1] .* u[1:nG1] .- p[3] .* u[1:nG1]

    # G2
    du[(nG1 + 2):end] .= p[2] .* u[(nG1 + 1):(length(u)-1)]
    du[(nG1 + 1):end] .+= -p[1] .* u[(nG1 + 1):end] .- p[4] .* u[(nG1 + 1):end]
end


""" Predicts the model given a set of parametrs. """
function predict(p, g1_0, g2_0, t, nG1, nG2)
    u0 = [ones(nG1)*g1_0/nG1  ones(nG2)*g2_0/nG2]
    prob = ODEProblem((a, b, c, d) -> ODEmodelFlex(a, b, c, d, nG1), u0, extrema(t), p)
    solution = solve(prob, Tsit5())
    return prob, solution
end

""" Calculates the cost function for a given set of parameters. """
function cost(p, g1_0, g2_0, g1, g2, nG1, nG2)
    t = range(0.0; stop = 95.5, length = 192)
    _, solution = predict(p, g1_0, g2_0, t, nG1, nG2)

    G1 = sum(solution(t, idxs=1:nG1), dims=1)
    G2 = sum(solution(t, idxs=nG1+1:nG1+nG2), dims=1)

    return sum((vec(G1) - g1).^2 + (vec(G2) - g2).^2)
end

""" Fit the ODE model to data. """
function ODEoptimizer(p::Array, i::Int, g1::Matrix, g2::Matrix, g1_0::Array, g2_0::Array, nG1, nG2)
    
    residuals(p) = cost(p, g1_0[i], g2_0[i], g1[:, i], g2[:, i], nG1, nG2)
    # lower and upper bounds for the parameters
    bound = collect(zip(zeros(4), ones(4)))

    # global optimization with black box optimization
    results_ode = bboptimize(residuals; SearchRange=bound, NumDimensions=4, TraceMode=:silent, MaxSteps=50000)

    return best_fitness(results_ode), best_candidate(results_ode)
end

""" Remake the problem by creating dual type. """
function remakeProblem(prob, p, SaveAt)
  _prob = remake(prob;u0=convert.(eltype(p),prob.u0),p=p)
  solve(_prob,AutoTsit5(Rosenbrock23()), saveat=SaveAt;abstol=1e-6, reltol= 1e-6)
end

""" Turing for 2-eq ODE. """
function turingODE(params_ode, g1_0, g2_0, i)
    t = range(0.0; stop = 95.5, length = 192)
    tp = collect(t)

    prob, sol = predict(params_ode, g1_0, g2_0, t, 2, 2)
    newsol = zeros(2, 192)
    newsol[1, :] = sol(t, idxs=1).u
    newsol[2, :] = sol(t, idxs=2).u

    @model bayesODE(prob, x, tp, params_ode) = begin
      alpha ~ truncated(Normal(0.5, 0.2), 0.0, 1.0)
      beta ~ truncated(Normal(0.5, 0.2), 0.0, 1.0)
      gamma1 ~ truncated(Normal(0.5, 0.2), 0.0, 1.0)
      gamma2 ~ truncated(Normal(0.5, 0.2), 0.0, 1.0)

      # gather parameters and solve equation
      p = [alpha, beta, gamma1 ,gamma2]
      sol_tmp = remakeProblem(prob, params_ode, tp)
          N = length(tp)

#           fill_length = length(tp) - length(sol_tmp.u)

#           for i in 1:fill_length
#             if eltype(sol_tmp.u) <: Number
#               push!(sol_tmp.u, Inf)
#             else
#               push!(sol_tmp.u, fill(Inf, size(sol_tmp[1])))
#             end

#           end

      for i in 1:N
        x[:,i] ~ MvNormal(sol_tmp.u[i], [0.01,0.01])
      end
    end
    chain = sample(bayesODE(prob, newsol, tp, params_ode), NUTS(0.65), 2000)
    return chain
end

function ode_plotIt4(params::Vector{Float64}, g1::Matrix, g2::Matrix, g1_0::Array, g2_0::Array, pop, i::Int, title::String, legend::Any)
    """ Given estimated parameters for each trial, solve the DDE model plot the predicted curve 
    for # of cells in G1, G2, or total, along with their corresponding real data,
    for a longer time which is 2 times of the original time (~195 hours)
    """
    t = range(0.0; stop = 95.5, length = 192)
    t_new = LinRange(0.0, 195.5, 292)
    _, solution = predict(params, g1_0[i], g2_0[i], t_new, 2, 2)

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
    rl = [ode_plotIt4(params_ode[:, ii], g1_l, g2_l, g1_0_l, g2_0_l, pop_l, ii, "", false) for ii in 1:7]
    r8 = ode_plotIt4(params_ode[:, 8], g1_l, g2_l, g1_0_l, g2_0_l, pop_l, 8, "", :topleft)
    plot(rl..., r8, layout = (2,4))
    plot!(size=(800, 400), layout = (4,2), dpi=200)
    ylims!((0.0, 120.0))
end

