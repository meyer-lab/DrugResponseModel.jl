"""
        In this file we want to estimate parameters of an ODE model describing the number of cells in G1 or G2 phase of the cell cycle 
"""

""" Make the transition matrix. """
function ODEjac(p, dt::Real, nG1::Int, nG2::Int)
    # p = [alpha, beta, gamma1, gamma2, nG1, nG2]
    A = diagm( 0 => vec([-ones(nG1)*(p[3]+p[1]); -ones(nG2)*(p[4]+p[2])]),
              -1 => vec([ones(nG1)*p[1];  ones(nG2-1)*p[2]]))
    A[1, end] = 2*p[2]

    rmul!(A, dt)
    A = LinearAlgebra.exp!(A)

    return A
end


""" Predicts the model given a set of parametrs. """
function predict(p, g1_0::Real, g2_0::Real, t, nG1::Int, nG2::Int)
    # Some assumptions
    @assert t[1] == 0.0

    v = [ones(nG1)*p[5]*(g1_0+g2_0)/nG1; ones(nG2)*(1.0-p[5])*(g1_0+g2_0)/nG2]
    A = ODEjac(p, t[2], nG1, nG2)

    G1 = Vector{eltype(p)}(undef, length(t))
    G2 = Vector{eltype(p)}(undef, length(t))

    for ii in 1:length(G1)
        G1[ii] = sum(v[1:nG1])
        G2[ii] = sum(v[nG1+1:nG1+nG2])

        v = A*v
    end

    return G1, G2
end


""" Calculates the cost function for a given set of parameters. """
function cost(p, g1_0, g2_0, g1, g2, nG1::Int, nG2::Int)
    v = [ones(nG1)*p[5]*(g1_0+g2_0)/nG1; ones(nG2)*(1.0-p[5])*(g1_0+g2_0)/nG2]
    temp = similar(v)
    A = ODEjac(p, 0.5, nG1, nG2)

    cost = 0.0
    for ii in 1:length(g1)
        @inbounds cost += (sum(v[1:nG1]) - g1[ii])^2
        @inbounds cost += (sum(v[nG1+1:nG1+nG2]) - g2[ii])^2

        @inbounds LinearAlgebra.mul!(temp, A, v)
        @inbounds copyto!(v, temp)
    end

    return cost
end


""" Fit the ODE model to data. """
function ODEoptimizer(i::Int, g1::Matrix, g2::Matrix, g1_0::Array, g2_0::Array)
    residuals(p) = cost(p, g1_0[i], g2_0[i], g1[:, i], g2[:, i], Int(floor(p[6])), Int(floor(p[7])))
    # lower and upper bounds for the parameters
    lower = [0.0, 0.0, 0.0, 0.0, 0.0, 1, 1]
    upper = [3.0, 3.0, 3.0, 3.0, 1.0, 70, 70]
    bound = collect(zip(lower, upper))

    # global optimization with black box optimization
    results_ode = bboptimize(residuals; SearchRange=bound, NumDimensions=7, TraceMode=:silent, MaxSteps=50000)

    return best_fitness(results_ode), best_candidate(results_ode)
end


""" Given estimated parameters for each trial, solve the DDE model plot the predicted curve 
    for number of cells in G1, G2, or total, along with their corresponding real data,
    for a longer time which is 2 times of the original time (~195 hours) """
function ode_plotIt(params::Vector, g1::Matrix, g2::Matrix, g1_0::Array, g2_0::Array, pop, i::Int, title::String, legend::Any)
    t = LinRange(0.0, 95.5, 192)
    t_new = LinRange(0.0, 150, 200)
    G1, G2 = predict(params, g1_0[i], g2_0[i], t_new, Int(floor(params[6])), Int(floor(params[7])))

    plot(t_new, G1, label = "G1 est", dpi = 150, xlabel = "time [hours]", ylabel = "# of cells", lw=2.0, alpha = 0.6, color=:green)
    plot!(t, g1[:, i], label = "G1", dpi = 150, markersize = 1.0, color=:darkgreen)
    plot!(t_new, G2, label = "G2 est", legend=legend, legendfontsize=6, fg_legend = :transparent, lw=2.0, alpha = 0.6, color=:sienna)
    plot!(t, g2[:, i], label = "G2", dpi = 150, markersize = 1.0, color=:darkorange)
    plot!(t_new, G1 .+ G2, label = "total est", dpi = 150, lw=2.0, alpha = 0.6, color=:hotpink)
    plot!(t, pop[!, i], label = "total", dpi = 150, markersize = 1.0, color=:indigo)
    plot!( annotation=[ (75,90, text(title, 12)) ])
end


""" Plot the data and curves for all concentrations. """
function ODEplot_all(params_ode, g1_l::Matrix, g2_l::Matrix, g1_0_l::Array, g2_0_l::Array, pop_l)
    # plotting the fitted curves
    rl = [ode_plotIt(params_ode[:, i], g1_l, g2_l, g1_0_l, g2_0_l, pop_l, i, "", false) for i in 1:7]
    r8 = ode_plotIt(params_ode[:, 8], g1_l, g2_l, g1_0_l, g2_0_l, pop_l, 8, "", :topleft)
    plot(rl..., r8, layout = (2,4))
    plot!(size=(900, 400), layout = (4,2), dpi=200)
    ylims!((0.0, 120.0))
end
