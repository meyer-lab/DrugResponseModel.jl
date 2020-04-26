"""
        In this file we want to estimate parameters of an ODE model describing the number of cells in G1 or G2 phase of the cell cycle 
"""

""" Make the transition matrix. """
function ODEjac(p::Vector{T}, nG1::Int, nG2::Int, nD1::Int, nD2::Int)::SparseMatrixCSC{T, Int64} where {T}
    # p = [alpha, beta, gamma1, gamma2, nG1, nG2, nD1, nD2]
    if nD1 == 0
        D1 = T[]
        diagD1 = T[]
    elseif nD1 == 1
        D1 = [0.0]
        diagD1 = -ones(nD1) * p[3]
    else
        D1 = [0.0; ones(nD1 - 1) * p[3]]
        diagD1 = -ones(nD1) * p[3]
    end

    if nD2 == 0
        D2 = T[]
        diagD2 = T[]
    elseif nD2 == 1
        D2 = [0.0]
        diagD2 = -ones(nD2) * p[4]
    else
        D2 = [0.0; ones(nD2 - 1) * p[4]]
        diagD2 = -ones(nD2) * p[4]
    end

    v1 = [-ones(nG1) * (p[3] + p[1]); -ones(nG2) * (p[4] + p[2]); diagD1; diagD2]
    v2 = [ones(nG1) * p[1]; ones(nG2 - 1) * p[2]; D1; D2]
    A = spdiagm(0 => v1, -1 => v2)

    A[1, nG1 + nG2] = 2 * p[2]
    if nD1 > 0
        A[nG1 + nG2 + 1, 1:nG1] = p[3] * ones(1, nG1)
    end
    if nD2 > 0
        A[nG1 + nG2 + nD1 + 1, (nG1 + 1):(nG1 + nG2)] = p[4] * ones(1, nG2)
    end

    return A
end


""" Predicts the model given a set of parametrs. """
function predict(p::Vector{T}, g_0, t) where {T}
    @assert length(p) == 9
    # Convert parameters to phase numbers
    nG1 = Int(floor(p[6]))
    nG2 = Int(floor(p[7]))
    nD1 = Int(floor(p[8]))
    nD2 = Int(floor(p[9]))

    if g_0 isa Real
        if nD1 == 0
            D1 = T[]
        else
            D1 = zeros(nD1)
        end

        if nD2 == 0
            D2 = T[]
        else
            D2 = zeros(nD2)
        end

        g_0 = [ones(nG1) * p[5] * g_0 / nG1; ones(nG2) * (1.0 - p[5]) * g_0 / nG2; D1; D2]
    end

    A = ODEjac(p, nG1, nG2, nD1, nD2)

    prob = ODEProblem((du, u, p, t) -> mul!(du, A, u), g_0, maximum(t))
    integrator = init(prob, VCABM(); save_on = false)

    G1 = Vector{T}(undef, length(t))
    G2 = Vector{T}(undef, length(t))

    ii = 1
    for (v, ttt) in TimeChoiceIterator(integrator, t)
        G1[ii] = sum(view(v, 1:nG1)) + sum(view(v, (nG1 + nG2 + 1):(nG1 + nG2 + nD1)))
        G2[ii] = sum(view(v, (nG1 + 1):(nG1 + nG2))) + sum(view(v, (nG1 + nG2 + nD1 + 1):(nG1 + nG2 + nD1 + nD2)))

        if ii == length(t)
            return abs.(G1), abs.(G2), v
        end

        ii += 1
    end
end


""" Calculates the cost function for a given set of parameters. """
function cost(p, g1, g2)
    t = LinRange(0.0, 0.5 * length(g1), length(g1))
    G1, G2, vecOut = predict(p, g1[1] + g2[1], t)

    return norm(G1 - g1) + norm(G2 - g2)
end


""" Given estimated parameters for each trial, solve the DDE model plot the predicted curve 
    for number of cells in G1, G2, or total, along with their corresponding real data,
    for a longer time which is 2 times of the original time (~195 hours) """
function ode_plotIt(params::Vector, g1::Matrix, g2::Matrix, pop, i::Int, title::String, legend::Any, ymax; tnew=0)
    t = LinRange(0.0, 0.5 * length(g1[:, 1]), length(g1[:, 1]))
    if tnew=0
        t_new = LinRange(0.0, 200, 400)
    else
        t_new = LinRange(0.0, tnew, 2*tnew)
    end
    G1, G2 = predict(params, g1[1] + g2[1], t_new)

    plot(
        t_new,
        G1,
        label = "G1 est",
        xlabel = "time [hours]",
        ylabel = "# of cells",
        xguidefontsize = 8,
        yguidefontsize = 8,
        lw = 2.0,
        alpha = 0.6,
        color = :green,
    )
    plot!(t, g1[:, i], label = "G1", markersize = 1.0, color = :darkgreen)
    plot!(t_new, G2, label = "G2 est", legend = legend, legendfontsize = 4, fg_legend = :transparent, lw = 2.0, alpha = 0.6, color = :sienna)
    plot!(t, g2[:, i], label = "G2", markersize = 1.0, color = :darkorange)
    plot!(t_new, G1 .+ G2, label = "total est", lw = 2.0, alpha = 0.6, color = :hotpink)
    plot!(t, pop[:, i], label = "total", markersize = 1.0, color = :indigo)
    plot!(annotation = [(100, ymax, text(title, 8))])
    ylims!((0.0, ymax))
end


""" Plot the data and curves for all concentrations. """
function ODEplot_all(params_ode, g1_l::Matrix, g2_l::Matrix, pop_l, conc::Array{Float64, 1})
    # plotting the fitted curves
    rl = [ode_plotIt(params_ode[:, i], g1_l, g2_l, pop_l, i, string(conc[i], " nM"), false, 80.0) for i = 1:4]
    r2 = [ode_plotIt(params_ode[:, i], g1_l, g2_l, pop_l, i, string(conc[i], " nM"), false, 40.0) for i = 5:7]
    r8 = ode_plotIt(params_ode[:, 8], g1_l, g2_l, pop_l, 8, string(conc[8], " nM"), :topleft, 40.0)
    plot(rl..., r2..., r8, layout = (2, 4))
end

function plotPercentage(params::Vector, g1::Matrix, g2::Matrix, pop, i::Int, title::String, legend::Any, ymax)
    t = LinRange(0.0, 0.5 * length(g1[:, 1]), length(g1[:, 1]))
    t_new = LinRange(0.0, 120, 200)
    G1, G2 = predict(params, g1[1] + g2[1], t_new)

    plot(
        t_new,
        100.0 * G1 ./ (G1 .+ G2),
        label = "G1 perc",
        xlabel = "time [hours]",
        ylabel = "% of cells",
        xguidefontsize = 8,
        yguidefontsize = 8,
        lw = 2.0,
        alpha = 0.6,
        color = :green,
    )
    plot!(t, 100.0 * g1[:, i] ./ pop[:, i], label = "G1", markersize = 1.0, color = :darkgreen)
    plot!(
        t_new,
        100.0 * G2 ./ (G1 .+ G2),
        label = "G2 perc",
        legend = legend,
        legendfontsize = 4,
        fg_legend = :transparent,
        lw = 2.0,
        alpha = 0.6,
        color = :sienna,
    )
    plot!(t, 100.0 * g2[:, i] ./ pop[:, i], label = "G2", markersize = 1.0, color = :darkorange)
    ylims!((0.0, ymax))
end

function ODEplot_allPerc(params_ode, g1_l::Matrix, g2_l::Matrix, pop_l, conc::Array{Float64, 1})
    # plotting the fitted curves
    rl = [plotPercentage(params_ode[:, i], g1_l, g2_l, pop_l, i, string(conc[i], " nM"), false, 110.0) for i = 1:7]
    r8 = plotPercentage(params_ode[:, 8], g1_l, g2_l, pop_l, 8, string(conc[8], " nM"), :topleft, 110.0)
    plot(rl..., r8, layout = (2, 4))
end
