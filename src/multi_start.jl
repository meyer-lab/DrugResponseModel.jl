function optimize_hill(conc::Vector, g1::Matrix, g2::Matrix, high, low; maxstep = 100000)

    f(x) = residHill(x, conc, g1, g2)
    g!(x, G) = grad_helper!(x, G, f, 10:13)

    # low = [minimum(conc), 1e-9, 1e-9, 0.1, 1e-9, 1e-9, 0.0, 0.0, 0.25, 3, 5, 0, 0]
    # high = [maximum(conc), 1.0, 1.0, 10.0, 1.0, 1.0, 3.0, 3.0, 0.75, 50, 50, 50, 50]

    return optimize_helper(f, g!, low, high, maxstep)
end


function multistart_minimization(multistart_method::TikTak, local_method,
                                 minimization_problem, conc, g1, g2)
    @unpack quasirandom_N, initial_N, θ_min, θ_max, θ_pow = multistart_method
    quasirandom_points = sobol_starting_points(minimization_problem, quasirandom_N)
    initial_points = _keep_lowest(quasirandom_points, initial_N)

    function _step(visited_minimum, (i, initial_point))
        θ = _weight_parameter(multistart_method, i)
        x = @. (1 - θ) * initial_point.location + θ * Optim.minimizer(visited_minimum)
        results = optimize_hill(conc, g1, g2, high, low)
        Optim.minimum(results) < Optim.minimum(visited_minimum) ? results : visited_minimum
    end
    
    foldl(_step, enumerate(initial_points); init = first(initial_points))
end
