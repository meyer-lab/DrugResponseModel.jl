"""
This file fits Hill function to the parameters
"""

""" This functions takes in hill parameters for all the concentrations and calculates
DDE parameters, passes them to residual function and based off of these, optimizes the model
and estimates hill parameters. """
function residHill(hillParams::Vector, concentrations::Vector, g1::Matrix, g2::Matrix)
    res = 0.0
    params = getODEparams(hillParams, concentrations)

    # Solve each concentration separately
    for ii = 1:length(concentrations)
        resTemp = cost(params[1:9, ii], g1[:, ii], g2[:, ii])

        res += resTemp
    end

    return res
end


function grad_helper!(G, f, range, params)
    ForwardDiff.gradient!(G, f, params)
    costCenter = f(params)

    # Handle the integer-valued parameters
    for ii = range
        pp = copy(params)
        pp[ii] += 1.0
        costPlus = f(pp)
        pp[ii] -= 2.0
        costMin = f(pp)
        pt = floor(params[ii])
        poly = fit([pt-1, pt, pt+1], [costMin, costCenter, costPlus])
        polyD = derivative(poly)
        G[ii] = polyD(params[ii])
    end
end


function optimize_helper(f, g!, low, high)
    initial_x = low + (high - low) / 2.0

    ls = LineSearches.BackTracking()
    options = Optim.Options(outer_iterations = 2, show_trace = true, iterations = 3)
    results = optimize(f, g!, low, high, initial_x, Fminbox(LBFGS(linesearch = ls)), options)

    return Optim.minimum(results), Optim.minimizer(results)
end


""" Hill optimization function. """
function optimize_hill(conc_l::Vector, g1::Matrix, g2::Matrix; maxstep = 1E5)
    f(hillParams) = residHill(hillParams, conc_l, g1, g2)
    g!(G, params) = grad_helper!(G, f, 10:13, params)

    low = [minimum(conc_l), 1e-9, 1e-9, 0.1, 1e-9, 1e-9, 0.0, 0.0, 0.25, 3, 5, 0, 0]
    high = [maximum(conc_l), 1.0, 1.0, 10.0, 1.0, 1.0, 3.0, 3.0, 0.75, 50, 50, 50, 50]

    return optimize_helper(f, g!, low, high)
end

""" A function to convert the estimated hill parameters back to ODE parameters. """
function getODEparams(p::Vector, concentrations::Vector{Float64})
    effects = Matrix{eltype(p)}(undef, 9, length(concentrations))

    # Scaled drug effect
    xx = 1.0 ./ (1.0 .+ (p[1] ./ (concentrations .+ eps())) .^ p[4])

    # [EC50, left, right, steepness]
    effects[1, :] = p[2] .+ (p[3] - p[2]) .* xx
    effects[2, :] = p[5] .+ (p[6] - p[5]) .* xx
    effects[3, :] = p[7] .* xx
    effects[4, :] = p[8] .* xx
    effects[5, :] .= p[9]
    effects[6, :] .= floor(p[10])
    effects[7, :] .= floor(p[11])
    effects[8, :] .= floor(p[12])
    effects[9, :] .= floor(p[13])

    return effects
end
