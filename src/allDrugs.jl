""" In this file we fit all the drugs att once. """

""" This function """
function getODEparamsAll(p::Array{Float64,1}, concentrations::Array{Float64,2})
    effects = zeros(9, length(concentrations[:, 1]), 5)

    k = 1
    # Scaled drug effect
    for i = 1:5
        xx = 1.0 ./ (1.0 .+ (p[k] ./ concentrations[:, i]) .^ p[k+1])
        effects[3, :, i] = p[k+4] .* xx # G1 death rate
        effects[4, :, i] = p[k+5] .* xx # G2 death rate

        if length(p) == 41
            effects[1, :, i] = p[36] .+ (p[k+2] - p[36]) .* xx # G1 prog. rate
            effects[2, :, i] = p[37] .+ (p[k+3] - p[37]) .* xx # G2 prog. rate
            effects[5, :, i] .= p[k+6] # percentage in G1
            k += 7
        elseif length(p) == 37
            effects[1, :, i] = p[31] .+ (p[k+2] - p[31]) .* xx
            effects[2, :, i] = p[32] .+ (p[k+3] - p[32]) .* xx
            k += 6
        end
    end

    if length(p) == 41
        effects[6:9, :, :] .= p[38:41, CartesianIndex(), CartesianIndex()]
    elseif length(p) == 37
        # percentage in G1, nG1, nG2, nD1, nD2
        effects[5:9, :, :] .= p[33:37, CartesianIndex(), CartesianIndex()]
    end

    return effects
end


function residHillAll(hillParams::Vector, concentrations::Matrix, g1::Array, g2::Array)
    res = 0.0

    # Solve for all drugs
    t = 1
    for j = 1:5
        hill = [
            hillParams[t],
            hillParams[36],
            hillParams[t+2],
            hillParams[t+1],
            hillParams[37],
            hillParams[t+3],
            hillParams[t+4],
            hillParams[t+5],
            hillParams[t+6],
            hillParams[38],
            hillParams[39],
            hillParams[40],
            hillParams[41],
        ]
        t += 7
        res += residHill(hill, concentrations[:, j], g1[:, :, j], g2[:, :, j])
    end

    return res
end


function optim_all(concs::Array{Float64,2}, g1::Array{Float64,3}, g2::Array{Float64,3})
    f(hillParams) = residHillAll(hillParams, concs, g1, g2)

    function g!(G, hillParams)
        ForwardDiff.gradient!(G, f, hillParams)
        costCenter = f(hillParams)

        # Handle the integer-valued parameters
        for ii = 37:40
            pp = copy(hillParams)
            pp[ii] += 1.0
            costPlus = f(pp)
            pp[ii] -= 2.0
            costMin = f(pp)
            pt = floor(hillParams[ii])
            poly = fit([pt-1, pt, pt+1], [costMin, costCenter, costPlus])
            polyD = derivative(poly)
            G[ii] = polyD(hillParams[ii])
        end

        println(G)
    end

    lP = [minimum(concs), 0.01, 0.05, 0.05, 0.00001, 0.00001, 0.3]
    low = vcat(lP, lP, lP, lP, lP, 1e-9, 1e-9, 3, 3, 2, 2)
    hP = [maximum(concs), 1.0, 1.0, 1.0, 0.1, 0.1, 0.7]
    high = vcat(hP, hP, hP, hP, hP, 1.0, 1.0, 10, 25, 50, 50)
    initial_x = low + (high - low) / 2.0

    ls = LineSearches.BackTracking()
    options = Optim.Options(outer_iterations = 2, show_trace = true, iterations = 20)
    results = optimize(f, g!, low, high, initial_x, Fminbox(LBFGS(linesearch = ls)), options)

    return Optim.minimum(results), Optim.minimizer(results)
end

""" To find IC50 or IC90 for each drug, separately."""
function find_IC(population, which)
    population = abs.(population[189, :, :] .- which * population[189, 1, :])
    IC = argmin(population, dims=2)
    return IC[1], IC[2], IC[3], IC[4], IC[5] # returns the argument
end
