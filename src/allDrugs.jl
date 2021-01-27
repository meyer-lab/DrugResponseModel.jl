""" In this file we fit all the drugs att once. """

""" This function returns the 41-length estimated parameters. """
function return_param41()
    p = [
        97.7434,
        1.41726,
        0.0501761,
        0.930205,
        0.111494,
        0.123604,
        0.552216,
        92.9045,
        1.40024,
        0.996396,
        0.0513087,
        0.521831,
        0.6535,
        0.566578,
        15.5317,
        2.3689,
        0.592433,
        0.999986,
        0.0283363,
        0.286975,
        0.503328,
        3.96929,
        4.62768,
        0.0512281,
        0.307528,
        0.549714,
        0.378717,
        0.50959,
        63.4248,
        0.976052,
        0.16582,
        0.740009,
        0.0572609,
        0.0776912,
        0.534201,
        0.734513,
        0.375555,
        16.8387,
        12.3945,
        30.176,
        14.5352,
    ]
    p
end




function getODEparamsPairs(p, conc)
    effects = zeros(eltype(p), 9, length(conc[:, 1]), 2)
    k = 1
    # Scaled drug effect
    for i = 1:2
        xx = 1.0 ./ (1.0 .+ (p[k] ./ conc[:, i]) .^ p[k + 1])
        effects[1, :, i] = p[k + 2] .+ (p[k + 3] - p[k + 2]) .* xx # G1 prog. rate
        effects[2, :, i] = p[k + 4] .+ (p[k + 5] - p[k + 4]) .* xx # G2 prog. rate
        effects[3, :, i] = p[k + 6] .* xx # G1 death rate
        effects[4, :, i] = p[k + 7] .* xx # G2 death rate
        effects[5, :, i] .= p[k + 8] # percentage in G1
        effects[6, :, i] .= p[k + 9]
        effects[7, :, i] .= p[k + 10]
        effects[8, :, i] .= p[k + 11]
        effects[9, :, i] .= p[k + 12]
        
        k += 13
    end
    return effects
end

function residHillAll(hP, concentrations::Matrix, g1::Array, g2::Array, v, u)
    res = 0.0

    # Solve for all drugs
    t = 1
    k = v
    for j = 1:2
        hill = hP[[t, t + 2, t + 3, t + 1, t + 4, t + 5, t + 6, t+ 7, t + 8, t+ 9, t + 10, t + 11, t + 12]]
        res += DrugResponseModel.residHill(hill, concentrations[:, k], g1[:, :, k], g2[:, :, k])
        t += 13
        k = u
    end

    return res
end



""" This function """
function getODEparamsAll(p, concentrations::Array{Float64, 2})
    effects = zeros(eltype(p), 9, length(concentrations[:, 1]), 5)

    k = 1
    # Scaled drug effect
    for i = 1:5
        xx = 1.0 ./ (1.0 .+ (p[k] ./ concentrations[:, i]) .^ p[k + 1])
        effects[3, :, i] = p[k + 4] .* xx # G1 death rate
        effects[4, :, i] = p[k + 5] .* xx # G2 death rate

        if length(p) == 41
            effects[1, :, i] = p[36] .+ (p[k + 2] - p[36]) .* xx # G1 prog. rate
            effects[2, :, i] = p[37] .+ (p[k + 3] - p[37]) .* xx # G2 prog. rate
            effects[5, :, i] .= p[k + 6] # percentage in G1
            k += 7
        elseif length(p) == 37
            effects[1, :, i] = p[31] .+ (p[k + 2] - p[31]) .* xx
            effects[2, :, i] = p[32] .+ (p[k + 3] - p[32]) .* xx
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


function residHillAll(hP, concentrations::Matrix, g1::Array, g2::Array)
    res = 0.0

    # Solve for all drugs
    t = 1
    for j = 1:5
        hill = hP[[t, 36, t + 2, t + 1, 37, t + 3, t + 4, t + 5, t + 6, 38, 39, 40, 41]]
        res += residHill(hill, concentrations[:, j], g1[:, :, j], g2[:, :, j])
        t += 7
    end

    return res
end

""" Organize Hill parameters for each drug in a 2D array. """
function Hill_p_eachDr(p)
    HillP = Matrix{eltype(p)}(undef, 7, 5)
    # each column: [EC50, steepness, max_g1_prog., max_g2_prog., max_g1_death, max_g2_death, %G1]
    j = 1
    for i = 1:5
        HillP[:, i] .= p[j:(j + 6)]
        j += 7
    end
    HillP
end

function optim_all(concs::Array{Float64, 2}, g1::Array{Float64, 3}, g2::Array{Float64, 3}; maxiter = 100000)
    f(x) = residHillAll(x, concs, g1, g2)
    g!(G, x) = grad_helper!(G, f, 38:41, x)

    lP = [minimum(concs), 0.01, 0.05, 0.05, 0.00001, 0.00001, 0.3]
    low = vcat(lP, lP, lP, lP, lP, 1e-9, 1e-9, 3, 3, 2, 2)
    hP = [maximum(concs), 1.0, 1.0, 1.0, 0.1, 0.1, 0.7]
    high = vcat(hP, hP, hP, hP, hP, 1.0, 1.0, 10, 25, 50, 50)

    return optimize_helper(f, g!, low, high, maxiter)
end

""" Takes in the 41 long Hill params and the index corresponding to the drug of interest, outputs the 9 long params at EC50. """
function EC50_params(p, i)
    d = DrugResponseModel.Hill_p_eachDr(p)
    # returns the following at EC50: [g1_prog., g2_prog, g1_death, g2_death, g1%, nG1, nG2, nD1, nD2]
    return append!([p[36] + (d[3, i] - p[36]) / 2, p[37] + (d[4, i] - p[37]) / 2, d[5, i] / 2, d[6, i] / 2, d[7, i]], p[38:41])
end
