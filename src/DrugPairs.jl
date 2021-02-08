""" This file handles functions for pairs of drugs. """

function residHillAll(hP, pControl::Vector, concentrations::Matrix, g1::Array, g2::Array, v::Int, u::Int)
    res = 0.0

    # Solve for all drugs
    t = 1
    k = v
    for j = 1:2
    hill = hP[[t, t + 1, 23, t + 2, 24, t + 3, 25, t + 4, 26, t + 5, t + 6, t + 7, t + 8, t + 9, t + 10]]
        res += DrugResponseModel.residHill(hill, pControl, concentrations[:, k], g1[:, :, k], g2[:, :, k])
        t += 11
        k = u
    end

    return res
end

function optim_helper()
    # sharing control condition (min_a1, min_a2, min_b1,min_b2)
    low_piece = [1.0, 0.01, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 0.35]
    high_piece = [500.0, 10.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.65]

    low = vcat(low_piece, low_piece, 1e-9, 1e-9, 1e-9, 1e-9) # 26 params
    high = vcat(high_piece, high_piece, 2.0, 2.0, 2.0, 2.0)
    return low, high
end

""" Function to optimize pairs of drugs. """
function BBoptim_DrugPairs(concs::Array{Float64, 2}, pControl::Vector, g1::Array{Float64, 3}, g2::Array{Float64, 3}, i::Int, j::Int; maxiter = 200000)
    hillCostAll(hillParams) = residHillAll(hillParams, pControl, concs, g1, g2, i, j)
    low, high = optim_helper()

    return optimize_helper(hillCostAll, low, high, maxiter)
end

function getODEparamspairs(p, conc)
    effects = zeros(eltype(p), 9, length(conc[:, 1]), 2)
    k = 1
    # Scaled drug effect
    for i = 1:2
        xx = 1.0 ./ (1.0 .+ (p[k] ./ conc[:, i]) .^ p[k + 1])

        # [EC50, left, right, steepness]
        effects[1, :, i] = p[23] .+ (p[k + 2] - p[23]) .* xx
        effects[2, :, i] = p[24] .+ (p[k + 3] - p[24]) .* xx
        effects[3, :, i] = p[25] .+ (p[k + 4] - p[25]) .* xx
        effects[4, :, i] = p[26] .+ (p[k + 5] - p[26]) .* xx
        effects[5, :, i] = p[k + 6] .* xx
        effects[6, :, i] = p[k + 7] .* xx
        effects[7, :, i] = p[k + 8] .* xx
        effects[8, :, i] = p[k + 9] .* xx
        effects[9, :, i] .= p[k + 10]
        
        k += 11
    end
    return effects
end