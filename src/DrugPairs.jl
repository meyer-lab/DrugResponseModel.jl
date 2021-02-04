""" This file handles functions for pairs of drugs. """

function residHillAll(hP, concentrations::Matrix, g1::Array, g2::Array, v::Int, u::Int, case_num::Int)
    res = 0.0

    # Solve for all drugs
    t = 1
    k = v
    for j = 1:2
    hill = hP[[t, t + 1, 23, t + 2, 24, t + 3, 25, t + 4, 26, t + 5, t + 6, t + 7, t + 8, t + 9, t + 11]]
        res += DrugResponseModel.residHill(hill, concentrations[:, k], g1[:, :, k], g2[:, :, k])
        t += 11
        k = u
    end

    return res
end

function optim_helper(case_num::Int)
    # sharing control condition (min_a1, min_a2, min_b1,min_b2)
    low_piece = [0.5, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 0.25]
    high_piece = [1000.0, 10.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 0.75]

    low = vcat(low_piece, low_piece, 1e-9, 1e-9, 1e-9, 1e-9) # 26 params
    high = vcat(high_piece, high_piece, 3.0, 3.0, 3.0, 3.0)
    return low, high
end

""" Function to optimize pairs of drugs. """
function BBoptim_DrugPairs(concs::Array{Float64, 2}, g1::Array{Float64, 3}, g2::Array{Float64, 3}, case_num::Int, i::Int, j::Int; maxiter = 100000)
    hillCostAll(hillParams) = residHillAll(hillParams, concs, g1, g2, i, j, case_num)
    low, high = optim_helper(case_num)

    return optimize_helper(hillCostAll, low, high, maxiter)
end
