""" This file handles functions for pairs of drugs. """

function residHillAll(hP, concentrations::Matrix, g1::Array, g2::Array, v::Int, u::Int, case_num::Int)
    res = 0.0

    # Solve for all drugs
    t = 1
    k = v
    for j = 1:2
        if case_num == 1
            hill = hP[[t, t + 2, t + 3, t + 1, t + 4, t + 5, t + 6, t + 7, t + 8, t + 9, t + 10, t + 11, t + 12]]
            res += DrugResponseModel.residHill(hill, concentrations[:, k], g1[:, :, k], g2[:, :, k])
            t += 13
        elseif case_num == 2
            hill = hP[[t, t + 2, t + 3, t + 1, t + 4, t + 5, t + 6, t + 7, t + 8, 19, 20, 21, 22]]
            res += DrugResponseModel.residHill(hill, concentrations[:, k], g1[:, :, k], g2[:, :, k])
            t += 9
        elseif case_num == 3
            hill = hP[[t, 23, t + 2, t + 1, 24, t + 3, t + 4, t + 5, t + 6, t + 7, t + 8, t + 9, t + 10]]
            res += DrugResponseModel.residHill(hill, concentrations[:, k], g1[:, :, k], g2[:, :, k])
            t += 11
        end
        k = u
    end

    return res
end

function optim_helper(concs, case_num::Int, i::Int, j::Int)
    low_piece = [1.0, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 0.25, 2, 2, 2, 2]
    high_piece = [5.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 0.75, 50, 50, 50, 50]

    if case_num == 1 # all parameters separated, nothing shared
        low = vcat(concs[1, i], low_piece, concs[1, j], low_piece)
        high = vcat(concs[7, i], high_piece, concs[7, j], high_piece)

    elseif case_num == 2 # only subphase numbers are shared
        low = vcat(concs[1, i], low_piece[1:8], concs[1, j], low_piece[1:8], 2, 2, 2, 2)
        high = vcat(concs[7, i], high_piece[1:8], concs[7, j], high_piece[1:8], 50, 50, 50, 50)

    elseif case_num == 3 # only control condition is shared
        low = vcat(concs[1, i], low_piece[1:2], low_piece[4], low_piece[6:end], concs[1, j], low_piece[1:2], low_piece[4:end], 1e-9, 1e-9)
        high = vcat(concs[7, i], high_piece[1:2], high_piece[4], high_piece[6:end], concs[7, j], high_piece[1:2], high_piece[4:end], 3.0, 3.0)
    end

    return low, high
end

""" Function to optimize pairs of drugs. """
function BBoptim_DrugPairs(concs::Array{Float64, 2}, g1::Array{Float64, 3}, g2::Array{Float64, 3}, case_num::Int, i::Int, j::Int; maxiter = 100000)

    hillCostAll(hillParams) = residHillAll(hillParams, concs, g1, g2, i, j, case_num)
    low, high = optim_helper(concs, case_num, i, j)

    results_ode = bboptimize(
        hillCostAll;
        SearchRange = collect(zip(low, high)),
        NumDimensions = length(low),
        TraceMode = :verbose,
        TraceInterval = 100,
        MaxSteps = maxiter,
    )

    return best_fitness(results_ode), best_candidate(results_ode)
end
