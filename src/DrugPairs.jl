""" This file handles functions for pairs of drugs. """

function residHillpairs(hP, concentrations::Matrix, g1::Array, g2::Array, v::Int, u::Int)
    res = 0.0

    # Solve for all drugs
    t = 1
    k = v # drug1 index
    for j = 1:2
    hill = hP[[t, t + 1, 29, t + 2, 30, t + 3, 31, t + 4, 32, t + 5, 33, t + 6, 34, t + 7, t + 8, t + 9, t + 10, t + 11, t + 12, t + 13]]
        res += residHill(hill, concentrations[:, k], g1[:, :, k], g2[:, :, k])
        t += 14
        k = u # drug 2 index
    end

    return res
end

function optim_helper()
    # sharing control condition (min_a1, min_a2, min_b1,min_b2)
    low_piece = [1.0, 0.01, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9]
    high_piece = [500.0, 10.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0]

    low = vcat(low_piece, low_piece, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9) # 26 params
    high = vcat(high_piece, high_piece, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0)
    return low, high
end

""" Function to optimize pairs of drugs. """
function BBoptim_DrugPairs(concs::Array{Float64, 2}, g1::Array{Float64, 3}, g2::Array{Float64, 3}, i::Int, j::Int; maxiter = 200000)
    hillCostAll(hillParams) = residHillpairs(hillParams, concs, g1, g2, i, j)
    low, high = optim_helper()

    return optimize_helper(hillCostAll, low, high, maxiter)
end


""" This function takes in the joint estimated parameters and returns both, individually, in one matrix. """
function estimated_params(concs, key)
    # 1. lpt and dox:
    p1 = [
        32.7512,
        1.16153,
        0.0253825,
        0.0440005,
        0.494171,
        0.680342,
        0.0148913,
        0.00367079,
        0.00946718,
        1.06803e-9,
        274.888,
        0.937185,
        2.0,
        0.13032,
        1.75981,
        2.0,
        1.78088e-9,
        1.19101e-9,
        6.15794e-9,
        0.445464,
        0.549706,
        2.0,
        0.212641,
        0.917995,
    ]
    # 2. lpt and gem:
    p2 = [
        36.1786,
        1.1477,
        0.14733,
        0.0370134,
        0.587787,
        0.813863,
        0.0201046,
        0.00570633,
        2.4371e-8,
        8.42837e-8,
        7.56176,
        1.98424,
        2.0,
        0.253205,
        0.197333,
        0.889442,
        4.7358e-9,
        3.37734e-9,
        0.0265278,
        3.09149e-8,
        0.432345,
        0.412341,
        0.341737,
        0.669241,
    ]
    # 3. lpt and tax:
    p3 = [
        33.3703,
        1.19243,
        0.0331291,
        0.0442161,
        0.513058,
        0.678679,
        0.0157238,
        0.00394221,
        0.00692375,
        1.16954e-9,
        2.97956,
        2.88671,
        1.71067,
        0.0638964,
        2.0,
        0.336735,
        9.77487e-9,
        0.0430892,
        1.09105e-8,
        0.00231518,
        0.498326,
        2.0,
        0.203987,
        0.725018,
    ]
    # 4. lpt and palbo:
    # 5. dox and gem:
    p5 = [
        19.5133,
        2.01514,
        0.0221494,
        0.456167,
        0.390407,
        2.0,
        0.0205179,
        0.0069432,
        1.01511e-9,
        0.14444,
        5.79555,
        1.60752,
        0.300264,
        0.492173,
        0.190679,
        0.990212,
        1.0469e-9,
        1.06128e-9,
        0.0279224,
        1.18084e-9,
        0.441948,
        2.0,
        0.206729,
        0.939593,
    ]
    # 6. dox and tax:
    p6 = [
        19.8052,
        1.61367,
        0.0272342,
        0.501662,
        0.474784,
        2.0,
        0.0242529,
        0.0488431,
        1.02986e-9,
        0.0805864,
        3.80863,
        2.51815,
        0.0918445,
        1.64435e-9,
        2.0,
        0.0720509,
        1.51462e-9,
        0.0288849,
        0.0319411,
        0.0186419,
        0.494724,
        2.0,
        0.201878,
        0.718836,
    ]
    # 7. dox and palbo:
    # 8. gem and tax:
    p8 = [
        9.7291,
        1.9724,
        0.515963,
        1.46484,
        1.16859,
        2.0,
        0.0371012,
        4.15166e-7,
        1.69266e-7,
        0.128604,
        2.97266,
        2.94866,
        1.8291,
        0.0649484,
        2.0,
        0.361753,
        7.58154e-8,
        0.0446008,
        1.75339e-7,
        1.55382e-7,
        0.504308,
        2.0,
        0.205093,
        0.838779,
    ]
    # 9. gem and palbo:
    # 10. tax and palbo:
    dictionary =
        Dict("1,2" => p1, "1,3" => p2, "1,4" => p3, "1,5" => p4, "2,3" => p5, "2,4" => p6, "2,5" => p7, "3,4" => p8, "3,5" => p9, "4,5" => p10)
    return DrugResponesModel.getODEparams(dictionary[key], concs)

end
