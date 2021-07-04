""" In this file we fit all the drugs att once. """
function residHillAll(x, concentrations::Matrix, g1::Array, g2::Array)
    res = 0.0
    hP = convertParams(x)

    # Solve for all drugs
    t = 1
    for j = 1:5
        hill = hP[[t:(t + 13); 71:76]]
        res += residHill(hill, concentrations[:, j], g1[:, :, j], g2[:, :, j])
        t += 14
    end

    return res
end

""" Organize Hill parameters for each drug in a 2D array. """
function Hill_p_eachDr(p)
    HillP = Matrix{eltype(p)}(undef, 14, 5)
    # each column: [EC50, steepness, max_g1,1_prog., max_g1,2_prog., max_g2,1_prog., max_g2,2_prog., max_g11_death, max_g12_death, max_g21_death, max_g22_death]
    j = 1
    for i = 1:5
        HillP[:, i] .= p[j:(j + 13)]
        j += 14
    end
    HillP
end

function convertParams(p)
    ps = zeros(eltype(p), length(p))
    ps[[1, 2, 15, 16, 29, 30, 43, 44, 57, 58]] .= p[[1, 2, 15, 16, 29, 30, 43, 44, 57, 58]]
    ps[3:8] .= p[3:8] .* p[9:14]
    ps[9:14] .= p[3:8] .* (1 .- p[9:14])
    ps[17:22] .= p[17:22] .* p[23:28]
    ps[23:28] .= p[17:22] .* (1 .- p[23:28])
    ps[31:36] .= p[31:36] .* p[37:42]
    ps[37:42] .= p[31:36] .* (1 .- p[37:42])
    ps[45:50] .= p[45:50] .* p[51:56]
    ps[51:56] .= p[45:50] .* (1 .- p[51:56])
    ps[59:64] .= p[59:64] .* p[65:70]
    ps[65:70] .= p[59:64] .* (1 .- p[65:70])
    ps[71:76] .= p[71:76]
    ps
end

function optim_all(concs::Array{Float64, 2}, g1::Array{Float64, 3}, g2::Array{Float64, 3}; maxiter = 800000)
    f(x) = residHillAll(x, concs, g1, g2)

    lP = [minimum(concs); 0.01; 5e-6 * ones(12)]
    low = vcat(lP, lP, lP, lP, lP, 5e-6, 5e-6, 5e-6, 5e-6, 5e-6, 5e-6)
    hP = [maximum(concs); 100.0; 10.0 * ones(6); ones(6)]
    high = vcat(hP, hP, hP, hP, hP, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0)

    return optimize_helper(f, low, high, maxiter)
end

""" Takes in the Hill params and the index corresponding to the drug of interest, outputs the 9 long params at EC50. """
function EC50_params(p, i)
    d = Hill_p_eachDr(p)
    # returns the following at EC50: [g1_prog., g2_prog, g1_death, g2_death, g1%]
    return append!([
        p[71] + (d[3, i] - p[71]) / 2,
        p[72] + (d[4, i] - p[72]) / 2,
        p[73] + (d[5, i] - p[73]) / 2,
        p[74] + (d[6, i] - p[74]) / 2,
        p[75] + (d[7, i] - p[75]) / 2,
        p[76] + (d[8, i] - p[76]) / 2,
        d[9, i] / 2,
        d[10, i] / 2,
        d[11, i] / 2,
        d[12, i] / 2,
        d[13, i] / 2,
        d[14, i] / 2,
    ])
end

function parameters()
    # ps = [64.7017, 1.23308, 0.0497766, 0.0354078, 0.298186, 0.352937, 5.99025, 5.51624, 0.996764, 0.111, 0.912474, 0.882585, 0.999991, 0.999993, 10.3459, 60.3182, 0.0235755, 0.189566, 0.208661, 0.168131, 5.67416, 5.50529, 0.999703, 0.856264, 0.999872, 0.999973, 0.485208, 0.999991, 4.28723, 1.96154, 0.369573, 0.206118, 0.164357, 9.99424, 0.762876, 5.5222, 0.999996, 0.999991, 0.999991, 0.99997, 0.802116, 0.99997, 1.8953, 3.22331, 0.134315, 0.54119, 0.180243, 0.62034, 5.99009, 5.51988, 0.999972, 0.69186, 0.999996, 0.86843, 0.999971, 0.999985, 38.1342, 1.10607, 0.0920182, 0.464887, 0.270508, 9.98713, 5.98948, 5.52591, 0.999994, 0.999995, 0.899942, 0.999983, 0.999992, 0.999995, 0.188227, 4.06045, 0.160831, 9.99678, 5.99233, 5.52016]
    ps = [41.0344, 1.19049, 0.0218886, 0.0227231, 0.0196711, 0.365798, 0.633352, 0.352679, 0.356918, 0.962542, 0.0300771, 0.999995, 0.999998, 0.971031, 129.837, 0.872128, 1.99988, 0.162396, 1.92442, 1.99744, 1.05315, 0.343968, 1.0, 1.0, 0.996255, 0.990439, 0.633404, 0.894326, 4.3579, 1.84243, 0.0258195, 0.361174, 0.279478, 0.409747, 0.0779829, 0.499509, 5.02441e-6, 1.0, 1.0, 1.0, 0.910216, 1.0, 2.01498, 3.39502, 0.174019, 0.140168, 1.99998, 1.9872, 0.74008, 0.313305, 0.584252, 0.999994, 0.959432, 0.999999, 0.999999, 1.0, 28.9282, 1.17123, 0.283458, 0.107746, 1.96149, 1.99992, 0.709239, 0.49967, 1.0, 0.999998, 0.979032, 0.988758, 0.999917, 0.933098, 2.0, 0.201884, 1.91476, 1.97491, 0.665105, 0.249237]
    return convertParams(ps)
end
