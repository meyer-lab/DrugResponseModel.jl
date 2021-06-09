""" In this file we fit all the drugs att once. """
function residHillAll(hP, concentrations::Matrix, g1::Array, g2::Array)
    res = 0.0

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


function optim_all(concs::Array{Float64, 2}, g1::Array{Float64, 3}, g2::Array{Float64, 3}; maxiter = 800000)
    f(x) = residHillAll(x, concs, g1, g2)

    lP = [minimum(concs); 0.01; 5e-3 * ones(12)]
    low = vcat(lP, lP, lP, lP, lP, 5e-3, 5e-3, 5e-3, 5e-3, 5e-3, 5e-3)
    hP = [maximum(concs); 10.0; 2.0 * ones(12)]
    high = vcat(hP, hP, hP, hP, hP, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0)

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
    return ps = [
            48.1186,
            1.17017,
            0.0344367,
            0.153542,
            0.343799,
            0.803251,
            1.99644,
            0.848162,
            0.009062,
            2.40783e-6,
            2.12734e-5,
            0.0538466,
            5.83618e-5,
            2.4157e-5,
            15.5429,
            1.43629,
            0.367771,
            0.0080823,
            0.454744,
            0.245463,
            0.218754,
            1.2806,
            5.99234e-6,
            0.0239985,
            1.27509e-6,
            4.87889e-6,
            0.0179952,
            0.146385,
            4.99849,
            1.90953,
            0.510128,
            0.220032,
            0.221798,
            1.60443,
            0.240652,
            0.950461,
            2.01533e-7,
            5.76098e-7,
            8.91899e-7,
            9.83827e-7,
            0.0734884,
            4.48427e-6,
            3.62437,
            2.74308,
            0.0939696,
            0.884336,
            1.02885,
            0.175659,
            0.312402,
            0.386282,
            1.71926e-6,
            0.373011,
            5.02502e-6,
            1.37135e-6,
            0.127517,
            2.9824e-6,
            37.9028,
            1.13736,
            0.0930045,
            0.591583,
            0.456692,
            1.10714,
            1.31904,
            1.3567,
            4.92751e-7,
            6.53767e-7,
            0.0295615,
            5.40729e-6,
            5.87169e-6,
            0.0568699,
            0.212288,
            1.28386,
            0.310092,
            1.44867,
            1.99996,
            0.471678,
        ]
end
