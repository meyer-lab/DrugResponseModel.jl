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


function optim_all(concs::Array{Float64, 2}, g1::Array{Float64, 3}, g2::Array{Float64, 3}; maxiter = 600000)
    f(x) = residHillAll(x, concs, g1, g2)

    lP = [minimum(concs); 0.01; 5e-3 * ones(12)]
    low = vcat(lP, lP, lP, lP, lP, 5e-3, 5e-3, 5e-3, 5e-3, 5e-3, 5e-3)
    hP = [maximum(concs); 10.0; 2.0 * ones(12)]
    high = vcat(hP, hP, hP, hP, hP, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5)

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
        51.0122,
        1.19478,
        0.0123853,
        0.197453,
        0.783039,
        6.53136e-5,
        1.35692e-6,
        0.284673,
        0.00521293,
        3.69958e-7,
        0.00913979,
        0.0258875,
        3.04229e-6,
        0.00527735,
        18.4107,
        1.38004,
        0.288625,
        9.6902e-9,
        0.787761,
        1.02151,
        1.99999,
        0.106618,
        4.35605e-9,
        0.0478454,
        1.22383e-7,
        1.04499e-7,
        0.381662,
        2.39835e-9,
        4.75582,
        1.78552,
        0.481014,
        0.404215,
        0.471125,
        0.187735,
        1.99999,
        0.255864,
        1.35294e-9,
        7.07919e-9,
        1.74332e-9,
        0.0672485,
        4.87662e-8,
        4.45473e-9,
        7.0734,
        2.47932,
        0.066145,
        5.62597e-8,
        1.94036,
        2.0,
        2.0,
        0.00866935,
        1.22435e-9,
        9.23547e-7,
        2.0,
        2.14921e-7,
        1.23361e-7,
        0.0174862,
        36.8515,
        1.11516,
        0.0806277,
        0.726529,
        1.92473,
        1.99999,
        1.97768,
        0.319934,
        2.65382e-9,
        6.12668e-9,
        0.0197645,
        1.06389e-6,
        5.28303e-8,
        0.0308013,
        0.196915,
        2.0,
        1.92313,
        2.0,
        1.99921,
        0.199044,
    ]
end
