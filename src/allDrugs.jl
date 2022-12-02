""" In this file we fit all the drugs att once. """
# import Pkg; Pkg.instantiate()
# Pkg.activate(".")
# using DrugResponseModel
# using Plots, LinearAlgebra, Statistics
# Plots.scalefontsizes(0.7)

# concs, _, g1s1, g2s1 = load(189, 1);
# _, _, g1s2, g2s2 = load(189, 2);
# _, _, g1s3, g2s3 = load(189, 3);

# G1S = cat(g1s1, g1s2, g1s3, dims=4);
# G2S = cat(g2s1, g2s2, g2s3, dims=4);

# g1m = mean(G1S, dims=4)[:, :, :, 1];
# g2m = mean(G2S, dims=4)[:, :, :, 1];

function residHillAll(hP, concentrations, g1::Array, g2::Array)
    res = 0.0

    # Solve for all drugs
    t = 1
    l = length(hP)
    for j = 1:size(g1)[3]
        hill = hP[[t:(t + 17); l-7:l]]
        res += residHill(hill, concentrations[j], g1[:, :, j], g2[:, :, j])
        t += 18
    end
    return res
end

""" Organize Hill parameters for each drug in a 2D array. """
function Hill_p_eachDr(p)
    HillP = Matrix{eltype(p)}(undef, 18, 6)
    # each column: [EC50, steepness, max_g1,1_prog., max_g1,2_prog., max_g2,1_prog., max_g2,2_prog., max_g11_death, max_g12_death, max_g21_death, max_g22_death]
    j = 1
    for i = 1:6
        HillP[:, i] .= p[j:(j + 17)]
        j += 18
    end
    HillP
end

function optim_all(concs, g1::Array{Float64, 3}, g2::Array{Float64, 3}; maxiter = 1000000)
    f(x) = residHillAll(x, concs, g1, g2)

    lP = [0.01; 1e-9 * ones(16)]
    low = vcat(minimum(concs[1]), lP, minimum(concs[2]), lP, minimum(concs[3]), lP, minimum(concs[4]), lP, minimum(concs[5]), lP, minimum(concs[6]), lP, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9)
    hP = [50.0; 4.0 * ones(16)]
    high = vcat(maximum(concs[1]), hP, maximum(concs[2]), hP, maximum(concs[3]), hP, maximum(concs[4]), hP, maximum(concs[5]), hP, maximum(concs[6]), hP, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0)

    return optimize_helper(f, low, high, maxiter)
end


function parameters()
    return ps = [
        34.6303,
        1.35959,
        0.360926,
        0.020735,
        2.14239e-6,
        0.00873742,
        0.291827,
        6.23808e-6,
        0.263093,
        1.86469,
        0.0638457,
        5.26201e-7,
        0.0357197,
        3.38634e-7,
        1.71977e-5,
        4.20002e-7,
        0.00385455,
        2.82233e-5,
        19.8599,
        1.48917,
        1.99969,
        0.179443,
        7.98855e-7,
        0.00299272,
        1.75263,
        1.99766,
        0.114934,
        1.87553,
        3.65302e-5,
        1.13364e-7,
        0.0543858,
        7.83016e-7,
        0.222225,
        0.0277672,
        6.12458e-7,
        0.144828,
        5.17319,
        1.8511,
        2.0,
        0.324095,
        0.503695,
        0.293661,
        0.416597,
        0.220342,
        0.261569,
        1.88773,
        3.2247e-6,
        1.82269e-7,
        1.2028e-6,
        1.09096e-7,
        8.27941e-8,
        0.0768722,
        2.0273e-7,
        2.94502e-6,
        2.55094,
        2.78168,
        0.986578,
        0.0695119,
        0.353046,
        0.0467065,
        0.220959,
        1.99977,
        0.318245,
        0.0387224,
        1.01961e-5,
        1.36537e-7,
        2.60407e-6,
        1.63094e-7,
        0.136062,
        6.90303e-6,
        1.4639e-7,
        0.0362207,
        35.9225,
        1.10826,
        1.99851,
        0.0419418,
        0.769692,
        0.26942,
        1.43254,
        1.54905,
        0.299407,
        1.86994,
        2.81656e-5,
        5.32358e-7,
        1.687e-6,
        1.45854e-6,
        2.97086e-5,
        0.0353744,
        0.0254184,
        2.32567e-5,
        1.99713,
        0.10733,
        1.99997,
        0.988148,
        1.74038,
        1.99998,
        0.20115,
        1.84866,
    ]
end

""" Takes in the Hill params and the index corresponding to the drug of interest, outputs the 9 long params at EC50. """
function EC50_params(p, i)
    d = Hill_p_eachDr(p)
    # returns the following at EC50: [g1_prog., g2_prog, g1_death, g2_death, g1%]
    return append!([
        p[91] + (d[3, i] - p[91]) / 2,
        p[92] + (d[4, i] - p[92]) / 2,
        p[93] + (d[5, i] - p[93]) / 2,
        p[94] + (d[6, i] - p[94]) / 2,
        p[95] + (d[7, i] - p[95]) / 2,
        p[96] + (d[8, i] - p[96]) / 2,
        p[97] + (d[9, i] - p[97]) / 2,
        p[98] + (d[10, i] - p[98]) / 2,
        d[11, i] / 2,
        d[12, i] / 2,
        d[13, i] / 2,
        d[14, i] / 2,
        d[15, i] / 2,
        d[16, i] / 2,
        d[17, i] / 2,
        d[18, i] / 2,
    ])
end

function parameters3()
    return ps = [
        47.5584,
        1.69019,
        1.92052,
        2.5885,
        3.87702e-6,
        0.0131544,
        0.357148,
        0.344382,
        2.52123,
        0.305383,
        1.33891e-6,
        1.03542e-6,
        0.0362224,
        5.19544e-8,
        2.11948e-7,
        2.26656e-8,
        5.32603e-7,
        0.0214524,
        4.68652,
        1.52275,
        1.92869,
        2.59219,
        2.50409,
        0.229096,
        0.632373,
        0.22838,
        7.87773e-6,
        0.301683,
        6.52621e-7,
        9.73155e-7,
        7.6348e-7,
        7.98396e-8,
        1.1389e-6,
        1.01559e-8,
        0.103743,
        1.18543e-7,
        36.9579,
        2.03863,
        1.92538,
        2.59073,
        2.50325,
        0.0703355,
        0.990774,
        0.238629,
        1.97474,
        0.305146,
        1.66371e-5,
        1.21966e-5,
        5.80293e-7,
        1.02364e-8,
        7.4474e-7,
        0.00203467,
        2.32727e-6,
        0.0140572,
        1.92073,
        2.58886,
        2.49959,
        0.121378,
        1.98288,
        0.565664,
        3.0,
        0.156844,
    ]
end
