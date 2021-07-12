""" In this file we fit all the drugs att once. """
function residHillAll(hP, concentrations::Matrix, g1::Array, g2::Array)
    res = 0.0

    # Solve for all drugs
    t = 1
    for j = 1:5
        hill = hP[[t:(t + 17); 91:98]]
        for i=3:10
            res += 20 * (maximum([0, (hill[i] - hill[i + 16])]))^2
        end
        res += residHill(hill, concentrations[:, j], g1[:, :, j], g2[:, :, j])
        t += 18
    end
    return res
end

""" Organize Hill parameters for each drug in a 2D array. """
function Hill_p_eachDr(p)
    HillP = Matrix{eltype(p)}(undef, 14, 5)
    # each column: [EC50, steepness, max_g1,1_prog., max_g1,2_prog., max_g2,1_prog., max_g2,2_prog., max_g11_death, max_g12_death, max_g21_death, max_g22_death]
    j = 1
    for i = 1:5
        HillP[:, i] .= p[j:(j + 17)]
        j += 18
    end
    HillP
end


function optim_all(concs::Array{Float64, 2}, g1::Array{Float64, 3}, g2::Array{Float64, 3}; maxiter = 800000)
    f(x) = residHillAll(x, concs, g1, g2)

    lP = [minimum(concs); 0.01; 5e-9 * ones(16)]
    low = vcat(lP, lP, lP, lP, lP, 5e-9, 5e-9, 5e-9, 5e-9, 5e-9, 5e-9, 5e-9, 5e-9)
    hP = [maximum(concs); 10.0; 2.0 * ones(16)]
    high = vcat(hP, hP, hP, hP, hP, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0)

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
    # return ps = [
    #     51.0122,
    #     1.19478,
    #     0.0123853,
    #     0.197453,
    #     0.783039,
    #     6.53136e-5,
    #     1.35692e-6,
    #     0.284673,
    #     0.00521293,
    #     3.69958e-7,
    #     0.00913979,
    #     0.0258875,
    #     3.04229e-6,
    #     0.00527735,
    #     18.4107,
    #     1.38004,
    #     0.288625,
    #     9.6902e-9,
    #     0.787761,
    #     1.02151,
    #     1.99999,
    #     0.106618,
    #     4.35605e-9,
    #     0.0478454,
    #     1.22383e-7,
    #     1.04499e-7,
    #     0.381662,
    #     2.39835e-9,
    #     4.75582,
    #     1.78552,
    #     0.481014,
    #     0.404215,
    #     0.471125,
    #     0.187735,
    #     1.99999,
    #     0.255864,
    #     1.35294e-9,
    #     7.07919e-9,
    #     1.74332e-9,
    #     0.0672485,
    #     4.87662e-8,
    #     4.45473e-9,
    #     7.0734,
    #     2.47932,
    #     0.066145,
    #     5.62597e-8,
    #     1.94036,
    #     2.0,
    #     2.0,
    #     0.00866935,
    #     1.22435e-9,
    #     9.23547e-7,
    #     2.0,
    #     2.14921e-7,
    #     1.23361e-7,
    #     0.0174862,
    #     36.8515,
    #     1.11516,
    #     0.0806277,
    #     0.726529,
    #     1.92473,
    #     1.99999,
    #     1.97768,
    #     0.319934,
    #     2.65382e-9,
    #     6.12668e-9,
    #     0.0197645,
    #     1.06389e-6,
    #     5.28303e-8,
    #     0.0308013,
    #     0.196915,
    #     2.0,
    #     1.92313,
    #     2.0,
    #     1.99921,
    #     0.199044,
    # ]
    # return ps = [54.6595, 1.22536, 0.00817339, 0.137652, 0.135604, 0.00503426, 0.418209, 1.92837, 0.00586817, 0.00500449, 0.00500894, 0.005006, 0.00500197, 0.0415639, 18.4401, 1.42434, 0.249482, 0.00500245, 0.00500466, 1.68958, 0.127406, 1.94002, 0.00500112, 0.0490445, 0.0050006, 0.170381, 0.00500512, 0.00519855, 5.57593, 1.88777, 0.47617, 
    # 0.765006, 0.456894, 0.177639, 0.392646, 1.96138, 0.00500035, 0.00500038, 0.00500026, 0.0313934, 0.00500016, 0.00500459, 5.92172, 2.55462, 0.0477424, 0.0564463, 0.975153, 1.694, 0.0315034, 1.92996, 0.00500041, 0.0050786, 1.99989, 0.00504788, 0.0158125, 0.00500111, 26.5997, 1.13797, 0.0722648, 0.840981, 0.313468, 1.51079, 0.396128, 1.96376, 0.0050001, 0.00500265, 0.00500003, 0.00504917, 0.0144283, 0.00500107, 0.153677, 2.49999, 1.77488, 1.67971, 0.284039, 1.91376]
    # return ps = [45.4589, 1.92586, 0.0571684, 0.0168042, 0.0273424, 0.35549, 0.111112, 0.962342, 0.385335, 1.67885, 6.42138e-5, 0.0122718, 6.19759e-6, 3.22143e-5, 5.02488e-6, 1.46972e-5, 9.12192e-6, 0.105813, 14.79, 2.04263, 0.162507, 1.69772, 0.244912, 0.00604, 2.48685, 1.1905, 0.129931, 2.28785, 4.74559e-5, 0.000171909, 1.66239e-5, 0.0440693, 0.137718, 4.42172e-5, 8.79518e-6, 0.165207, 8.39994, 0.995218, 1.00763, 0.0393673, 2.49813, 2.499, 2.24026, 0.32845, 2.06678, 2.49882, 8.65157e-5, 1.19405e-5, 1.73906e-5, 0.000227363, 1.75223e-5, 0.0443169, 0.0886547, 0.103747, 3.13534, 3.5912, 0.108102, 0.780275, 0.122165, 0.0582937, 1.4578, 1.25933, 0.336067, 0.0844491, 1.65217e-5, 4.88924e-5, 8.54358e-6, 0.00326196, 0.413141, 9.99396e-6, 6.53298e-6, 0.0370135, 373.385, 9.99871, 8.01008e-5, 0.0494254, 0.0018814, 0.0219445, 2.49694, 1.25899, 0.123104, 1.7199, 0.246123, 2.49978, 2.49998, 1.25524, 0.247188, 2.27526, 6.91976e-6, 0.000132209, 0.563699, 2.06385, 1.47051, 0.798052, 0.598185, 2.28319, 0.907361, 0.764529]
    # return ps = [17.4572, 1.26195, 0.232028, 0.0447656, 0.0388862, 0.0911821, 0.252258, 2.91288, 2.9967, 1.96374, 2.46961e-5, 0.0236608, 1.24154e-5, 1.85966e-5, 0.00348492, 0.00516664, 0.000372484, 2.19017e-5, 24.0314, 1.36392, 0.222, 1.16523, 1.04524e-5, 0.00121498, 0.128687, 2.05857, 2.99962, 1.94971, 1.59343e-5, 9.1866e-5, 0.020808, 0.010247, 2.67206e-6, 0.000253706, 0.569355, 9.46428e-5, 131.435, 14.5948, 2.98234, 0.072218, 0.00217116, 2.83815, 0.00290212, 2.91134, 2.04009, 2.98192, 0.0015774, 0.000898055, 0.00137922, 0.0108708, 2.12793, 2.96006, 2.22106, 0.000349361, 2.78691, 2.70097, 0.164359, 0.207115, 0.0357208, 0.903378, 0.424601, 0.181006, 0.095295, 1.94148, 1.03816e-6, 8.05488e-6, 1.3503e-5, 0.666043, 2.08268e-5, 0.0216351, 0.0874651, 7.17821e-5, 21.1571, 1.17583, 0.184586, 1.17055, 0.0528577, 0.722757, 0.373982, 2.16658, 2.99395, 1.96009, 3.0769e-5, 0.00040568, 1.77357e-5, 0.00117944, 0.0166885, 0.0711195, 0.000126459, 1.93502e-5, 0.162448, 1.168, 0.27075, 2.9947, 0.231515, 2.91077, 2.98858, 1.92348]
    # return [41.4148, 1.19704, 2.74466, 7.55344e-5, 0.389681, 0.0117411, 0.00326188, 0.203705, 0.272936, 2.77613, 0.00273757, 0.0257183, 0.00206648, 9.10361e-6, 0.0332894, 0.014988, 0.00276915, 0.00200817, 30.2894, 1.32737, 2.77292, 2.18424, 2.09899, 0.025542, 7.8431e-6, 0.361595, 0.133812, 2.77817, 0.000217725, 0.000118328, 2.5783e-5, 2.35624e-7, 0.151818, 6.87167e-6, 9.48412e-6, 0.767885, 19.8025, 1.84845, 2.79601, 2.19004, 2.11248, 0.0398322, 1.5959, 0.450643, 0.103282, 2.76713, 2.2319, 5.68957e-5, 0.000211497, 2.77531e-6, 1.03374e-5, 3.15075e-5, 1.44032e-5, 0.14002, 6.85996, 2.47142, 0.777784, 0.701617, 0.657885, 0.0326617, 2.99115, 0.889788, 0.0106635, 2.77142, 2.63485e-5, 2.52547e-5, 2.93699e-5, 8.25194e-6, 2.99958, 5.83857e-5, 0.0177439, 2.55075e-5, 32.7892, 1.06486, 2.76285, 2.176, 2.0887, 0.0295006, 1.22962, 0.547603, 0.321995, 2.77235, 2.95002e-5, 0.000440076, 0.00112154, 1.69491e-5, 0.0360555, 0.00330994, 0.0235594, 0.000234048, 2.77165, 2.17894, 2.1003, 0.100506, 2.98554, 0.877331, 0.212145, 2.77002]
    return ps = [34.6303, 1.35959, 0.360926, 0.020735, 2.14239e-6, 0.00873742, 0.291827, 6.23808e-6, 0.263093, 1.86469, 0.0638457, 5.26201e-7, 0.0357197, 3.38634e-7, 1.71977e-5, 4.20002e-7, 0.00385455, 2.82233e-5, 19.8599, 1.48917, 1.99969, 0.179443, 7.98855e-7, 0.00299272, 1.75263, 1.99766, 0.114934, 1.87553, 3.65302e-5, 1.13364e-7, 0.0543858, 7.83016e-7, 0.222225, 0.0277672, 6.12458e-7, 0.144828, 5.17319, 1.8511, 2.0, 0.324095, 0.503695, 0.293661, 0.416597, 0.220342, 0.261569, 1.88773, 3.2247e-6, 1.82269e-7, 1.2028e-6, 1.09096e-7, 8.27941e-8, 0.0768722, 2.0273e-7, 2.94502e-6, 2.55094, 2.78168, 0.986578, 0.0695119, 0.353046, 0.0467065, 0.220959, 1.99977, 0.318245, 0.0387224, 1.01961e-5, 1.36537e-7, 2.60407e-6, 1.63094e-7, 0.136062, 6.90303e-6, 1.4639e-7, 0.0362207, 35.9225, 1.10826, 1.99851, 0.0419418, 0.769692, 0.26942, 1.43254, 1.54905, 0.299407, 1.86994, 2.81656e-5, 5.32358e-7, 1.687e-6, 1.45854e-6, 2.97086e-5, 0.0353744, 0.0254184, 2.32567e-5, 1.99713, 0.10733, 1.99997, 0.988148, 1.74038, 1.99998, 0.20115, 1.84866]
end
