""" In this file we fit all the drugs att once. """
function residHillAll(hP, concentrations::Matrix, g1::Array, g2::Array)
    res = 0.0

    # Solve for all drugs
    t = 1
    for j = 1:5
        hill = hP[[t, t + 1, 56, t + 2, 57, t + 3, 58, t + 4, 59, t + 5, t + 6, t + 7, t + 8, t + 9, t + 10]]
        res += residHill(hill, concentrations[:, j], g1[:, :, j], g2[:, :, j])
        t += 11
    end

    return res
end

""" Organize Hill parameters for each drug in a 2D array. """
function Hill_p_eachDr(p)
    HillP = Matrix{eltype(p)}(undef, 11, 5)
    # each column: [EC50, steepness, max_g1,1_prog., max_g1,2_prog., max_g2,1_prog., max_g2,2_prog., max_g11_death, max_g12_death, max_g21_death, max_g22_death, %G1]
    j = 1
    for i = 1:5
        HillP[:, i] .= p[j:(j + 10)]
        j += 11
    end
    HillP
end

function optim_all(concs::Array{Float64, 2}, g1::Array{Float64, 3}, g2::Array{Float64, 3}; maxiter = 100000)
    f(x) = residHillAll(x, concs, g1, g2)

    lP = [minimum(concs), 0.01, 0.05, 0.05, 0.05, 0.05, 0.00001, 0.00001, 0.00001, 0.00001, 0.3]
    low = vcat(lP, lP, lP, lP, lP, 1e-9, 1e-9, 1e-9, 1e-9)
    hP = [maximum(concs), 10.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.7]
    high = vcat(hP, hP, hP, hP, hP, 1.0, 1.0, 1.0, 1.0)

    return optimize_helper(f, low, high, maxiter)
end

""" Takes in the Hill params and the index corresponding to the drug of interest, outputs the 9 long params at EC50. """
function EC50_params(p, i)
    d = DrugResponseModel.Hill_p_eachDr(p)
    # returns the following at EC50: [g1_prog., g2_prog, g1_death, g2_death, g1%]
    return append!([p[56] + (d[3, i] - p[56]) / 2, p[57] + (d[4, i] - p[57]) / 2, p[58] + (d[5, i] - p[58]) / 2, p[59] + (d[6, i] - p[59]) / 2, d[7, i] / 2, d[8, i] / 2, d[9, i] / 2, d[10, i] / 2, d[11, i]])
end
