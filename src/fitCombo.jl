""" In this file the data will be organized for combinations to be fit to the model. """

########### 1. control, Palbociclib 50, Palbociclib 50 nM + lapatinibs [25 nM, 50 nM, 100 nM, 250 nM]
########### 2. control, Palbociclib 50, Palbociclib 50 nM + gemcitabines [5 nM, 10 nM, 17 nM, 30 nM]
########### 3. control, Gemcitabine 10, Gemcitabine 10 nM + palbociclibs [25 nM, 100 nM, 250 nM, 250 nM]
########### 4. control, Gemcitabine 10, Gemcitabine 10 nM + lapatinibs [25 nM, 50 nM, 100 nM, 250 nM]
########### 5. control, Doxorubicin 20, Doxorubicin 20 nM + gemcitabines [5 nM, 10 nM, 17 nM, 30 nM]
########### 6. control, Lapatinib  100, Lapatinib  100 nM + palbociclibs [25 nM, 100 nM, 250 nM, 250 nM]

GC = JLD.load("GC.jld")["GC"] # 2 x 193 x 6 x 6 which is G1/G2, timePoints, concentrations, trials
concens = zeros(6, 6)
concens[:, 1] = [0, 50, 75, 100, 150, 300]
concens[:, 2] = [0, 50, 55, 60, 67, 80]
concens[:, 3] = [0, 10, 35, 110, 260, 260]
concens[:, 4] = [0, 10, 35, 60, 110, 260]
concens[:, 5] = [0, 20, 25, 30, 37, 50]
concens[:, 6] = [0, 100, 125, 200, 350, 350]

function optim_all(concs::Array{Float64, 2}, g1::Array{Float64, 3}, g2::Array{Float64, 3}; maxiter = 500000)
    f(x) = DrugResponseModel.residHillAll(x, concs, g1, g2)

    lP = [minimum(concs); 0.01; 1e-9 * ones(12)]
    low = vcat(lP, lP, lP, lP, lP, lP, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9)
    hP = [maximum(concs); 10.0; 2 * ones(12)]
    high = vcat(hP, hP, hP, hP, hP, hP, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0)

    return DrugResponseModel.optimize_helper(f, low, high, maxiter)
end
