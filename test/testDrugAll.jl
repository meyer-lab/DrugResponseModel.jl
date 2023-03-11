@testset "Fit All drug at once tests" begin
    tensor, names, concs, conds = DrugResponseModel.hcc_all()

    # test the local optimization function
    DrugResponseModel.optim_all(concs, tensor[1, :, :, :], tensor[2, :, :, :], maxiter = 2)
end
