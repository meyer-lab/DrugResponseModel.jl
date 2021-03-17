""" Figure 2: drug combination. """

function figure2()

    # data import
    gt1, gt2 = DrugResponseModel.import_combination("AU01001");
    gt1_2, gt2_2 = DrugResponseModel.import_combination("AU01101");
    gt1_3, gt2_3 = DrugResponseModel.import_combination("AU00901");
    concs, _, g1s1, g2s1 = load(189, 1);
    _, _, g1s2, g2s2 = load(189, 2);
    _, _, g1s3, g2s3 = load(189, 3);

    g1m = (g1s1 .+ g1s2 .+ g1s3) ./ 3; # pure data
    g2m = (g2s1 .+ g2s2 .+ g2s3) ./ 3; # pure data

    GS1 = cat(gt1, gt1_2, gt1_3, dims = 4);
    GS2 = cat(gt2, gt2_2, gt2_3, dims = 4);

    meanGS1 = mean(GS1, dims = 4);
    meanGS2 = mean(GS2, dims = 4);
    meanGS2[:, :, 19] .= mean(cat(gt2[:, :, 19], gt2_2[:, :, 19], dims=3), dims=3)[:, :, 1]

    # params from fitting all 5 drugs at once
    p = [44.184, 1.24076, 0.0692788, 0.0460918, 0.3822, 0.854034, 0.605391, 0.771326, 0.0138293, 0.00183699, 0.000293753, 0.0127534, 0.00011816, 0.0142541, 60.6069, 0.899573, 1.99993, 0.0748216, 1.99326, 0.468332, 1.99864, 1.22536, 0.000141615, 0.0318616, 0.000216899, 8.80158e-7, 0.598489, 0.00110572, 6.68492, 2.05974, 1.99936, 0.167588, 0.507586, 0.316074, 0.248084, 0.826596, 1.6164e-5, 3.10987e-6, 3.55996e-5, 7.73526e-6, 0.0774056, 8.26708e-5, 3.34656, 2.83739, 0.0907361, 0.108245, 1.9758, 1.96985, 1.9993, 0.210137, 0.0690636, 1.30442e-5, 0.0767181, 0.00991078, 6.87891e-5, 1.45086e-5, 18.2253, 1.1841, 1.00505, 0.0735852, 1.97326, 0.783828, 0.45769, 1.99355, 0.0519941, 0.000533671, 0.00204743, 9.52975e-5, 5.23806e-5, 0.0677505, 0.339953, 0.403341, 0.802518, 0.470576, 1.298, 0.423103];
    efcs = getODEparams(p, concs);
    t = LinRange(0.0, 95.0, 189)

    # Bliss on Model
    LPT_PLB = DrugResponseModel.AllBliss_params(efcs[:, :, 1], efcs[:, :, 5])
    LPT_GEM = DrugResponseModel.AllBliss_params(efcs[:, :, 1], efcs[:, :, 3])
    g1_Bliss_model1, g2_Bliss_model1, _ = predict(LPT_PLB[:, 6, 5], LPT_PLB[:, 1, 1], t) # palbo 50 & lapatinib 100
    g1_Bliss_model2, g2_Bliss_model2, _ = predict(LPT_PLB[:, 6, 6], LPT_PLB[:, 1, 1], t) # palbo 100 & lapatinib 100
    g1_Bliss_model3, g2_Bliss_model3, _ = predict(LPT_GEM[:, 6, 6], LPT_GEM[:, 1, 1], t) # gem 10 % lapatinib 100
    g1_Bliss_model4, g2_Bliss_model4, _ = predict(LPT_GEM[:, 6, 7], LPT_GEM[:, 1, 1], t) # gem 30 % lapatinib 100

    # Experimental data
    # meanG2[2, :, 5] # palbo 50 & lpt 100
    # meanG2[2, :, 14] # palbo 100 & lpt 100
    # meanG2[2, :, 11] # gem 10 & lpt 100
    # meanG2[2, :, 19] # gem 30 & lpt 100

    # Bliss on data
    Bliss_cellnum = zeros(189, 8, 8, 10)
    for i=1:189
        Bliss_cellnum[i, :, :, :] .= blissCellNum(g1m, g2m, i)
    end

    p1 = plot(t, g1_Bliss_model1, label="G1 Model", xlabel="time [hr]", ylabel="cell number", title="palbo 50 nM & lapt 100 nM")
    plot!(t, g2_Bliss_model1, label="G2 Model")
    plot!(t, meanGS2[1, 1:189, 5], label="G1 data")
    plot!(t, meanGS2[2, 1:189, 5], label="G2 data")
    annotate!(-25.0, 57.0, text("a", :black, :left, Plots.font("Helvetica Bold", 15)))
    p2 = plot(t, g1_Bliss_model2, label="G1 Model", xlabel="time [hr]", ylabel="cell number", title="palbo 100 nM & lapt 100 nM")
    plot!(t, g2_Bliss_model2, label="G2 Model")
    plot!(t, meanGS2[1, 1:189, 14], label="G1 data")
    plot!(t, meanGS2[2, 1:189, 14], label="G2 data")
    annotate!(-25.0, 57.0, text("b", :black, :left, Plots.font("Helvetica Bold", 15)))
    p3 = plot(t, g1_Bliss_model3, label="G1 Model", xlabel="time [hr]", ylabel="cell number", title="gem 10 nM & lapt 100 nM")
    plot!(t, g2_Bliss_model3, label="G2 Model")
    plot!(t, meanGS2[1, 1:189, 11], label="G1 data")
    plot!(t, meanGS2[2, 1:189, 11], label="G2 data")
    annotate!(-25.0, 57.0, text("c", :black, :left, Plots.font("Helvetica Bold", 15)))
    p4 = plot(t, g1_Bliss_model4, label=["G1 Model" "G2 Model"], xlabel="time [hr]", ylabel="cell number", title="gem 30 nM & lapt 100 nM")
    plot!(t, g2_Bliss_model4, label="G2 Model")
    plot!(t, meanGS2[1, 1:189, 19], label="G1 data")
    plot!(t, meanGS2[2, 1:189, 19], label="G2 data")
    annotate!(-25.0, 57.0, text("d", :black, :left, Plots.font("Helvetica Bold", 15)))
    fig = plot(p1, p2, p3, p4, size=(900, 200), layout=(1, 4))
    savefig(fig, "figure2.svg")
end