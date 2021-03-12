""" Figure 1: time-course simulation and data. """

function figure1()

    setGadflyTheme()
    concs, popul1, g1s1, g2s1 = load(189, 1);
    _, popul2, g1s2, g2s2 = load(189, 2);
    _, popul3, g1s3, g2s3 = load(189, 3);

    # find G1 std and mean ***** data ******
    g1S = cat(g1s1, g1s2, g1s3, dims=4);
    g2S = cat(g2s1, g2s2, g2s3, dims=4);
    g111 = (100 * g1S) ./ (g1S .+ g2S)
    g1m = mean(g111, dims = 4); # mean G1
    g1std = std(g111, dims = 4);
    g1m_min = g1m .- g1std;
    g1m_max = g1m .+ g1std;

    # find total stdand mean
    tot = cat(popul1, popul2, popul3, dims=4)
    totmean = mean(tot, dims = 4);
    tots_min = totmean .- std(tot, dims = 4);
    tots_max = totmean .+ std(tot, dims = 4);

    p = [44.184, 1.24076, 0.0692788, 0.0460918, 0.3822, 0.854034, 0.605391, 0.771326, 0.0138293, 0.00183699, 0.000293753, 0.0127534, 0.00011816, 0.0142541, 60.6069, 0.899573, 1.99993, 0.0748216, 1.99326, 0.468332, 1.99864, 1.22536, 0.000141615, 0.0318616, 0.000216899, 8.80158e-7, 0.598489, 0.00110572, 6.68492, 2.05974, 1.99936, 0.167588, 0.507586, 0.316074, 0.248084, 0.826596, 1.6164e-5, 3.10987e-6, 3.55996e-5, 7.73526e-6, 0.0774056, 8.26708e-5, 3.34656, 2.83739, 0.0907361, 0.108245, 1.9758, 1.96985, 1.9993, 0.210137, 0.0690636, 1.30442e-5, 0.0767181, 0.00991078, 6.87891e-5, 1.45086e-5, 18.2253, 1.1841, 1.00505, 0.0735852, 1.97326, 0.783828, 0.45769, 1.99355, 0.0519941, 0.000533671, 0.00204743, 9.52975e-5, 5.23806e-5, 0.0677505, 0.339953, 0.403341, 0.802518, 0.470576, 1.298, 0.423103];
    efcs = getODEparams(p, concs);

    # ******* model simulations ********
    G1 = zeros(189, 7, 5)
    G2 = zeros(189, 7, 5)

    t = LinRange(0.0, 95.0, 189)
    for k=1:5 # drug number
        for i = 1:7 # concentration number
            G1[:, i, k], G2[:, i, k], _ = predict(efcs[:, i, k], efcs[:, 1, k], t)
        end
    end

    p0 = DrugResponseModel.plot_percTimeCourse(G1[:, :, 1], G2[:, :, 1], g1m_min[:, :, 1, 1], g1m_max[:, :, 1, 1], concs[:, 1], "Lapatinib")
    p1 = DrugResponseModel.plot_totalTimeCourse(G1[:, :, 1] .+ G2[:, :, 1], tots_min[:, :, 1, 1], tots_max[:, :, 1, 1], concs[:, 1], "Lapatinib")
    pp1 = DrugResponseModel.plot_G1_params(efcs, 1, "G1 progression", 3)
    p2 = DrugResponseModel.plot_percTimeCourse(G1[:, :, 2], G2[:, :, 2], g1m_min[:, :, 2, 1], g1m_max[:, :, 2, 1], concs[:, 2], "Doxorubicin")
    p3 = DrugResponseModel.plot_totalTimeCourse(G1[:, :, 2] .+ G2[:, :, 2], tots_min[:, :, 2, 1], tots_max[:, :, 2, 1], concs[:, 2], "Doxorubicin")
    pp2 = DrugResponseModel.plot_G1_params(efcs, 7, "G1 cell death", 1)
    p4 = DrugResponseModel.plot_percTimeCourse(G1[:, :, 3], G2[:, :, 3], g1m_min[:, :, 3, 1], g1m_max[:, :, 3, 1], concs[:, 3], "Gemcitabine")
    p5 = DrugResponseModel.plot_totalTimeCourse(G1[:, :, 3] .+ G2[:, :, 3], tots_min[:, :, 3, 1], tots_max[:, :, 3, 1], concs[:, 3], "Gemcitabine")
    pp3 = DrugResponseModel.plot_G2_params(efcs, 3, "G2 progression", 3)
    p6 = DrugResponseModel.plot_percTimeCourse(G1[:, :, 4], G2[:, :, 4], g1m_min[:, :, 4, 1], g1m_max[:, :, 4, 1], concs[:, 4], "Paclitaxel")
    p7 = DrugResponseModel.plot_totalTimeCourse(G1[:, :, 4] .+ G2[:, :, 4], tots_min[:, :, 4, 1], tots_max[:, :, 4, 1], concs[:, 4], "Paclitaxel")
    pp4 = DrugResponseModel.plot_G2_params(efcs, 9, "G2 cell death", 1)
    p8 = DrugResponseModel.plot_percTimeCourse(G1[:, :, 5], G2[:, :, 5], g1m_min[:, :, 5, 1], g1m_max[:, :, 5, 1], concs[:, 5], "Palbociclib")
    p9 = DrugResponseModel.plot_totalTimeCourse(G1[:, :, 5] .+ G2[:, :, 5], tots_min[:, :, 5, 1], tots_max[:, :, 5, 1], concs[:, 5], "Palbociclib")
    draw(SVG("figure1.svg", 10inch, 14inch), plotGrid((5, 3), [nothing, p0, p1, pp1, p2, p3, pp2, p4, p5, pp3, p6, p7, pp4, p8, p9]))
end