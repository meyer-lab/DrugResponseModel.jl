""" Figure S2 including model cartoon, time series simulations and parameters."""
# remember: to load the simple ODE params do: JLD.load("G1_simpleODE.jld")["data"]

function figureS2()

    concs, popul1, g1s1, g2s1 = load(189, 1);
    _, popul2, g1s2, g2s2 = load(189, 2);
    _, popul3, g1s3, g2s3 = load(189, 3);

    # find G1 std and mean ***** data ******
    g1S = cat(g1s1, g1s2, g1s3, dims=4);
    g2S = cat(g2s1, g2s2, g2s3, dims=4);
    g1m = mean(g1S, dims = 4); # mean G1
    g2m = mean(g2S, dims = 4); # mean G2

    G1 = JLD.load("G1_simpleODE.jld")["data"]
    G2 = JLD.load("G2_simpleODE.jld")["data"]

    p0 = plot(legend=false, grid=false, foreground_color_subplot=:white, top_margin=1.5cm)
    p1 = DrugResponseModel.plot_fig1(concs[:, 1], G1[:, :, 1], g1m[:, :, 1, 1], "Lapatinib", "G1", "b")
    p2 = DrugResponseModel.plot_fig1(concs[:, 1], G2[:, :, 1], g2m[:, :, 1, 1], "Lapatinib", "G2", "c")
    p3 = DrugResponseModel.plot_fig1(concs[:, 2], G1[:, :, 2], g1m[:, :, 2, 1], "Doxorubicin", "G1", "d")
    p4 = DrugResponseModel.plot_fig1(concs[:, 2], G2[:, :, 2], g2m[:, :, 2, 1], "Doxorubicin", "G2", "e")

    fig1 = plot(p0, p1, p2, p3, p4, size=(2000, 400), layout=(1,5))
    savefig(fig1, "figureS2.svg")
end
