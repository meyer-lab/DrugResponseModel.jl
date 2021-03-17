""" Figure 1 including model cartoon, time series simulations and parameters."""
# remember: to load the simple ODE params do: JLD.load("G1_simpleODE.jld")["data"]

function plot_fig1(concs, g1, g1data, tite, G, subPlabel)
    time = LinRange(0.0, 95.0, 189)
    
    p = plot(time, g1, lw=2, legend=:topleft, label=["control" "$(concs[2]) nM" "$(concs[3]) nM" "$(concs[4]) nM" "$(concs[5]) nM" "$(concs[6]) nM" "$(concs[7]) nM"], fg_legend = :transparent, palette =:PuBu_8, title = tite, titlefont = Plots.font("Helvetica", 12), legendfont = Plots.font("Helvetica", 9), guidefont=Plots.font("Helvetica", 12), xtickfont=Plots.font("Helvetica", 12),ytickfont=Plots.font("Helvetica", 12), xlabel = "time [hr]", ylabel = "$G cell number", bottom_margin=1.5cm, top_margin=1.25cm, left_margin=1.25cm, right_margin=1.25cm)
    plot!(time, g1data, lw=2, linestyle = :dot, label=["" "" "" "" "" "" ""])
    annotate!(-25.0, 57.0, text(subPlabel, :black, :left, Plots.font("Helvetica Bold", 15)))
    ylims!((0.0, 50))
    p
end

function figureS1()

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
    p1 = plot_fig1(concs[:, 1], G1[:, :, 1], g1m[:, :, 1, 1], "Lapatinib", "G1", "b")
    p2 = plot_fig1(concs[:, 1], G2[:, :, 1], g2m[:, :, 1, 1], "Lapatinib", "G2", "c")
    p3 = plot_fig1(concs[:, 2], G1[:, :, 2], g1m[:, :, 2, 1], "Doxorubicin", "G1", "d")
    p4 = plot_fig1(concs[:, 2], G2[:, :, 2], g2m[:, :, 2, 1], "Doxorubicin", "G2", "e")

    fig1 = plot(p0, p1, p2, p3, p4, size=(2000, 400), layout=(1,5))
    savefig(fig1, "figureS1.svg")
end
