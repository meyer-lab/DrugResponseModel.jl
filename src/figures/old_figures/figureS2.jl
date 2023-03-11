""" Figure S1 to show the data and model comparison in other 3 drugs. """



""" To plot the baseline Bliss overlayed with experimental data. """
function plot_figss()
    concs, popul1, g1s1, g2s1 = load(189, 1)
    expD = JLD.load("g.jld")["g"]
    C = output_Bliss_cellnum()

    D = (expD[1, 1:189, [1, 3, 4, 5, 6], [1, 2, 3, 4, 6]] .+ expD[2, 1:189, [1, 3, 4, 5, 6], [1, 2, 3, 4, 6]])
    p1 = DrugResponseModel.plot_fig1(concs[:, 1], C[:, :, 3, 4], D[:, :, 1], "palbo50 + lpts", "total", "A")
    p2 = DrugResponseModel.plot_fig1(concs[:, 3], C[:, :, 3, 9], D[:, :, 2], "palbo50 + gems", "total", "B")
    p3 = DrugResponseModel.plot_fig1(concs[:, 5], C[:, 3, :, 9], D[:, :, 3], "gem10 + palbos", "total", "C")
    p4 = DrugResponseModel.plot_fig1(concs[:, 1], C[:, :, 3, 2], D[:, :, 4], "gem10 + lpts", "total", "D")
    p5 = DrugResponseModel.plot_fig1(concs[:, 5], C[:, 4, :, 4], D[:, :, 5], "lpt100 + palbos", "total", "E")

    figureSS = plot(p1, p2, p3, p4, p5, size = (2000, 450), layout = (1, 5))
    savefig(figureSS, "figureSS.svg")
end
