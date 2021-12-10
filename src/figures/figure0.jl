""" plot the new data. """
function plot_data()
    g, c = DrugResponseModel.import_data()
    newg = DrugResponseModel.trim_data(g, c)
    tensor, conditions = DrugResponseModel.form_tensor(newg, c)
    tens = mean(tensor, dims=2)

    t = LinRange(0.0, 95, 189)
    p1 = scatter(t, tens[1, 1, :, :, 1], label = conditions[1], ylabel="G1 #")
    p2 = scatter(t, tens[2, 1, :, :, 1], label = conditions[1], ylabel="S/G2 #")

    p3 = scatter(t, tens[1, 1, :, :, 2], label = conditions[2], ylabel="G1 #")
    p4 = scatter(t, tens[2, 1, :, :, 2], label = conditions[2], ylabel="S/G2 #")

    p5 = scatter(t, tens[1, 1, :, :, 3], label = conditions[3], ylabel="G1 #")
    p6 = scatter(t, tens[2, 1, :, :, 3], label = conditions[3], ylabel="S/G2 #")

    p7 = scatter(t, tens[1, 1, :, :, 4], label = conditions[4], ylabel="G1 #")
    p8 = scatter(t, tens[2, 1, :, :, 4], label = conditions[4], ylabel="S/G2 #")

    p9 = scatter(t, tens[1, 1, :, :, 5], label = conditions[5], ylabel="G1 #")
    p10 = scatter(t, tens[2, 1, :, :,5], label = conditions[5], ylabel="S/G2 #")

    p11 = scatter(t, tens[1, 1, :, :, 6], label = conditions[6], ylabel="G1 #")
    p12 = scatter(t, tens[2, 1, :, :,6], label = conditions[6], ylabel="S/G2 #")

    fig1 = plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, size = (2000, 1300), layout = (3, 4))
    savefig(fig1, "figure0.svg")
end