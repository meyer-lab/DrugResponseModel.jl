""" Plot the 6 drugs from HCC dataset. """
function figure5()
    setGadflyTheme()

    tens, names, concs = DrugResponseModel.hcc_all()

    pp = []
    for i = 1:size(tens)[4]
        push!(pp, DrugResponseModel.Eachdrug(tens[1, :, :, i], concs[i], "G1", names[i]))
        push!(pp, DrugResponseModel.Eachdrug(tens[2, :, :, i], concs[i], "S/G2", names[i]))
    end

    pl = plotGrid((3, 4), [pp...];)
    return draw(SVG("figure1.svg", 16inch, 12inch), pl)
end