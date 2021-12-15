""" Plotting the model and data, fitted separately. """
global drugs = ["BEZ235", "Trametinib", "5FU", "AZD5438", "Panobinostat", "MG132", "Everolimus", "JQ1", "Bortezomib", "MK1775", "Cabozantinib"]

t = LinRange(0.0, 95, 189)

""" Takes the above parameters and returns the inferred ODE params. The output is used in figure3. """
function out_ODEparams()
    p1 = [4.559, 1.534, 6.577953772895624e-6, 0.060, 0.0876, 0.0002, 9.9002, 0.389, 9.959, 0.573, 0.0045, 3.11e-5, 8.92e-5, 0.0336, 0.500, 4.599678610545463e-5, 0.0006, 0.080, 0.056, 0.375, 0.281, 0.313, 9.914, 0.223, 9.958, 9.993];
    p2 = [16.341, 2.0267, 0.0697, 0.037, 4.732, 7.917, 0.333, 0.562, 9.88, 9.554, 0.0033, 6.14e-5, 0.002, 0.035, 0.0041, 0.00058, 0.00026, 0.00073, 9.661, 0.0411, 6.869, 9.483, 1.438, 0.250, 9.995, 9.580];
    p3 = [45.303, 1.228, 4.586, 0.082, 0.0309, 0.354, 4.793, 0.0916, 4.898, 0.486, 0.00022, 3.3071e-5, 1.821e-5, 4.13e-5, 3.126e-5, 7.669e-5, 0.0003, 1.25e-5, 9.1008, 0.0644, 0.8656, 0.272, 9.0082, 3.243, 5.628, 0.2105];
    p4 = [465.346, 2.833, 0.088, 0.00032, 0.0044, 0.0012, 0.339, 0.1199, 1.0501, 0.528, 0.0018, 0.0102, 0.00017, 0.00031, 0.00034, 0.00024, 0.0091, 0.0865, 0.0494, 1.557, 4.212, 0.7538, 7.355, 2.842, 6.619, 0.283];
    p5 = [4.9001, 3.907, 7.241, 3.13, 0.053, 0.138, 0.050, 9.502, 6.293, 0.357, 0.0053, 0.000525, 4.81e-5, 2.249, 0.0023, 0.137, 0.0044, 2.99e-5, 8.0712, 3.319, 0.046, 8.142, 0.712, 9.906, 7.0182, 0.368];
    p6 = [20.424, 7.59, 4.155e-5, 9.445, 0.0267, 1.0941e-5, 6.056, 7.66, 0.537, 9.998, 0.0654, 0.000151, 1.22e-6, 0.165, 0.0138, 1.123, 8.911e-8, 2.645e-5, 9.99, 9.752, 0.0357, 1.931, 9.974, 9.993, 0.219, 9.999];
    p7 = [0.0049, 0.108, 9.655, 0.206, 0.0076, 0.292, 9.255, 1.156, 0.327, 9.827, 0.000265, 0.004, 1.181e-6, 3.801e-5, 0.0093, 0.00047, 0.0016, 3.12e-5, 9.958, 0.0497, 0.165, 0.536, 9.289, 1.1475, 0.2946, 9.909];
    p8 = [298.17, 1.043, 0.00029, 0.019, 7.128, 0.683, 1.36, 0.449, 0.00125, 7.388, 0.0266, 0.000679, 0.542, 0.0016, 0.00049, 0.00354, 0.09281, 0.00769, 7.4303, 0.0408, 7.322, 4.69, 9.589, 0.2284, 9.976, 7.9261];
    p9 = [0.361, 0.439, 0.114, 0.564, 0.0026, 0.049, 3.676, 0.401, 0.787, 9.586, 3.334e-6, 0.000269, 0.127, 0.000439, 9.812, 9.916, 4.574e-5, 0.000555, 0.042429, 2.3973, 9.9912, 0.214, 4.454, 0.399, 0.6626, 9.602];
    p10 = [1468.105, 2.062, 0.00213, 0.121, 0.109, 0.6952, 0.515, 0.262, 0.01049, 0.127, 0.0198, 9.414e-5, 0.00139, 0.0292, 5.94, 2.732, 0.00852, 0.000142, 0.58902, 0.0473, 0.176, 3.393, 8.863, 2.21, 1.5485, 0.366];
    p11 = [3307.225, 2.825, 4.261, 0.1908, 0.0002, 6.151, 0.6462, 0.00569, 0.311, 0.0768, 0.0199, 2.106e-5, 0.0030, 0.00519, 0.0907, 0.00159, 0.000354, 0.0095, 9.792, 0.0706, 0.13645, 9.82027, 8.371, 0.828, 0.2809, 7.846];
    parameters = hcat(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11)
    pp1 = zeros(11, 16, 8, 1) # 16 params, 8 concentrations, 11 drugs: the value is the effect at concentration 5th
    for (i, drug) in enumerate(drugs)
        _, _, conc = DrugResponseModel.load_newData(drug)
        pp1[i, :, :, 1] = DrugResponseModel.getODEparams(parameters[:, i], conc)
    end
    return pp1[:, :, :, 1]
end

function Eachdrug_sim(G_sim, G_data, condition, g1g2, drug_name)
    df = [DataFrame(x=t, y=G_sim[:, i], y2=G_data[:, i], conc=condition[i]) for i = 1:8]
    DF = vcat(df...)
    p = Gadfly.plot(DF, layer(x="x", y="y", Geom.line, color="conc"), 
                        layer(x="x", y="y2", Geom.line, color="conc", Theme(line_style=[:dash])),
                        Guide.xlabel("time [hr]"), Guide.ylabel("$g1g2 Cell #"), Coord.Cartesian(ymin=-0.05,ymax=2.0), Guide.title(drug_name))
    return p
end

function figure2()
    setGadflyTheme()

    g, c = DrugResponseModel.import_data()
    newg = DrugResponseModel.trim_data(g, c)
    tensor, conditions = DrugResponseModel.form_tensor(newg, c)
    tens = mean(tensor, dims=2)[:, 1, :, :, :]
    pp = []
    pODE = out_ODEparams()
    for (i, drug) in enumerate(drugs)
        _, _, conc = DrugResponseModel.load_newData(drug)
        Gs = zeros(189, 8, 2)
        for j=1:8
            Gs[:, j, 1], Gs[:, j, 2], _ = DrugResponseModel.predict(pODE[i, :, j], pODE[i, :, 1], t)
        end
        push!(pp, Eachdrug_sim(Gs[:, :, 1], tens[1, :, :, i], conditions[i], "G1", drug))
        push!(pp, Eachdrug_sim(Gs[:, :, 2], tens[2, :, :, i], conditions[i], "S/G2", drug))
    end

    pl = plotGrid((6, 4), [pp..., nothing, nothing];)
    return draw(SVG("figure2.svg", 22inch, 18inch), pl)
end
