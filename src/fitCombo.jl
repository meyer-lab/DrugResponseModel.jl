""" In this file the data will be organized for combinations to be fit to the model. """

########### 1. control, Palbociclib 50, Palbociclib 50 nM + lapatinibs [25 nM, 50 nM, 100 nM, 250 nM]:
# p = [0.189371, 1.31225, 0.331636, 3.0, 3.0, 0.274225, 3.32364e-9, 1.48818e-9, 1.5106e-9, 2.04528e-9, 2.3648e-8, 0.0132698, 33.8183, 1.15343, 1.16521e-8, 0.0249639, 0.392687, 0.402105, 2.3868e-8, 0.338073, 0.00245421, 1.36638e-8, 1.62531e-9, 9.78796e-8, 0.0592202, 0.0143732, 0.228744, 2.5, 1.23386, 0.347899, 2.49999, 0.210224]
########### 2. control, Palbociclib 50, Palbociclib 50 nM + gemcitabines [5 nM, 10 nM, 17 nM, 30 nM]
########### 3. control, Gemcitabine 10, Gemcitabine 10 nM + palbociclibs [25 nM, 100 nM, 250 nM, 250 nM]
########### 4. control, Gemcitabine 10, Gemcitabine 10 nM + lapatinibs [25 nM, 50 nM, 100 nM, 250 nM]
########### 5. control, Doxorubicin 20, Doxorubicin 20 nM + gemcitabines [5 nM, 10 nM, 17 nM, 30 nM]
########### 6. control, Lapatinib  100, Lapatinib  100 nM + palbociclibs [25 nM, 100 nM, 250 nM, 250 nM]

""" To fit combinations. """
function RESID_through_bliss(p, concs, g1c, g2c)
    p_combo = my_helper(p, concs)
    t = LinRange(0.0, 0.5 * size(g1c, 1), size(g1c, 1))
    res = 0
    for i=1:6 # 
        res += predict(p_combo[:, i], p_combo[:, 1], t, g1c[:, i], g2c[:, i])[1]
    end
    res
end

function optimize_combo(conc::Vector, g1::Matrix, g2::Matrix; maxstep = 300000)

    f(x) = RESID_through_bliss(x, conc, g1, g2)

    # [EC50, k, max_a1, max_a2, max_b1,  max_b2, max_b3, max_b4, max_g11, max_g12, max_g21, max_g22, max_g23, max_g24, min_a1, min_a2, min_b1,  min_b2, min_b3, min_b4]
    low = [1e-9 * ones(12); minimum(conc); 1e-9 * ones(19)] # 32 parameters including 12 for single drug A, and 20 for the drugB
    high = [3 * ones(12); maximum(conc); 10.0; 2.5 * ones(18)]

    return optimize_helper(f, low, high, maxstep)
end

function my_helper(p, concs)

    p_ode = getODEparams(p[13:32], concs)
    p_combo = zeros(12, 6)
    p_combo[:, 1] = p_ode[:, 1, 1] # control

    for i=2:6
        p_combo[:, i] = DrugResponseModel.Bliss_params_unit(p[1:12], p_ode[:, i-1, 1], hcat(p_ode[:, 1, 1], p_ode[:, 1, 1]))
    end
    p_combo
end
