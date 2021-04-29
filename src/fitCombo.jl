""" In this file the data will be organized for combinations to be fit to the model. """

########### 1. control, Palbociclib 50, Palbociclib 50 nM + lapatinibs [25 nM, 50 nM, 100 nM, 250 nM]:
########### 2. control, Palbociclib 50, Palbociclib 50 nM + gemcitabines [5 nM, 10 nM, 17 nM, 30 nM]
########### 3. control, Gemcitabine 10, Gemcitabine 10 nM + palbociclibs [25 nM, 50 nM, 100 nM, 250 nM]
########### 4. control, Gemcitabine 10, Gemcitabine 10 nM + lapatinibs [25 nM, 50 nM, 100 nM, 250 nM]
########### 5. control, Doxorubicin 20, Doxorubicin 20 nM + gemcitabines [5 nM, 10 nM, 17 nM, 30 nM]
########### 6. control, Lapatinib  100, Lapatinib  100 nM + palbociclibs [25 nM, 50 nM, 100 nM, 250 nM]

""" To fit combinations. """
function RESID_through_bliss(p, g1c, g2c)
    p_combo = my_helper(p)
    t = LinRange(0.0, 0.5 * size(g1c, 1), size(g1c, 1))
    res = 0
    for j=1:5
        for i=1:6
            res += predict(p_combo[:, i, j], p_combo[:, 1, j], t, g1c[:, i, j], g2c[:, i, j])[1]
        end
    end
    res
end

function optimize_combo(g1::Array{Float64, 3}, g2::Array{Float64, 3}; maxstep = 300000)

    f(x) = RESID_through_bliss(x, g1, g2)

    # [EC50, k, max_a1, max_a2, max_b1,  max_b2, max_b3, max_b4, max_g11, max_g12, max_g21, max_g22, max_g23, max_g24, min_a1, min_a2, min_b1,  min_b2, min_b3, min_b4]
    low = [5.0; 0.1; 1e-5 * ones(12); 1.0; 0.1; 1e-5 * ones(12); 5.0; 0.1; 1e-5 * ones(18)] # 60 parameters including 12 for single drug A, and 20 for the drugB
    high = [500.0; 10.0; 2.5 * ones(12); 50.0; 10.0; 2.5 * ones(12); 500.0; 10.0; 2.5 * ones(18)]

    return optimize_helper(f, low, high, maxstep)
end

""" Unit function to calculate the bliss for 2 drugs at one specific concentration. """
function Bliss_unit(pp1, pp2, control)
    # pp1 and pp2 are 1D arrays of size 8, including 8 parameters for a single concentration.
    p1 = copy(pp1)
    p2 = copy(pp2)
    # normalization
    p1[1:6] .= 1.0 .- (pp1[1:6] ./ control[1:6]) # g1 and g2 prog. rates
    p1[7:12] .= pp1[7:12]                          # g1 and g2 death rates
    # drug B
    p2[1:6] .= 1.0 .- (pp2[1:6] ./ control[1:6])
    p2[7:12] .= pp2[7:12]

    c = Array{eltype(pp1), 1}(undef, 12)
    c[1:6] .= (1.0 .- (p1[1:6] .+ p2[1:6] .- p1[1:6] .* p2[1:6])) .* control[1:6]
    c[7:12] .= p1[7:12] .+ p2[7:12]

    c
end


function my_helper(p)

    # [12p_DOX20, 14p_lap, 14p_gem, 14p_palbo, 6p_control]
    concs = hcat([0, 25, 50, 100, 250], [0, 5, 10, 17, 30], [0, 25, 50, 100, 250])
    dox = p[1:12]
    p_ode = getODEparams(p, concs) # returns 12 x 5 x 3
    p_combo = zeros(12, 6, 5)
    p_combo[:, 1, :] .= p_ode[:, 1, 1] # control

    for i=2:6
        p_combo[:, i, 1] = Bliss_unit(p_ode[:, 3, 3], p_ode[:, i-1, 1], p_ode[:, 1, 1])
        p_combo[:, i, 2] = Bliss_unit(p_ode[:, 3, 3], p_ode[:, i-1, 2], p_ode[:, 1, 1])
        p_combo[:, i, 3] = Bliss_unit(p_ode[:, 2, 2], p_ode[:, i-1, 3], p_ode[:, 1, 1])
        p_combo[:, i, 4] = Bliss_unit(p_ode[:, 2, 2], p_ode[:, i-1, 1], p_ode[:, 1, 1])
        # p_combo[:, i, 5] = Bliss_unit(dox, p_ode[:, i-1, 2], p_ode[:, 1, 1])
        p_combo[:, i, 5] = Bliss_unit(p_ode[:, 4, 1], p_ode[:, i-1, 3], p_ode[:, 1, 1])
    end
    p_combo
end
