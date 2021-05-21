""" In this file the data will be organized for combinations to be fit to the model. """

### 1. control, Palbociclib 50, Palbociclib 50 nM + lapatinibs [25 nM, 50 nM, 100 nM, 250 nM]:
### 2. control, Palbociclib 50, Palbociclib 50 nM + gemcitabines [5 nM, 10 nM, 17 nM, 30 nM]
### 3. control, Gemcitabine 10, Gemcitabine 10 nM + palbociclibs [25 nM, 50 nM, 100 nM, 250 nM]
### 4. control, Gemcitabine 10, Gemcitabine 10 nM + lapatinibs [25 nM, 50 nM, 100 nM, 250 nM]
### 5. control, Lapatinib  100, Lapatinib  100 nM + palbociclibs [25 nM, 50 nM, 100 nM, 250 nM]
### 6. control, Lapatinib  100, Lapatinib  100 nM + gemcitabines [5 nM, 10 nM, 17 nM, 30 nM]
### 7. control, Doxorubicin 20, Doxorubicin 20 nM + gemcitabines [5 nM, 10 nM, 17 nM, 30 nM]
### 8. control, Paclitaxel 2, paclitaxel 2 nM + [palbo50, dox20, lap50, lap100, gem10 nM]

""" The residuals function that calculates cost function for combination data vs model. """
function RESID_through_bliss(p, g1c, g2c)
    t = LinRange(0.0, 0.5 * size(g1c, 1), size(g1c, 1))
    res = 0
    res += 1e4 * (maximum([0, mean(p[1:2]) - mean(p[67:68])]))^2
    res += 1e4 * (maximum([0, mean(p[3:6]) - mean(p[69:72])]))^2
    res += 1e4 * (maximum([0, mean(p[13:14]) - mean(p[67:68])]))^2
    res += 1e4 * (maximum([0, mean(p[15:18]) - mean(p[69:72])]))^2
    res += 1e4 * (maximum([0, mean(p[27:28]) - mean(p[67:68])]))^2
    res += 1e4 * (maximum([0, mean(p[29:32]) - mean(p[69:72])]))^2
    res += 1e4 * (maximum([0, mean(p[41:42]) - mean(p[67:68])]))^2
    res += 1e4 * (maximum([0, mean(p[43:46]) - mean(p[69:72])]))^2
    res += 1e4 * (maximum([0, mean(p[55:56]) - mean(p[67:68])]))^2
    res += 1e4 * (maximum([0, mean(p[57:60]) - mean(p[69:72])]))^2

    p_combo = my_helper(p)
    for j=1:8
        for i=1:6
            res += predict(p_combo[:, i, j], p_combo[:, 1, j], t, g1c[:, i, j], g2c[:, i, j])[1]
        end
    end
    res
end

""" Optimizer function that we run for fitting. """
function optimize_combo(g1::Array{Float64, 3}, g2::Array{Float64, 3}; maxstep = 700000)

    f(x) = RESID_through_bliss(x, g1, g2)

    # [EC50, k, max_a1, max_a2, max_b1,  max_b2, max_b3, max_b4, max_g11, max_g12, max_g21, max_g22, max_g23, max_g24, min_a1, min_a2, min_b1,  min_b2, min_b3, min_b4]
    low = [0.01 * ones(12); 0.01 * ones(12); 5.0; 0.1; 0.01 * ones(12); 1.0; 0.1; 0.01 * ones(12); 5.0; 0.1; 0.01 * ones(12); 0.2 * ones(6)] # 60 parameters including 12 for single drug A, and 20 for the drugB
    high = [2 * ones(12); 1.8 * ones(12); 500.0; 100.0; 1.8 * ones(12); 200.0; 20.0; 1.8 * ones(12); 500.0; 10.0; 2 * ones(12); 3.0 * ones(6)]

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

    @assert(all(c[7:12] .>= p1[7:12]))
    @assert(all(c[7:12] .>= p2[7:12]))
    c
end

""" This function converts estimated mix parameters into ODE parameters. """
function my_helper(p)

    # [12p_DOX20, 14p_lap, 14p_gem, 14p_palbo, 6p_control]
    concs = hcat([0, 25, 50, 100, 250], [0, 5, 10, 17, 30], [0, 25, 50, 100, 250])
    p_ode = zeros(12, 5, 5)
    p_ode[:, 2, 2] = p[1:12]
    p_ode[:, 2, 4] = p[13:24]
    p_ode[:, :, [1, 3, 5]] = getODEparams(p[25:end], concs) # returns 12 x 5 x 3
    p_combo = zeros(12, 6, 8)
    p_combo[:, 1, :] .= p_ode[:, 1, 1] # control

    for i=2:6
        p_combo[:, i, 1] = Bliss_unit(p_ode[:, 3, 5], p_ode[:, i-1, 1], p_ode[:, 1, 1])
        p_combo[:, i, 2] = Bliss_unit(p_ode[:, 3, 5], p_ode[:, i-1, 3], p_ode[:, 1, 1])
        p_combo[:, i, 3] = Bliss_unit(p_ode[:, 3, 3], p_ode[:, i-1, 5], p_ode[:, 1, 1])
        p_combo[:, i, 4] = Bliss_unit(p_ode[:, 3, 1], p_ode[:, i-1, 1], p_ode[:, 1, 1])
        p_combo[:, i, 5] = Bliss_unit(p_ode[:, 4, 1], p_ode[:, i-1, 5], p_ode[:, 1, 1])
        p_combo[:, i, 6] = Bliss_unit(p_ode[:, 4, 1], p_ode[:, i-1, 3], p_ode[:, 1, 1])
        p_combo[:, i, 7] = Bliss_unit(p_ode[:, 2, 2], p_ode[:, i-1, 3], p_ode[:, 1, 1])
    end
    p_combo[:, 2, 8] = Bliss_unit(p_ode[:, 2, 4], p_ode[:, 3, 5], p_ode[:, 1, 1]) # pax2 + palbo50
    p_combo[:, 3, 8] = Bliss_unit(p_ode[:, 2, 4], p_ode[:, 2, 2], p_ode[:, 1, 1]) # pax2 + 12p_DOX20
    p_combo[:, 4:5, 8] .= Bliss_unit(p_ode[:, 2, 4], p_ode[:, 3:4, 1], p_ode[:, 1, 1]) # pax2 + lap50, lap100
    p_combo[:, 6, 8] = Bliss_unit(p_ode[:, 2, 4], p_ode[:, 3, 3], p_ode[:, 1, 1]) # pax2 + gem10
    p_combo
end
