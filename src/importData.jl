"""
        Imports data works for the ODE model
"""
function get_data(path_g2::String, path_total::String; max = 189)
    # Import data all the trials for each drug
    perc = readdlm(path_g2, ','; skipstart = 1)
    total = readdlm(path_total, ','; skipstart = 1)

    # Clip to data of interest
    perc = convert(Array, perc[1:max, 2:9])
    total = convert(Array, total[1:max, 2:9])

    init_cells = 20.0

    # rescaling the experimental data assuming we have 20 initial cells for each trial
    gs = zeros(2, size(perc, 1), 8)
    population = init_cells * total
    gs[2, :, :] = 0.01 * population .* perc
    gs[1, :, :] = population - gs[2, :, :]

    # removing the peaks
    for i = 1:8
        gs[1, :, i] = svg_filter(gs[1, :, i])
        gs[2, :, i] = svg_filter(gs[2, :, i])
    end
    return gs
end

function import_combination(filename::String)
    basePath = joinpath(dirname(pathof(DrugResponseModel)), "..", "data")
    path_g2 = string("/", basePath, "/", filename, "_CYCLE.csv")
    path_total = string("/", basePath, "/", filename, "_CELL.csv")
    perc = readdlm(path_g2, ','; skipstart = 1)
    total = readdlm(path_total, ','; skipstart = 1)

    # Clip to data of interest
    perc = convert(Array{Float64, 2}, perc)
    total = convert(Array{Float64, 2}, total)
    # removing the peaks
    for i = 1:size(perc, 2)
        perc[:, i] = svg_filter(perc[:, i])
        perc[:, i] = svg_filter(perc[:, i])
        total[:, i] = svg_filter(total[:, i])
        total[:, i] = svg_filter(total[:, i])
    end
    init_cells = 20.0

    # rescaling the experimental data assuming we have 20 initial cells for each trial
    gs = zeros(3, size(perc, 1), size(perc, 2))
    population = init_cells * total
    gs[1, :, :] = 0.01 * population .* perc
    gs[2, :, :] = population - gs[1, :, :]
    gs[3, :, :] = gs[1, :, :] + gs[2, :, :]

    return gs[:, :, 2:25], gs[:, :, 27:50], perc, total
end

function setup_data(drug_name::String)
    basePath = joinpath(dirname(pathof(DrugResponseModel)), "..", "data")

    dfname = string(drug_name, ".csv")
    dfname2 = string(drug_name, "_pop.csv")

    if occursin("Lapatinib", drug_name)
        idx = 1
    elseif occursin("Doxorubicin", drug_name)
        idx = 2
    elseif occursin("Gemcitabine", drug_name)
        idx = 3
    elseif occursin("Paclitaxel", drug_name)
        idx = 4
    elseif occursin("Palbociclib", drug_name)
        idx = 1
    end

    #----------- import concentrations
    concentration = readdlm(joinpath(basePath, "concentrations.csv"), ','; skipstart = 1)
    conc_l = [Float64(concentration[idx, col]) for col = 2:9]
    conc_l[1] = 0.05

    #------------ import cell data
    gs = get_data(joinpath(basePath, dfname), joinpath(basePath, dfname2))

    return conc_l, gs[2, :, :], gs[1, :, :]
end

function load(max, repi)
    g1s = zeros(max, 8, 5)
    g2s = zeros(max, 8, 5)
    concentrations = zeros(8, 5)
    drugs = ["Lapatinib", "Doxorubicin", "Gemcitabine", "Paclitaxel", "Palbociclib"]

    for i = 1:5
        concentrations[:, i], g2s[:, :, i], g1s[:, :, i] = setup_data(string(drugs[i], repi))
    end

    return concentrations, g1s + g2s, g1s, g2s
end

struct SavitzkyGolayFilter{M,N} end
@generated function (::SavitzkyGolayFilter{M,N})(data::AbstractVector{T}) where {M, N, T}
            #Create Jacobian matrix
            J = zeros(2M+1, N+1)
            for i=1:2M+1, j=1:N+1
                J[i, j] = (i-M-1)^(j-1)
            end
            e₁ = zeros(N+1)
            e₁[1] = 1.0

            #Compute filter coefficients
            C = J' \ e₁

            #Evaluate filter on data matrix

            To = typeof(C[1] * one(T)) #Calculate type of output
            expr = quote
                n = size(data, 1)
                smoothed = zeros($To, n)
                @inbounds for i in eachindex(smoothed)
                    smoothed[i] += $(C[M+1])*data[i]
                end
                smoothed
            end

            for j=1:M
                insert!(expr.args[6].args[3].args[2].args, 1,
                    :(if i - $j ≥ 1
                        smoothed[i] += $(C[M+1-j])*data[i-$j]
                      end)
                )
                push!(expr.args[6].args[3].args[2].args,
                    :(if i + $j ≤ n
                        smoothed[i] += $(C[M+1+j])*data[i+$j]
                      end)
                )
            end

            return expr
end

svg_filter = SavitzkyGolayFilter{41, 3}()
