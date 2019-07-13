include("DDEmodel.jl")

hill(p, concentration) =  p[1] + ((p[2] - p[1]) / (1 + 10^((concentration - p[3])*p[4])))

function residHill(hillParams, concentrations, g1, g2)
    
    concentration_size = length(concentrations)
    data_size = length(g1[:, 1])
    total = zeros(2, data_size*concentration_size)
    for ii in 1:concentration_size
        alpha = hill(hillParams[1:4], concentrations[ii])
        beta = hill(hillParams[5:8], concentrations[ii])
        tau1 = hill(hillParams[9:12], concentrations[ii])
        tau2 = hill(hillParams[13:16], concentrations[ii])
        history = hill(hillParams[17:20], concentrations[ii])
        gamma1 = hill(hillParams[21:24], concentrations[ii])
        gamma2 = hill(hillParams[25:28], concentrations[ii])

        pp = [alpha, beta, tau1, tau2, history, gamma1, gamma2]
        residuals = resid(pp, ii, g1, g2)
        total[1, ((ii-1)*data_size +1):(ii*data_size)] = residuals[1, :]
        total[2, ((ii-1)*data_size +1):(ii*data_size)] = residuals[2, :]
    end 
    return total
end


function optimize_hill(concentrations, guess, low, high)
    residue(hillParams) =residHill(hillParams, concentrations, g1, g2)
    results_hill = optimize(residue, guess, LevenbergMarquardt(), lower = low, upper = high)
    return results_hill
end