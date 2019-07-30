include("DDEmodel.jl")
""" 
        This file contains Hill function, residuals of Hill based off of DDE, and optimization of it.
"""


# hillParams should be something like [EC50, alpha_min, alpha_max, alpha_b, EC50, beta_min, beta_max, beta_b...]
hill(p, concentration) =  p[2] + ((p[3] - p[2]) / (1 + 10^((concentration - p[1])*p[4])))

function residHill(hillParams, concentrations, g1, g2)
    """ This functions takes in hill parameters for all the concentrations and calculates
    DDE parameters, passes them to residual function and based off of these, optimizes the model
    and estimates hill parameters. """

    num_concentration = length(concentrations)
    data_size = length(g1[:, 1])
    total_res = zeros(2, data_size*num_concentration)

#     EC50  for all of the trials is hillParams[1]
    for ii in 1:num_concentration
        alpha = hill(append!([hillParams[1]], hillParams[2:4]), concentrations[ii])
        beta = hill(append!([hillParams[1]], hillParams[5:7]), concentrations[ii])
        tau1 = hill(append!([hillParams[1]], hillParams[8:10]), concentrations[ii])
        tau2 = hill(append!([hillParams[1]], hillParams[11:13]), concentrations[ii])
        gamma1 = hill(append!([hillParams[1], hillParams[14]], [0, hillParams[15]]), concentrations[ii])
        gamma2 = hill(append!([hillParams[1], hillParams[16]], [0, hillParams[17]]), concentrations[ii])

        pp = [alpha, beta, tau1, tau2, gamma1, gamma2]
        
        residuals = resid(pp, ii, g1, g2)

        # To append all of the residuals of each trial to the end of the last one, finally
        # having a matrix of 2 by 192*8 : the first row contains the residuals for G1 for 
        # all concentrations and the second row contains the residuals for all concentrations for G2. 
        total_res[1, ((ii-1)*data_size +1):(ii*data_size)] = residuals[1, :]
        total_res[2, ((ii-1)*data_size +1):(ii*data_size)] = residuals[2, :]
    end 
    return total_res
end

function optimize_hill(concentrations, guess, low, high)
    residue(hillParams) =residHill(hillParams, concentrations, g1, g2)
    results_hill = optimize(residue, guess, LevenbergMarquardt(), lower = low, upper = high)
    return results_hill
end
