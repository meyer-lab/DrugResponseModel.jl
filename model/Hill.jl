using LeastSquaresOptim, DifferentialEquations, DelayDiffEq, DiffEqBase
include("DDEmodel.jl")
""" 
        This file contains Hill function, residuals of Hill based off of DDE, and optimization of it.
"""


# hillParams should be something like [EC50, alpha_min, alpha_max, alpha_b, EC50, beta_min, beta_max, beta_b...]
hill(p, concentration) =  p[2] + ((p[3] - p[2]) / (1 + 10^((concentration - p[1])*p[4])))

function residHill(hillParams, concentrations, g1, g2, g1_0, g2_0)
    """ This functions takes in hill parameters for all the concentrations and calculates
    DDE parameters, passes them to residual function and based off of these, optimizes the model
    and estimates hill parameters. """
    residues = 0.0

    # EC50  for all of the trials is hillParams[1]
    for ii in 1:length(concentrations)
        alpha = hill(append!([hillParams[1]], hillParams[2:4]), concentrations[ii])
        beta = hill(append!([hillParams[1]], hillParams[5:7]), concentrations[ii])
        tau1 = hill(append!([hillParams[1]], hillParams[8:10]), concentrations[ii])
        tau2 = hill(append!([hillParams[1]], hillParams[11:13]), concentrations[ii])
        gamma1 = hill(append!([hillParams[1], hillParams[14]], [0, hillParams[15]]), concentrations[ii])
        gamma2 = hill(append!([hillParams[1], hillParams[16]], [0, hillParams[17]]), concentrations[ii])

        pp = [alpha, beta, tau1, tau2, gamma1, gamma2]
        residues += ddesolve(g1, g2, g1_0, g2_0, pp, ii)
    end 

    return residues
end


function optimize_hill(concentrations, guess, low, high, g1, g2, g1_0, g2_0)
    residue(hillParams) = residHill(hillParams, concentrations, g1, g2, g1_0, g2_0)
    results_hill = optimize(residue, guess, LevenbergMarquardt(), lower = low, upper = high)
    return results_hill
end