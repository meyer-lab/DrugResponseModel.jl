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

        # collecting all the DDE model parameters
        pp = [alpha, beta, tau1, tau2, gamma1, gamma2]
        # calculating the residulas for this set of parameters
        residues += ddesolve(g1, g2, g1_0, g2_0, pp, ii)
    end 
    return residues
end

# optimization function for the Hill model
function optimize_hill(hillParams, concentrations, g1, g2, g1_0, g2_0, low, high, num_steps)
    # changing the objective function to be compatible with bboptimize
    residue(hillParams) = residHill(hillParams, concentrations, g1, g2, g1_0, g2_0)
    res = bboptimize(residue; SearchRange=collect(zip(low, high)), TraceMode=:compact, MaxSteps=num_steps, Method = :adaptive_de_rand_1_bin_radiuslimited)
    return res, best_candidate(res)

end

function getDDEparams(p, concentrations)
    effects = zeros(6, 8)
    for i in 1:8
        effects[1, i] = hill(p[1:4], concentrations[i])
        effects[2, i] = hill(append!([p[1]], p[5:7]), concentrations[i])
        effects[3, i] = hill(append!([p[1]], p[8:10]), concentrations[i])
        effects[4, i] = hill(append!([p[1]], p[11:13]), concentrations[i])
        effects[5, i] = hill(append!([p[1], p[14]], [0, p[15]]), concentrations[i])
        effects[6, i] = hill(append!([p[1], p[16]], [0, p[17]]), concentrations[i])
    end
    return effects
end
