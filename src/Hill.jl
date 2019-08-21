include("DDEmodel.jl")
""" 
        This file contains Hill function, residuals of Hill based off of DDE, and optimization of it.
"""


# hillParams should be something like [EC50, min, max, b]
hill(p, concentration) =  p[2] + ((p[3] - 1*p[2]) / (1 + 10^((concentration - 1*p[1])*p[4])))

function residHill(hillParams, concentrations, g1, g2, g1_0, g2_0)
    """ This functions takes in hill parameters for all the concentrations and calculates
    DDE parameters, passes them to residual function and based off of these, optimizes the model
    and estimates hill parameters. """
    residues = 0.0

    # EC50  for all of the trials is hillParams[1]
    # sharpness parameter (b) for all the trials is hillParams[2]
    # the two in-between-parameters are the min and max of the hill curve
    # for gamma1 and gamma2, the min has been set to 0, so totally we have 12 parameters to estimate
    for ii in 1:length(concentrations)
        alpha = hill(append!([hillParams[1], hillParams[3], hillParams[4]], [hillParams[2]]), concentrations[ii])
        beta = hill(append!([hillParams[1], hillParams[5], hillParams[6]], [hillParams[2]]), concentrations[ii])
#         tau1 = hill(append!([hillParams[1], hillParams[7], hillParams[8]], [hillParams[2]]), concentrations[ii])
#         tau2 = hill(append!([hillParams[1], hillParams[9], hillParams[10]], [hillParams[2]]), concentrations[ii])
        gamma1 = hill(append!([hillParams[1], hillParams[7], 0.0], [hillParams[2]]), concentrations[ii])
        gamma2 = hill(append!([hillParams[1], hillParams[8], 0.0], [hillParams[2]]), concentrations[ii])
        coef = hillParams[9]

        # collecting all the DDE model parameters
        pp = [alpha, beta, gamma1, gamma2, coef]
        # calculating the residulas for this set of parameters
        residues += ddesolve(g1, g2, g1_0, g2_0, pp, ii)
    end 
    return residues
end

# optimization function for the Hill model
function optimize_hill(guess, concentrations, g1, g2, g1_0, g2_0, num_steps)
    # changing the objective function to be compatible with bboptimize
    residue(hillParams) = residHill(hillParams, concentrations, g1, g2, g1_0, g2_0)
    # lower bound
#     low = [30.0, 0.01, 0.006, 0.001, 0.02, 0.04, 6.0, 30.0, 8.0, 6.0, 0.0004, 0.03]
#     # upper bound
#     high = [250.0, 10.0, 0.008, 0.03, 0.08, 0.05, 10.0, 40.0, 15.0, 10.0, 0.003, 0.05]
    # lower bound
    low = [30.0, 0.01, 0.006, 0.001, 0.02, 0.04, 0.0004, 0.03, 0.0001]
    # upper bound
    high = [300.0, 10.0, 0.008, 0.03, 0.08, 0.05, 0.003, 0.05, 1.0]
    res = bboptimize(residue; SearchRange=collect(zip(low, high)), TraceMode=:compact, MaxSteps=num_steps, TraceInterval=50, Method = :adaptive_de_rand_1_bin_radiuslimited)
    return res, best_candidate(res)

end

function getDDEparams(p, concentrations)
    effects = zeros(5, 8)
    for i in 1:8
        effects[1, i] = hill(append!([p[1], p[3]], [p[4], p[2]]), concentrations[i])
        effects[2, i] = hill(append!([p[1], p[5]], [p[6], p[2]]), concentrations[i])
#         effects[3, i] = hill(append!([p[1], p[7]], [p[8], p[2]]), concentrations[i])
#         effects[4, i] = hill(append!([p[1], p[9]], [p[10], p[2]]), concentrations[i])
        effects[3, i] = hill(append!([p[1], p[7]], [0, p[2]]), concentrations[i])
        effects[4, i] = hill(append!([p[1], p[8]], [0, p[2]]), concentrations[i])
        effects[5, i] = p[9]
    end
    return effects
end
