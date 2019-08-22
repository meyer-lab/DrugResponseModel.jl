include("DDEmodel.jl")
""" 
        This file contains Hill function, residuals of Hill based off of DDE, and optimization of it.
"""


# hillParams should be something like [EC50, min, max, b]
hill(p, concentration) =  p[2] + ((p[3] - p[2]) / (1 + 10^((concentration - p[1])*p[4])))

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
        tau1 = hill(append!([hillParams[1], hillParams[7], hillParams[8]], [hillParams[2]]), concentrations[ii])
        tau2 = hill(append!([hillParams[1], hillParams[9], hillParams[10]], [hillParams[2]]), concentrations[ii])
        gamma1 = hill(append!([hillParams[1], hillParams[11], 0.0], [hillParams[2]]), concentrations[ii])
        gamma2 = hill(append!([hillParams[1], hillParams[12], 0.0], [hillParams[2]]), concentrations[ii])

        # collecting all the DDE model parameters
        pp = [alpha, beta, tau1, tau2, gamma1, gamma2]
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
    low = [30.0, 0.0001, 0.006, 0.001, 0.02, 0.04, 6.0, 30.0, 8.0, 6.0, 0.0004, 0.03]
    # upper bound
    high = [250.0, 0.01, 0.008, 0.03, 0.08, 0.05, 10.0, 40.0, 15.0, 10.0, 0.003, 0.05]
    res = bboptimize(residue; SearchRange=collect(zip(low, high)), TraceMode=:compact, MaxSteps=num_steps, TraceInterval=50, Method = :adaptive_de_rand_1_bin_radiuslimited)
    return res, best_candidate(res)

end

function getDDEparams(p, concentrations)
    effects = zeros(6, 8)
    for i in 1:8
        effects[1, i] = hill(append!([p[1], p[3]], [p[4], p[2]]), concentrations[i])
        effects[2, i] = hill(append!([p[1], p[5]], [p[6], p[2]]), concentrations[i])
        effects[3, i] = hill(append!([p[1], p[7]], [p[8], p[2]]), concentrations[i])
        effects[4, i] = hill(append!([p[1], p[9]], [p[10], p[2]]), concentrations[i])
        effects[5, i] = hill(append!([p[1], p[11]], [0, p[2]]), concentrations[i])
        effects[6, i] = hill(append!([p[1], p[12]], [0, p[2]]), concentrations[i])
    end
    return effects
end

### to fit Hill curve to the estimated rates
# this fits a single parameter for all four drugs, so for every single parameters, i.e., alpha, beta, tau1, tau2, gamma, we will need one call of this function to get all the plots for all our drugs.
# k is an integer $\in$ {1, 2, 3, 4, 5, 6} representing {alpha, beta, tau1, tau2, gamma1, gamma2}, respctively.
function fit_hill(parameters, k, conc_l, conc_d, conc_g, conc_t)
    p0 = [100.0, 1.0, 0.007, 0.01]
    fit_l = curve_fit(hill, conc_l, parameters[1, k, :], p0)
    fit_d = curve_fit(hill, conc_d, parameters[2, k, :], p0)
    fit_g = curve_fit(hill, conc_g, parameters[3, k, :], p0)
    fit_t = curve_fit(hill, conc_t, parameters[4, k, :], p0)
    return (coef(fit_l), coef(fit_d), coef(fit_g), coef(fit_t))
end
