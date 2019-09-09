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
    low = [50.0, 0.1, 0.005, 0.04, 0.05, 0.01, 26.0, 5.0, 10.0, 7.0, 0.01, 0.1]
    # upper bound
    high = [250, 100.0, 0.02, 0.1, 0.2, 0.03, 36.0, 10.0, 16.0, 11.0, 0.04, 1.0]

    println("global optimization begins ...")
    res = bboptimize(residue; SearchRange=collect(zip(low, high)), MaxSteps=num_steps, TraceInterval=50, Method =:adaptive_de_rand_1_bin_radiuslimited)
    new_guess = best_candidate(res)
    optim_res = optimize(residue, low, high, new_guess, Fminbox(), Optim.Options(show_trace=true))
    return Optim.minimum(optim_res), Optim.minimizer(optim_res)
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
