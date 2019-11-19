include("DDEmodel.jl")

""" 
        This file contains Hill function, residuals of Hill based off of DDE, and optimization of it.
"""

hill(p::Array{Float64}, concentration::Float64) = p[2] + ((p[3]-p[2])/(1 + ((p[1])/(concentration))^p[4]))

function residHill(hillParams::Array{Float64}, concentrations::Array{Float64}, g1::Matrix{Float64}, g2::Matrix{Float64}, g1_0::Array{Float64}, g2_0::Array{Float64})
    """ This functions takes in hill parameters for all the concentrations and calculates
    DDE parameters, passes them to residual function and based off of these, optimizes the model
    and estimates hill parameters. """
    residues = 0.0
    times = range(0.0; stop = 95.5, length = 192)

    # EC50  for all of the trials is hillParams[1]
    # sharpness parameter (b) for all the trials is hillParams[2]
    # the two in-between-parameters are the min and max of the hill curve
    # for gamma1 and gamma2, the min has been set to 0, so totally we have 12 parameters to estimate
    for ii in 1:length(concentrations)
        alpha = hill(append!([hillParams[1], hillParams[4], hillParams[3]], [hillParams[2]]), concentrations[ii])
        beta = hill(append!([hillParams[1], hillParams[6], hillParams[5]], [hillParams[2]]), concentrations[ii])
        tau1 = hill(append!([hillParams[1], hillParams[8], hillParams[7]], [hillParams[2]]), concentrations[ii])
        tau2 = hill(append!([hillParams[1], hillParams[10], hillParams[9]], [hillParams[2]]), concentrations[ii])
        gamma1 = hill(append!([hillParams[11], 0.0, hillParams[12]], [hillParams[2]]), concentrations[ii])
        gamma2 = hill(append!([hillParams[13], 0.0, hillParams[14]], [hillParams[2]]), concentrations[ii])

        # collecting all the DDE model parameters
        pp = log.([alpha, beta, tau1, tau2, gamma1, gamma2])
        # calculating the residulas for this set of parameters
        alg, prob, data = ddesolve(collect(times), g1, g2, g1_0, g2_0, pp, ii)
        obj = build_loss_objective(prob, alg, L2Loss(times, data);
                               prob_generator=prob_generator,
                               verbose_opt=false)
        res = obj(pp)
        residues += res
    end 
    return residues
end

# optimization function for the Hill model
function optimize_hill(guess::Array{Float64}, concentrations::Array{Float64}, g1::Matrix{Float64}, g2::Matrix{Float64}, g1_0::Array{Float64}, g2_0::Array{Float64}, num_steps::Int, lowEC50::Float64, highEC50::Float64)
    """ A function to do the optimization given the lower and upper bound of estimation space. """
    # changing the objective function to be compatible with bboptimize
    residue(hillParams) = residHill(hillParams, concentrations, g1, g2, g1_0, g2_0)
    # lower bound
    low = [lowEC50, 0.1, 0.005, 0.04, 0.005, 0.01, 6.0, 5.0, 6.0, 5.0, lowEC50, 0.00001, lowEC50, 0.00001]
    # upper bound
    high = [highEC50, 10.0, 0.1, 0.1, 0.2, 0.03, 30.0, 25.0, 30.0, 25.0, highEC50, 0.05, highEC50, 0.01]

    res = bboptimize(residue; SearchRange=collect(zip(low, high)), MaxSteps=num_steps, TraceInterval=100, Method =:adaptive_de_rand_1_bin_radiuslimited)
    new_guess = best_candidate(res)
    return best_fitness(res), new_guess
end

function getDDEparams(p::Array{Float64}, concentrations::Array{Float64})
    """ A function to convert the estimated hill parameters back to DDE parameters. """
    effects = zeros(6, 8)
    for i in 1:8
        effects[1, i] = hill(append!([p[1], p[4]], [p[3], p[2]]), concentrations[i])
        effects[2, i] = hill(append!([p[1], p[6]], [p[5], p[2]]), concentrations[i])
        effects[3, i] = hill(append!([p[1], p[8]], [p[7], p[2]]), concentrations[i])
        effects[4, i] = hill(append!([p[1], p[10]], [p[9], p[2]]), concentrations[i])
        effects[5, i] = hill(append!([p[11], 0.0, p[12]], [p[2]]), concentrations[i])
        effects[6, i] = hill(append!([p[13], 0.0, p[14]], [p[2]]), concentrations[i])
    end
    return effects
end

function ParamForBliss(p)
    """ To calculate Bliss independence drug effect
    we assume delays are constant, death rates are additive,
    and will keep the alpha and beta intact."""
    par = zeros(3,8)
    par[1,:] = p[1,:] # alpha stays the same
    par[2,:] = p[2,:] # beta stays the same
    par[3,:] = p[5,:] + p[6,:] # additivity assumption for death rates
    return par
end

function BlissCombination(p1::Matrix{Float64}, p2::Matrix{Float64})
    """ A function to calculate Bliss independence for drug combination assuming
    the two drugs hit different pathways and they effect completely independently. """

    param1 = ParamForBliss(p1)
    param2 = ParamForBliss(p2)
    Effect = zeros(8,8,3)
    for j in 1:8
        for k in 1:8
            Effect[j,k,:] .= param1[:,j] .+ param2[:,k] .- param1[:,j] .* param2[:,k]
            end
        end
    return Effect
end

function DDEcombinationParam(Effects, p1::Matrix, p2::Matrix)
   """ To put together parameters of the combination, into the DDE parameter format, 
    To be able to plot the number of cells over time with the combined effect. """
    ddeParam = zeros(8,8,5)
    ddeParam[:,:,1:2] = Effects[:,:,1:2] # alpha and beta
    ddeParam[:,:,3] = ((p1[3,4] + p2[3,4])/2)*ones(8,8) # I'm keeping the forth delay as the consant one
    ddeParam[:,:,4] = ((p1[4,4] + p2[4,4])/2)*ones(8,8) # I'm keeping the forth delay as the consant one
    ddeParam[:,:,5] = Effects[:,:,3] # death rate
    return ddeParam
end

