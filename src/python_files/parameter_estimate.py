""" Here is a function to estimate the Gamma distribution parameters for a given data set."""

import numpy as np


def GAnalytical(g):  # for G1 and G2
    """
    This function estimates two parameters of Gamma distribution 'a' and 'b' analytically, given data.

    Args:
    -----
        g (1D np.array): an array holding the data points

    Returns:
    --------
        a_hat_new (float): estimated shape parameter of Gamma distribution.
        b_hat_new (float): estimated shape parameter of Gamma distriubtion.
    """

    # calculate required mean and log_mean of data
    tau_mean = np.mean(g)
    tau_logmean = np.log(tau_mean)
    tau_meanlog = np.mean(np.log(g))

    # initialization step
    a_hat0 = 0.5 / (tau_logmean - tau_meanlog)  # shape
    b_hat0 = tau_mean / a_hat0  # scale
    psi_0 = np.log(a_hat0) - 1 / (2 * a_hat0)  # psi is the derivative of log of gamma function, which has been approximated as this term
    psi_prime0 = 1 / a_hat0 + 1 / (a_hat0 ** 2)  # this is the derivative of psi
    assert a_hat0 != 0, "the first parameter has been set to zero!"

    # updating the parameters
    for i in range(100):
        a_hat_new = (a_hat0 * (1 - a_hat0 * psi_prime0)) / (1 - a_hat0 * psi_prime0 + tau_meanlog - tau_logmean + np.log(a_hat0) - psi_0)
        b_hat_new = tau_mean / a_hat_new

        a_hat0 = a_hat_new
        psi_prime0 = 1 / a_hat0 + 1 / (a_hat0 ** 2)
        psi_0 = np.log(a_hat0) - 1 / (2 * a_hat0)
        psi_prime0 = 1 / a_hat0 + 1 / (a_hat0 ** 2)

        if np.abs(a_hat_new - a_hat0) <= 0.01:
            return [a_hat_new, b_hat_new]
        else:
            pass

    assert np.abs(a_hat_new - a_hat0) <= 0.01, "a_hat has not converged properly, a_hat_new - a_hat0 = {}".format(np.abs(a_hat_new - a_hat0))

    result = [a_hat_new, b_hat_new]
    return result
