import pandas as pd
import numpy as np
import scipy.stats as sp
import matplotlib.pyplot as plt

def readData(path):
    df = pd.ExcelFile(path)

    CONTROL = pd.read_excel(df, 'CONTROL').values
    LAP25 = pd.read_excel(df, 'LAP25').values
    LAP50 = pd.read_excel(df, 'LAP50').values
    GEM5 = pd.read_excel(df, 'GEM5').values
    GEM10 = pd.read_excel(df, 'GEM10').values
    TAX2 = pd.read_excel(df, 'TAX2').values
    TAX3 = pd.read_excel(df, 'TAX3').values

    return CONTROL, LAP25, LAP50, GEM5, GEM10, TAX2, TAX3


def polish(CONTROL):
    trial = []
    for i in range(3):
        tmp1 = CONTROL[~np.isnan(CONTROL[:,i]), i]
        tmp2 = tmp1[tmp1 != 0]
        trial.append(tmp2)
    return trial

def estimate(trial):
    parameters = []
    for ind, x in enumerate(trial):
        shape, _, scale = sp.rv_continuous.fit(sp.gamma, x, floc=0)
        parameters.append([np.around(shape), round(scale, 3)])

    return parameters

def plotDist(control, trial, cont_pr, parameters, label):
    titles = ['G1', 'S-G2', 'total']

    plt.figure(figsize=(18,8), dpi=200)
    for ind, x in enumerate(trial):
        plt.subplot(1,3,(ind+1))
        plt.hist(x, alpha=0.8, density=True, bins=30, label=label)
        plt.hist(control[ind], alpha=0.8, density=True, bins=30, label="control")
        plt.title(titles[ind])
        rv = sp.erlang(parameters[ind][0],0,parameters[ind][1])
        xx = np.linspace(0, max(x))
        plt.plot(xx, rv.pdf(xx), lw=2)
        plt.text(50, 0.07, 'shape: %d \n scale: %.3f' % (parameters[ind][0],parameters[ind][1]), fontsize=14)
        plt.text(50, 0.09, 'control shape: %d \n scale: %.3f' % (cont_pr[ind][0],cont_pr[ind][1]), fontsize=14)
        plt.ylim([0.0, 0.12])
        plt.xlim([0.0, 90.0])
        plt.legend()
