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

    return CONTROL[0:-1,0:2], LAP25[0:-1,0:2], LAP50[0:-1,0:2], GEM5[0:-1,0:2], GEM10[0:-1,0:2], TAX2[0:-1,0:2], TAX3[0:-1,0:2]

def polish(CONTROL):
    # remove NaNs and zeros
    trial = []
    for i in range(2):
        tmp1 = CONTROL[~np.isnan(CONTROL[:,i]), i]
        tmp2 = tmp1[tmp1 != 0]
        trial.append(tmp2)
    return trial

def estimate(CONTROL):
    trial = polish(CONTROL)
    # estimate
    parameters = []
    for ind, x in enumerate(trial):
        shape, _, scale = sp.rv_continuous.fit(sp.gamma, x, floc=0)
        parameters.append([np.around(shape), round(scale, 3)])

    return parameters

def plotDist(Control, drug, cont_pr, parameters, label):
    titles = ['G1', 'S-G2']

    trial = polish(drug)
    control = polish(Control)
    plt.figure(figsize=(12,6), dpi=200)
    for ind, x in enumerate(trial):
        plt.subplot(1,2,(ind+1))
        plt.hist(control[ind], alpha=0.4, density=True, label="control", color="slategrey")
        plt.hist(x, alpha=0.4, density=True, label=label, color="orchid")
        plt.title(titles[ind])
        rv = sp.erlang(parameters[ind][0],0,parameters[ind][1])
        rv2 = sp.erlang(cont_pr[ind][0],0,cont_pr[ind][1])
        xx = np.linspace(0, max(x))
        plt.plot(xx, rv2.pdf(xx), lw=2, color="slategrey")
        plt.plot(xx, rv.pdf(xx), lw=2, color="orchid")
        plt.text(50, 0.07, 'shape: %d \n scale: %.3f' % (parameters[ind][0],parameters[ind][1]), fontsize=14)
        plt.text(50, 0.09, 'control shape: %d \n scale: %.3f' % (cont_pr[ind][0],cont_pr[ind][1]), fontsize=14)
        plt.ylim([0.0, 0.12])
        plt.xlim([0.0, 90.0])
        plt.xlabel("phase durations [hrs]")
        plt.ylabel("probability")
        plt.legend()
    
def pTotal(contP, Drug):
    shapeG1 = contP[0][0]
    shapeG2 = contP[1][0]
    scale_contG1 = contP[1][0]
    scale_contG2 = contP[1][1]
    drug = polish(Drug)

    scale_drugRange1 = np.linspace(0.1*scale_contG1, 10*scale_contG1, 100)
    scale_drugRange2 = np.linspace(0.1*scale_contG2, 10*scale_contG2, 100)
    cnt1 = []
    cnt2 = []
    for ind, sc in enumerate(scale_drugRange1):
        cnt1.append(sp.gamma.logpdf(drug[0], a=shapeG1, loc=0, scale=sc).sum())
        cnt2.append(sp.gamma.logpdf(drug[1], a=shapeG2, loc=0, scale=scale_drugRange2[ind]).sum())
    drugParams = [[round(shapeG1), scale_drugRange1[np.argmax(cnt1)]], [round(shapeG2), scale_drugRange2[np.argmax(cnt2)]]]
    return drugParams