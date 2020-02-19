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

def plotDist(control, trial, cont_pr, parameters, label):
    titles = ['G1', 'S-G2']

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
    
def pTotal(Control, Drug):
    control = polish(Control)
    drug = polish(Drug)

    xbarCg1 = np.sum(control[0]) / len(control[0])
    xbarCg2 = np.sum(control[1]) / len(control[1])
    xbarDg1 = np.sum(drug[0]) / len(drug[0])
    xbarDg2 = np.sum(drug[1]) / len(drug[1])

    sh1 = []
    sh2 = []
    for k in range(1,100):
        tmc1 = sp.gamma.logpdf(control[0], a=k, loc=0, scale=xbarCg1/k).sum()
        tmc2 = sp.gamma.logpdf(control[1], a=k, loc=0, scale=xbarCg2/k).sum()
        tmd1 = sp.gamma.logpdf(drug[0], a=k, loc=0, scale=xbarDg1/k).sum()
        tmd2 = sp.gamma.logpdf(drug[1], a=k, loc=0, scale=xbarDg2/k).sum()
        sh1.append(tmc1 + tmd1) # g1
        sh2.append(tmc2 + tmd2) # g2
        
    controlParams = [[np.argmax(sh1), xbarCg1/np.argmax(sh1)], [np.argmax(sh2), xbarCg2/np.argmax(sh2)]]
    drugParams = [[np.argmax(sh1), xbarDg1/np.argmax(sh1)], [np.argmax(sh2), xbarDg2/np.argmax(sh2)]]
    return control, drug, controlParams, drugParams