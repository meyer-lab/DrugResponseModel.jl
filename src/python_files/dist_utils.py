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

    NCg1 = len(control[0])
    NCg2 = len(control[1])
    NDg1 = len(drug[0])
    NDg2 = len(drug[1])
    # estimate scale for control in G1 and G2
    x_lnxCg1 = [x * np.log(x) for x in control[0]]
    lnxCg1 = [np.log(x) for x in control[0]]
    x_lnxCg2 = [x * np.log(x) for x in control[1]]
    lnxCg2 = [np.log(x) for x in control[1]]
    scaleCg1 = ((1 + 1e-10) / (NCg1 ** 2 + 1e-10)) * (NCg1 * (sum(x_lnxCg1)) - (sum(lnxCg1)) * (sum(control[0])))
    scaleCg2 = ((1 + 1e-10) / (NCg2 ** 2 + 1e-10)) * (NCg2 * (sum(x_lnxCg2)) - (sum(lnxCg2)) * (sum(control[1])))

    # estimate scale for drug in G1 and G2
    x_lnxDg1 = [x * np.log(x) for x in drug[0]]
    lnxDg1 = [np.log(x) for x in drug[0]]
    x_lnxDg2 = [x * np.log(x) for x in drug[1]]
    lnxDg2 = [np.log(x) for x in drug[1]]
    scaleDg1 = ((1 + 1e-10) / (NDg1 ** 2 + 1e-10)) * (NDg1 * (sum(x_lnxDg1)) - (sum(lnxDg1)) * (sum(drug[0])))
    scaleDg2 = ((1 + 1e-10) / (NDg2 ** 2 + 1e-10)) * (NDg2 * (sum(x_lnxDg2)) - (sum(lnxDg2)) * (sum(drug[1])))

    sh1 = []
    sh2 = []
    for k in range(1,100):
        tmc1 = sp.gamma.logpdf(control[0], a=k, loc=0, scale=scaleCg1).sum()
        tmc2 = sp.gamma.logpdf(control[1], a=k, loc=0, scale=scaleCg2).sum()
        tmd1 = sp.gamma.logpdf(drug[0], a=k, loc=0, scale=scaleDg1).sum()
        tmd2 = sp.gamma.logpdf(drug[1], a=k, loc=0, scale=scaleDg2).sum()
        sh1.append(tmc1 + tmd1) # g1
        sh2.append(tmc2 + tmd2) #g2
        
    controlParams = [[np.argmax(sh1), scaleCg1], [np.argmax(sh2), scaleCg2]]
    drugParams = [[np.argmax(sh1), scaleDg1], [np.argmax(sh2), scaleDg2]]
    return control, drug, controlParams, drugParams