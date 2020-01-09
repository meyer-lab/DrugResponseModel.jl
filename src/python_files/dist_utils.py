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

def plotDist(TRIAL):
    trial = polish(TRIAL)

    plt.figure(figsize=(18,8))
    titles = ['G1', 'S-G2', 'total']
    for ind, x in enumerate(trial):
        plt.subplot(1,3,(ind+1))
        plt.hist(x, alpha=0.8, density=True, bins=30)
        plt.title(titles[ind])
        shape, loc, scale = sp.rv_continuous.fit(sp.gamma, x, floc=0)
        print(np.around(shape), round(scale, 3))
        rv = sp.erlang(np.around(shape),loc,scale)
        xx = np.linspace(0, max(x))
        plt.plot(xx, rv.pdf(xx), lw=2)
        plt.text(50, 0.07, 'shape %d,\n scale %.3f' % (np.around(shape), round(scale, 3)), fontsize=12)
        plt.ylim([0.0, 0.12])
        plt.xlim([0.0, 90.0])
