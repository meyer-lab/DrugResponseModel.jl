import matplotlib.pyplot as plt
import numpy as np
N = 5

atOnce = (125.14822260315185, 126.85670180327187, 131.99718777067014, 116.21571145670654, 113.42232819941934);
null = (129.82445945984998, 161.8148149055424, 161.54913149878365, 105.1436332568337, 125.1031959713262);
sep = (101.79621387717243, 108.54398578263337, 107.21235662578536, 92.3456289085338, 99.6779558796112)
ind = np.arange(N) 
width = 0.25    
plt.figure(dpi=200)
plt.bar(ind, atOnce, width, label='at once')
plt.bar(ind+(2*width), sep, width, label='separate')
plt.bar(ind + width, null, width,
    label='null')
plt.ylabel('SSE')
plt.title('SSE for all concentrations of a drug')
plt.xticks(ind + width / 2, ('lapatinib', 'doxorubicin', 'gemcitabine', 'paclitaxel', 'palbociclib'))
plt.legend(loc='best')
plt.show()