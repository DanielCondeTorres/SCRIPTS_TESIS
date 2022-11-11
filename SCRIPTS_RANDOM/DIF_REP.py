# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 07:15:37 2022

@author: DANIEL
"""

import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt


# Set the figure size
un_ns = [11.83,7.97,10.016]#POPC, CANCER, BACTERIA 1 ns
dos_ns=[5.43,3.98,5],#POPC, CANCER, BACTERIA 2 ns
cinco_ns=[4.15,3.1,4.08]#POPC, CANCER, BACTERIA 5 ns
diez_ns=[3.84,0.919494,1.142073]#POPC, CANCER, BACTERIA 10 ns
veinte_ns=[3.75,2.9,3.42]


un_ns_std=[1.28,0.46,0.51]
dos_ns_std=[0.34,0.15,0.25]
cinco_ns_std=[0.25,0.18,0.25]
diez_ns_std=[0.017924,0.048838,0.054638]
veinte_ns_std=[0.22,0.091,0.17]
X=np.arange(3)
fig = plt.figure()


plt.bar(1, cinco_ns[0], color = "green", yerr=cinco_ns_std[0]*1.96,width = 0.05,label='POPC')
plt.bar(1.05,cinco_ns[1], color = "orange", yerr=cinco_ns_std[1]*1.96,width = 0.05)
plt.bar(1.1, cinco_ns[2], color = "blue", yerr=cinco_ns_std[2]*1.96,width = 0.05,label='BACTERIA')
barWidth = 0.05

#plt.xticks([r + 0.224 for r in range(len(X))], ['2ns',  '5ns','10ns'],fontsize=18)
plt.yticks(fontsize=35)
plt.xticks(fontsize=0)
plt.xlim(0.97,1.3)
plt.ylim(0,15)
plt.xlabel("Window",fontsize=1)
plt.ylabel(r'D $(10^{-7} cm^{2}/s)$' ,fontsize=38)
plt.show()
plt.close()
