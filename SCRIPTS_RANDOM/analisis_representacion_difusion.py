import numpy as np
import matplotlib.pyplot as plt
import sys


# parameters file
#bacteria=sys.argv[1]
bacteria=open('bacteria_20ns.ldc')
#cancer=sys.argv[2]
cancer=open('cancer_20ns.ldc')
#popc=sys.argv[3]
popc=open('popc_20ns.ldc')
def RMSF(archivo):
    pos=[]
    P=[]
    for line in archivo:
        if '#' in line or '@' in line:
            pass
        else:
            pos.append(float(line.split()[0]))
            P.append(float(line.split()[2]))
    return pos,P
pos_bacteria,P_bacteria=RMSF(bacteria)
pos_cancer,P_cancer=RMSF(cancer)
pos_popc,P_popc=RMSF(popc)

plt.plot(pos_bacteria,P_bacteria,linewidth=4.5,label='Bacteria')
plt.plot(pos_cancer,P_cancer,linewidth=4.5,label='Cancer')
plt.plot(pos_cancer,P_popc,label='Mammal',linewidth=4.5)
plt.legend(fontsize=28)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.xlabel('Distance (nm)',fontsize=38)
plt.ylabel('P',fontsize=38)
plt.xlim(0,60)
plt.ylim(0,0.16)
plt.show()
#plt.close()
