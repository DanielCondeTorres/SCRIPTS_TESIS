import numpy as np 
import matplotlib.pyplot as plt
import sys

# parameters file
bacteria=sys.argv[1]
bacteria=open(bacteria)
cancer=sys.argv[2]
cancer=open(cancer)
popc=sys.argv[3]
popc=open(popc)

def RMSD(archivo):
    tiempo=[]
    RMSD=[]
    for line in archivo:
        if '#' in line or '@' in line:
            pass
        else:
            tiempo.append(float(line.split()[0])*10**-6)
            RMSD.append(float(line.split()[1]))
    return tiempo,RMSD
tiempo_bacteria,RMSD_bacteria=RMSD(bacteria)
tiempo_cancer,RMSD_cancer=RMSD(cancer)
tiempo_popc,RMSD_popc=RMSD(popc)
def media_movil(lista,orden=500):
    media_movil=[]
    i=0
    print(len(lista))
    while (i<orden-1):
            media_movil.append(np.mean(lista[0:i]))
            print(np.mean(lista[0:i]))
            i=i+1
    array=np.convolve(np.array(lista),np.ones(orden),'valid')/orden
    array=array.tolist()
    media_movil.extend(array)
    return media_movil
media_movil_bacteria=media_movil(RMSD_bacteria)
media_movil_cancer=media_movil(RMSD_cancer)
media_movil_popc=media_movil(RMSD_popc)
plt.plot(tiempo_bacteria,media_movil_bacteria,ls='dotted',linewidth=4,label='Bacteria')
plt.plot(tiempo_cancer,media_movil_cancer,ls=(0,(1,10)),label='Cancer',linewidth=6.5)
plt.plot(tiempo_popc,media_movil_popc,ls='dashed',label='Mammal',linewidth=4.5)
plt.legend(fontsize=28)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.xlabel(r'Time ($\mu s$)',fontsize=38)
plt.ylabel('RMSD (nm)',fontsize=38)
#plt.ylabel(r'Area ($nm^{2}$)',fontsize=38)
plt.xlim(0,12)
plt.show()
plt.close()

