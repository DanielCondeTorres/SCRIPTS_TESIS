import numpy as np
from math import acos
import pandas as pd
print('created by Daniel Conde Torres: El Puto Amo')
HILLS=np.loadtxt('HILLS', comments='#')
tiempo=HILLS[:,0]
D=HILLS[:,1]
coseno=HILLS[:,2]
col_4=HILLS[:,3]
col_5=HILLS[:,4]
col_6=HILLS[:,5]
col_7=HILLS[:,6]
#columns = ['time','D.z','angle','sigma_D.z','sigma_angle','height','biasf']
#df=pd.read_fwf('HILLS2.dat',comment='#!',names=columns)
#for row in df.iterrows():
#    print('LINEA: ',row)
#print('COSENO: ',df['angle'])
#df['angle']=np.arccos(df['angle'])
D1=[];tiempo1=[];col_41=[];coseno1=[];col_51=[];col_61=[];col_71=[]
for i in range(len(D)):
    if D[i]<0:
        coseno1.append(coseno[i])
        D1.append(D[i]);tiempo1.append(tiempo[i]);col_41.append(col_4[i]);col_51.append(col_5[i]);col_61.append(col_6[i]);col_71.append(col_7[i])
    else:
        pass
HILL=open('HILLS','r')
with open('HILLS_angulo', 'w') as f:
    cont=0
    for linea in HILL:
        f.write(linea)
        cont=cont+1
        if cont==3:
            break
    for count in range(len(col_71)):
        f.write('{:>23.1f}{:>23.15f}{:>23.15f}{:>23.2f}{:>23.2f}{:>23.15f}{:>23.0f} '.format(tiempo1[count],D1[count],coseno1[count],col_41[count],col_51[count],col_61[count],col_71[count]))
        f.write("\n")
