import numpy as np
from math import acos
import pandas as pd
HILLS=np.loadtxt('HILLS', comments='#')
tiempo=HILLS[:,0]
D=HILLS[:,1]
coseno=HILLS[:,2]
col_4=HILLS[:,3]
col_5=HILLS[:,4]
col_6=HILLS[:,5]
col_7=HILLS[:,6]
transformacion=np.arccos(coseno)
#columns = ['time','D.z','angle','sigma_D.z','sigma_angle','height','biasf']
#df=pd.read_fwf('HILLS2.dat',comment='#!',names=columns)
#for row in df.iterrows():
#    print('LINEA: ',row)
#print('COSENO: ',df['angle'])
#df['angle']=np.arccos(df['angle'])
HILL=open('HILLS','r')
with open('HILLS_DANI', 'w') as f:
    cont=0
    for linea in HILL:
        f.write(linea)
        cont=cont+1
        if cont==3:
            break
    for count in range(len(col_7)):
        f.write('{:>23.1f}{:>23.15f}{:>23.15f}{:>23.2f}{:>23.2f}{:>23.15f}{:>23.0f} '.format(tiempo[count],D[count],transformacion[count],col_4[count],col_5[count],col_6[count],col_7[count]))
        f.write("\n")
