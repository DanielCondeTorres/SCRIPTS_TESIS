import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib.ticker import MaxNLocator, IndexFormatter


# parameters file
bacteria=sys.argv[1]
bacteria=open(bacteria)
cancer=sys.argv[2]
cancer=open(cancer)
popc=sys.argv[3]
popc=open(popc)
def RMSF(archivo):
    ATOM=[]
    RMSF=[]
    for line in archivo:
        if '#' in line or '@' in line:
            pass
        else:
            ATOM.append(float(line.split()[0]))
            RMSF.append(float(line.split()[1]))
    return ATOM,RMSF
ATOM_bacteria,RMSF_bacteria=RMSF(bacteria)
ATOM_cancer,RMSF_cancer=RMSF(cancer)
ATOM_popc,RMSF_popc=RMSF(popc)
print('ATOM',ATOM_popc)
#DICCIONARIO
d1 = {
  "GLY": "G",
  "ALA": "A",
  "ARG": "R",
  "ASN": "N",
  "ASP": "D",
  "CYS": "C",
  "GLU": "E",
  "GLN": "Q",
  "HIS": "H",
  "ILE": "I",
  "SEE": "U",
  "LEU": "L",
  "LYS": "K",
  "MET": "M",
  "PHE": "F",
  "PRO": "P",
  "SER": "S",
  "THR": "T",
  "TRP": "W",
  "TYR": "Y",
  "VAL": "V",
  "PYL": "O",
}

EJE_X=[]
LISTA_ATOMS=[]
pdb=sys.argv[4]
pdb=open(pdb)
for line in pdb:
	if not 'ATOM' in line:
		pass
	else:
		LISTA_ATOMS.append(str(line.split()[3]))
	
count=0
print(len(LISTA_ATOMS))
EJE_X.append('a')
for elemento in LISTA_ATOMS:
        count=count+1
        try:
            if count==30:
                        print('ha')
            if count==1 or count==2 or count==4 or count==5 or count==8 or count==12 or count==14 or count==18 or count==20 or count==21 or count==24 or count==27 or count==31 or count==32 or count==40 or count==45 or count==42 or count==43 or count==35 or count==36 or count==47 or count==49 or count==51:
                        print('count',count)
                        print('ELEMENTO',d1[elemento])
                        EJE_X.append(d1[elemento])
        except:
                continue
import matplotlib.ticker as ticker
print(EJE_X)
fig, ax = plt.subplots()
#plt.xlim(1,55)
ax.set_xticklabels(EJE_X)
ax.xaxis.set_major_locator(ticker.MultipleLocator(2.4))
#ax.xaxis.set_major_locator(MaxNLocator(26))
print(len(RMSF_bacteria))
print(len(EJE_X))
#############
#plt.scatter(ATOM_cancer,RMSF_bacteria)
plt.plot(ATOM_cancer,RMSF_bacteria,ls='dotted',linewidth=4,label='Bacteria')
plt.plot(ATOM_cancer,RMSF_cancer,ls=(0,(4,10)),marker='*',linewidth=6.5,label='Cancer')
plt.plot(ATOM_cancer,RMSF_popc,ls='dashed',label='Mammal',linewidth=4.5)
plt.legend(fontsize=28)
plt.xticks(fontsize=25)
plt.yticks(fontsize=30)
plt.xlabel(r'ATOM',fontsize=38)
plt.ylabel('RMSF (nm)',fontsize=38)
plt.show()
#plt.close()

