import numpy as np
import matplotlib.pyplot as plt
Numero_beads=12
numero_beads=list(range(1,13))
print('p',numero_beads)
for i in numero_beads:
    file_in = open('aa_ref_1b3e_autopsf.pdb')
    bead=[]
    pos=[]
    neg=[]
    neutro=[]
    apolar=[]
    aromatico=[]
    ########################

    next(file_in)
    for line in file_in:
        if 'CA' in line:
            if float(line.split()[10])==i:
                bead.append(str(line.split()[3]))
                if str(line.split()[3])==('GLY' or 'MET' or 'ALA' or 'PRO' or 'VAL' or 'ILE' or 'LEU'):
                    apolar.append(str(line.split()[3]))
                elif str(line.split()[3])==('PHE' or 'TRP'):
                    aromatico.append(str(line.split()[3]))
                elif str(line.split()[3])==('SER' or 'THR' or 'TYR' or 'CYS' or 'GLN' or 'ASN'):
                    neutro.append(str(line.split()[3]))
                elif str(line.split()[3])==('ASP' or 'GLU'):
                    neg.append(str(line.split()[3]))
                elif str(line.split()[3])==('LYS' or 'ARG' or 'HIS'):
                    pos.append(str(line.split()[3]))
    file_in.close()
    print('BEAD:',i)
    print('RESIDUOS',bead)
    print('RESIDUOS POLARES:',len(pos)+len(neg)+len(neutro))
    print('RESIDUOS APOLARES: ',len(apolar))
    print('RESIDUOS AROMATICOS: ',len(aromatico))
    print('NEGATIVOS: ',len(neg))
    print('POSITIVOS: ',len(pos))
    print('#######################################')
