import numpy as np
import matplotlib.pyplot as plt
import optparse

def main() :
    parser = optparse.OptionParser(version='%prog version 1.0')
    parser.add_option('-f', '--infile', help='input pdb file)', action='store')
    options, arguments = parser.parse_args()
#****************************************************************************
    file_in = open(str(options.infile),'r')
    bead=[]
    pos=[]
    neg=[]
    neutro=[]
    apolar=[]
    aromatico=[]
    #next(file_in)
    for line in file_in:
        if 'CA' in line or 'BB' in line:
            bead.append(str(line.split()[3]))
            if str(line.split()[3])in ['GLY','MET','ALA','PRO','VAL','ILE','LEU']:
                apolar.append(str(line.split()[3]))
            elif str(line.split()[3]) in ['PHE','TRP']:
                aromatico.append(str(line.split()[3]))
            elif str(line.split()[3])in ['SER','THR','TYR','CYS','GLN','ASN']:
                neutro.append(str(line.split()[3]))
            elif str(line.split()[3]) in ['ASP','GLU']:
                neg.append(str(line.split()[3]))
            elif str(line.split()[3]) in ['LYS','ARG','HIS']:
                pos.append(str(line.split()[3]))
        if 'ENDMDL' in line:
            break
    file_in.close()
    print('BEAD:',len(bead))
    print('RESIDUOS',bead)
    print('RESIDUOS POLARES:',len(pos)+len(neg)+len(neutro))
    print('RESIDUOS APOLARES: ',len(apolar))
    print('RESIDUOS AROMATICOS: ',len(aromatico))
    print('NEGATIVOS: ',len(neg))
    print('POSITIVOS: ',len(pos))
    print('Polar Neutro: ',len(neutro))
    print('#######################################')
if __name__=="__main__" :
    main()
~                   
