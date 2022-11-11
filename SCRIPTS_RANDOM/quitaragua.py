#!/usr/bin/env python2.7

import optparse
import numpy as np

def main() :
    parser = optparse.OptionParser(version='%prog version 1.0')
    parser.add_option('-f', '--infile', help='input pdb file)', action='store')
    options, arguments = parser.parse_args()

#***********************************************************************************************

    infile =open(str(options.infile),"r")
    w=[]
    for line in infile:
        if line.strip():
            fields=line.split()
            if (len(fields)>1) and fields[0]=="ATOM":
                w.append(fields)
                print('fields',fields)
			
	#print 'w',w
    print('w',w)
    matrix=np.array(w)
    infile.close()
    infile =open(str(options.infile),"w")
    ww=[]; print('matriz', matrix)#;print('AGUAAAAAAAAAAAA:',w)
    fila,col=matrix.shape
    i=0
    print ('fila', fila)
    while i < fila:
        if matrix[i][2]== 'W' and float(matrix[i][7]) >73 and float(matrix[i][7]) < 124:
            i=i+3
        else:
            ww.append(matrix[i][:])
            print 
            i=i+1
            print ('i==',i)
            print ('filaa=',fila) 


    print ('wwwww',ww)

    k=0
    matri=np.array(ww)
    f,c=matri.shape
    for i in range(f):
        if matri[i][2]=='W':
            k=k+1
    print ('kkkk',k)
    for i in range(f):
        cfinal="{:6s}{:5s} {:^4s} {:3s}  {:4s}    {:8.6s}{:8.6s}{:8.7s}{:6.4s}{:6.4s} \n".format(matri[i][0],matri[i][1],matri[i][2],matri[i][3],matri[i][4],matri[i][5],matri[i][6],matri[i][7],matri[i][8],matri[i][9])
        infile.write(cfinal)	
    infile.close()
	

if __name__=="__main__" :
    main()

#if (fields[3]=="SOL" and  float(fields[7]) > 44.23 and float(fields[7]) < 96.79):
#					pass
	#			else:


