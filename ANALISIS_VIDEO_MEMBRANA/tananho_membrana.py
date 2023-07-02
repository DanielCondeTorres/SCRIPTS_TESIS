#!/usr/bin/env python2.7

import optparse
from numpy import *
def main() :
	parser = optparse.OptionParser(version='%prog version 1.0')
	parser.add_option('-f', '--infile', help='input pdb file)', action='store')
	options, arguments = parser.parse_args()

#***********************************************************************************************
	#from numpy import *
	infile =open(str(options.infile),"r")
	z=[]
	x=[]
	y=[]
	for line in infile:
		if line.strip() :
			fields=line.split()
			if (len(fields)>1) and fields[0]=="ATOM":
				#if (fields[3]=="POPC"):
                            if (fields [3]=="CHOL" or fields[3]=="DOPC" or fields[3]=="DOPE" or fields[3]=="DOPS" or fields[3]=="DPSM" or fields[3]=="POPC" or fields[3]=="POPS"):
                                z.append(float(fields[7]))
                                x.append(float(fields[5]))
                                y.append(float(fields[6]))
	infile.close()
	print( 'zmax', max(z))
	print( 'zmin', min(z))
	print( 'xmax', max(x))
	print( 'xmin', min(x))
	print( 'ymax', max(y))
	print( 'ymin', min(y))
if __name__=="__main__" :
    main()
