import numpy as np
from os.path import exists
import matplotlib.pyplot as plt
fes=open('fes.dat')
xx=[]
yy=[]
zz=[]
nfesfiles=0
for i in range(1,1000): # Detect number of 1D fes files
    if exists("fes_"+str(i)+".dat"): nfesfiles=nfesfiles+1
    else: break
n1Dplots=int(nfesfiles*10/100)


for line in fes:
	if '#' in line:#para saltar as liñas comentadas
		pass
	else:
		xx.append(float(line.split()[0]))
		yy.append(float(line.split()[1]))
		zz.append(float(line.split()[2]))


for indices in range(nfesfiles,nfesfiles-n1Dplots,-1):
	file_in = open('fes_'+str(indices)+'.dat')
	x=[]
	y=[]
	z=[]
	for line in file_in:
		if '#' in line:#para saltar as liñas comentadas
			pass
		else:
			x.append(float(line.split()[0]))
			y.append(float(line.split()[1]))
			z.append(float(line.split()[2])) 
	plt.xlim(0,30);plt.plot(x,y,linestyle='--')
	plt.ylim(0,150)
	plt.xlabel('Position (nm)',fontsize=38)
plt.legend(fontsize=28)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.ylabel('Energy (kJ/mol)',fontsize=38)
plt.legend(fontsize=28)






print(y)


##################################
#-----INCERTIDUMBRE----------#####
##################################
a=np.loadtxt('CV1_BA.dat')
energia=a[:,1]
pos=a[:,0]
print(len(pos))
print(len(y))
inc=a[:,2]
plt.plot(pos,energia,color='green')
plt.fill_between(pos, np.array(energia)-np.array(inc),np.array(energia)+np.array(inc),alpha=0.3,color='green')







plt.show()
