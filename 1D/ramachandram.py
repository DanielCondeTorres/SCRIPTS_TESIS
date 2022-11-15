# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 02:00:48 2020

@author: danie
"""



from MDAnalysis.analysis.dihedrals import Ramachandran
import matplotlib.pylab as plt
import MDAnalysis as mda
#from MDAnalysis.analysis import arg
print('aaaa')
import MDAnalysis.analysis.hbonds
import MDAnalysis.analysis.rms
u=mda.Universe('descargar.gro','descargar.xtc')
print('u',u)
proteina=u.select_atoms('protein')
print(proteina)
#=u.select_atoms("resid 0")
#print('pasppsapsa',r)
#R = Ramachandran(u).run()
r=Ramachandran(u.select_atoms('protein')).run()
print('vai',r)
fig, ax = plt.subplots(figsize=plt.figaspect(1))
r.plot(ax=ax, color='k', marker='s',ref=True)

plt.ylabel(r'$\psi$',fontsize=30)
plt.xlabel(fr'$\phi$ ',fontsize=30)
plt.xticks(fontsize=28)
plt.yticks(fontsize=28)
plt.show()
plt.close()

