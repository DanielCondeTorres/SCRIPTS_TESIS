# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 18:44:45 2020

@author: Daniel
"""
#import netCDF4
from scipy import *
import numpy as np
import matplotlib.pylab as plt
import MDAnalysis as mda
from MDAnalysis.analysis import align
import MDAnalysis.analysis.rms
from scipy.io import netcdf_file





import scipy.io.netcdf
import numpy as np

import MDAnalysis

import netCDF4

#######################################################################
#########################RMSF USANDO DE REFERENCIA PDBBANK#############
#######################################################################
archivos=['GROMOS_D','GROMOS_HMR','AMBER_D','GROMOS_DEU','CHARMM_D','GROMOS_QUAD']
labels=['GROMOS','HMR','AMBER','Deuterium','CHARMM','Quadium']
colors=['blue','k','purple','green','red','orange']
fig, ax = plt.subplots(figsize =(10, 7))
count=0
for elemento in archivos:
    print(elemento)
    u=mda.Universe('../../'+str(elemento)+'/conc.pdb','../../'+str(elemento)+'/conc.xtc',in_memory=True)
#########################################################
    protein=u.select_atoms('protein')#seleccionamos nuestro sistema

    average = align.AverageStructure(u, u, select='protein and (name CA or name CH3)',
                                 ref_frame=0).run()
    ref = average.results.universe
    protein=u.select_atoms('protein and (name CA or name CH3)')
    aligner = align.AlignTraj(u, ref,
                          select='protein and (name CA or name CH3)',
                          in_memory=True).run()

    R=MDAnalysis.analysis.rms.RMSF(protein)
    A=R.run()
    t=np.arange(1,len(A.rmsf)+1,1)
    t=np.arange(1,len(A.rmsf)+1,1)
    from matplotlib.ticker import MaxNLocator
# Creating plot
    ax.plot(t,A.rmsf,marker='o',color=colors[count],label=labels[count])
    count+=1
#ax.set_xlabel('X-axis') 
#ax.set_ylabel('Y-axis') 
  
ax.xaxis.set_major_locator(MaxNLocator(18))
EJE_X=['0','ACE','LEU','ARG','ALA','LEU','ALA','ALA','LEU','ALA','ARG','ALA','ALA','ALA','ALA','ALA','ALA','ARG']
ax.set_xticklabels(EJE_X)
ax.legend(fontsize=22,ncol=3)
ax.tick_params(axis ='y',labelsize=20)
ax.tick_params(axis ='x', labelsize=20)
ax.tick_params(axis ='both', labelsize=20)
plt.ylabel(r'RMSF (${\AA}$)',fontsize=30)
plt.xlabel(r'C$\alpha$',fontsize=30)
plt.show()
plt.close()

ref=mda.Universe('out.pdb')
R = MDAnalysis.analysis.rms.RMSD(u, ref,select="backbone")                              # NMP
R.run()

