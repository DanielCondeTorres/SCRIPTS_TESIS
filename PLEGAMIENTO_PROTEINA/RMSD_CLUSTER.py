import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from MDAnalysis.analysis import contacts, rms
from matplotlib.cm import get_cmap
from MDAnalysis.analysis import dihedrals
from MDAnalysis.lib.distances import calc_angles
# Cargar el archivo xtc y el archivo pdb
from MDAnalysis.analysis.rms import rmsd as mdarmsd
import warnings

from MDAnalysis.tests.datafiles import waterPSF, waterDCD
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis




from MDAnalysis.analysis.dihedrals import Ramachandran
import matplotlib.pylab as plt
import MDAnalysis as mda
import MDAnalysis.analysis.hbonds
import MDAnalysis.analysis.rms




def explicit_heat_smooth(prices: np.array,
                         t_end: float = 3.0) -> np.array:
    '''
    Smoothen out a time series using a simple explicit finite difference method.
    The scheme uses a first-order method in time, and a second-order centred
    difference approximation in space. The scheme is only numerically stable
    if the time-step 0<=k<=1.
    
    The prices are fixed at the end-points, so the interior is smoothed.
    Parameters
    ----------
    prices : np.array
        The price to smoothen
    t_end : float
        The time at which to terminate the smootheing (i.e. t = 2)
        
    Returns
    -------
    P : np.array
        The smoothened time-series
    '''
    
    k = 0.1 # Time spacing
    
    # Set up the initial condition
    P = prices
    
    t = 0
    while t < t_end:
        # Solve the finite difference scheme for the next time-step
        P = k*(P[2:] + P[:-2]) + P[1:-1]*(1-2*k)
        
        # Add the fixed boundary conditions since the above solves the interior
        # points only
        P = np.hstack((
            np.array([prices[0]]),
            P,
            np.array([prices[-1]]),
        ))
        t += k

    return P



def RAMACHANDRAM_PLOT(archivo_gro,archivo_xtc):
    u=mda.Universe(archivo_gro,archivo_xtc)
    proteina=u.select_atoms('protein')
    r=Ramachandran(u.select_atoms('protein')).run()
    #print('psi:',r.results.angles[0][:, 0])
    return r
'''fig, ax = plt.subplots(figsize=plt.figaspect(1))
r.plot(ax=ax, color='k', marker='s',ref=True)

plt.ylabel(r'$\psi$',fontsize=30)
plt.xlabel(fr'$\phi$ ',fontsize=30)
plt.xticks(fontsize=28)
plt.yticks(fontsize=28)
plt.show()
plt.close()'''







warnings.filterwarnings('ignore')

def analisis_cluster(archivo_gro,archivo_xtc,referencia):
    u = mda.Universe(archivo_gro, archivo_xtc)
    ref=mda.Universe(referencia)
    # Cargar el archivo xtc y pdb con MDAnalysis
    # Seleccionar los átomos de interés (por ejemplo, átomos de carbono alfa)
    seleccion = u.select_atoms('backbone')
    seleccion_ref=ref.select_atoms('backbone')
    # Crear una matriz de características basada en el porcentaje de hélice alfa y el radio de giro
    caracteristicas = []
    tiempos = []
    almacenar_helice=[]
    i=0
    r=RAMACHANDRAM_PLOT(archivo_gro,archivo_xtc)
    for ts in u.trajectory:
        # Calcular los ángulos dihedrales para determinar la presencia de hélices alfa
        phi=r.results.angles[i][:, 0]
        psi=r.results.angles[i][:, 1]
        #print('a',phi[i])
        i=i+1
        # Calcular el porcentaje de hélice alfa
        helices = np.logical_and(phi > -150, phi < -30) & np.logical_and(psi > -120, psi < 30)  # Rangos de ángulos para considerar hélices alfa
        porcentaje_helix = np.sum(helices) / float(len(helices)) * 100
        # Calcular el radio de giro
        radio_giro = seleccion.radius_of_gyration()
        # Agregar las características a la matriz
        caracteristicas.append([porcentaje_helix, radio_giro])
        almacenar_helice.append(porcentaje_helix)
        tiempos.append(ts.time)
    helix_structure=explicit_heat_smooth(np.array(almacenar_helice),len(almacenar_helice))
    caracteristicas = np.array(caracteristicas)
    tiempos = np.array(tiempos)
    # Realizar el análisis de clustering por k-means con 6 clusters
    n_clusters = 6
    kmeans = KMeans(n_clusters=n_clusters)
    kmeans.fit(caracteristicas)
    	# Obtener las etiquetas de clustering asignadas a cada estructura
    etiquetas_cluster = kmeans.labels_
    # Calcular el RMSD de cada estructura respecto a la primera estructura
    rmsd_valores = []
    for ts in u.trajectory:
    #    rmsd_valores.append(mdarmsd(seleccion.positions, seleccion.positions, superposition=True))
        rmsd_valores.append(rms.rmsd(seleccion.positions,seleccion_ref.positions , superposition=True))
    rmsd_valores = np.array(rmsd_valores)
    return tiempos, rmsd_valores,n_clusters,etiquetas_cluster,helix_structure,r





carpetas=['GROMOS_QUAD','GROMOS_DEU','GROMOS_HMR','GROMOS_D','AMBER_D','CHARMM_D']
lista_de_simulaciones=['Quadium','Deuterium','HMR','GROMOS','AMBER','CHARMM']#lista de simulaciones a representar
fig, axs = plt.subplots(3, 2,sharex=True, sharey=True)#son 6 simulaciones, por lo que 3*2
fig.subplots_adjust(wspace=0, hspace=0)

figs, axsw = plt.subplots(3, 2,sharex=True, sharey=True)#son 6 simulaciones, por lo que 3*2
figs.subplots_adjust(wspace=0, hspace=0)







figrama, axsrama = plt.subplots(3, 2,sharex=True, sharey=True)#son 6 simulaciones, por lo que 3*2
#figrama.subplots_adjust(wspace=0, hspace=0)





for ff in range(len(lista_de_simulaciones)):
    tiempos,rmsd_valores,n_clusters,etiquetas_cluster,helice,ramachandran=analisis_cluster(str(carpetas[ff])+'/conc.pdb',str(carpetas[ff])+'/conc.xtc',str(carpetas[ff])+'/ref.pdb')
    colores = ['r', 'g', 'b', 'c', 'm', 'y']  # Colores para cada cluster
    # Plotear los puntos del RMSD con colores según el cluster
    tiempos=tiempos/1000
    if ff == 0 or  ff == 1:
        axs[0,ff].plot(tiempos,rmsd_valores,color='k')
        axsw[0,ff].plot(tiempos,helice,color='k')
        ramachandran.plot(ax=axsrama[0,ff] , color='k', marker='s',ref=True)

        titulo_1=axs[0,ff];titulo_2=axsw[0,ff]
        for i in range(n_clusters):
            indices_cluster = np.where(etiquetas_cluster == i)[0]
            porcentajes=len(indices_cluster)/len(rmsd_valores)*100
            percentaje=np.round(porcentajes,2)
            axs[0,ff].scatter(tiempos[indices_cluster], rmsd_valores[indices_cluster], c=colores[i], label=f'Cluster {i+1}'+' '+str(percentaje)+'%')
            axs[0,ff].legend(fontsize=10,ncol=3,markerscale=2.,loc='lower right');axs[0,ff].tick_params(labelsize=18)

    elif ff == 2 or ff == 3:
        axs[1,ff-2].plot(tiempos,rmsd_valores,color='k')
        axsw[1,ff-2].plot(tiempos,helice,color='k')
        ramachandran.plot(ax=axsrama[1,ff-2] , color='k', marker='s',ref=True)


        titulo_1=axs[1,ff-2];titulo_2=axsw[1,ff-2]
        for i in range(n_clusters):
            indices_cluster = np.where(etiquetas_cluster == i)[0]
            porcentajes=len(indices_cluster)/len(rmsd_valores)*100
            percentaje=np.round(porcentajes,2)
            axs[1,ff-2].scatter(tiempos[indices_cluster], rmsd_valores[indices_cluster], c=colores[i], label=f'Cluster {i+1}'+' '+str(percentaje)+'%')
            axs[1,ff-2].legend(fontsize=10,ncol=3,markerscale=2.,loc='lower right');axs[1,ff-2].tick_params(labelsize=18)
    else:
        axs[2,ff-4].plot(tiempos,rmsd_valores,color='k')
        axsw[2,ff-4].plot(tiempos,helice,color='k')
        ramachandran.plot(ax=axsrama[2,ff-4] , color='k', marker='s',ref=True)

        titulo_1=axs[2,ff-4];titulo_2=axsw[2,ff-4]
        for i in range(n_clusters):
            indices_cluster = np.where(etiquetas_cluster == i)[0]
            porcentajes=len(indices_cluster)/len(rmsd_valores)*100
            percentaje=np.round(porcentajes,2)
            axs[2,ff-4].scatter(tiempos[indices_cluster], rmsd_valores[indices_cluster], c=colores[i], label=f'Cluster {i+1}'+' '+str(percentaje)+'%')
            axs[2,ff-4].legend(fontsize=10,ncol=3,markerscale=2.,loc='lower right');axs[2,ff-4].tick_params(labelsize=18)
    titulo_2.tick_params(labelsize=20)
    titulo_1=titulo_1.set_title(str(lista_de_simulaciones[ff]),y=0.6,x=0.5);titulo_1.set_position([.08,0.85]);titulo_1.set_bbox(dict(facecolor='white', edgecolor='black', pad=5.0))
    titulo_2=titulo_2.set_title(str(lista_de_simulaciones[ff]),y=0.6,x=0.5);titulo_2.set_position([.08,0.85]);titulo_2.set_bbox(dict(facecolor='white', edgecolor='black', pad=5.0))

fig.text(0.5, 0.04, r'Time (ns)', ha='center', fontsize=22)
fig.text(0.07, 0.5, r'RMSD (${\AA}$)', va='center', rotation='vertical', fontsize=22)


figs.text(0.5, 0.04, r'Time (ns)', ha='center', fontsize=22)
figs.text(0.07, 0.5, r'%$\mathbf{(\alpha)-helix}$', va='center', rotation='vertical', fontsize=22)
plt.minorticks_on()
plt.show()



figs_hidrogen, axs_hidrogeno = plt.subplots(3, 2,sharex=True, sharey=True)#son 6 simulaciones, por lo que 3*2
figs_hidrogen.subplots_adjust(wspace=0, hspace=0)


import MDAnalysis
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import (
  HydrogenBondAnalysis as HBA)

def calcular_enlaces_hidrogeno(pdb_file,xtc_file):
    u = mda.Universe(pdb_file, xtc_file)
    hbonds = HBA(universe=u,donors_sel='name O or name N',hydrogens_sel='name H or name HN',acceptors_sel='name O or name N',d_a_cutoff=3.5)
    hbonds.run()
    hbonds_smooth=explicit_heat_smooth(hbonds.count_by_time(),len(hbonds.count_by_time()))

    return hbonds.times,hbonds_smooth


# Calcular el número de enlaces de hidrógeno a lo largo de la trayectoria
for ff in range(len(lista_de_simulaciones)):
    tempo,enlaces_H=calcular_enlaces_hidrogeno(str(carpetas[ff])+'/conc.pdb',str(carpetas[ff])+'/conc.xtc')
    if ff == 0 or  ff == 1:
        axs_hidrogeno[0,ff].plot(tempo/1000,enlaces_H,color='k')
        titulo_1=axs_hidrogeno[0,ff]
    elif ff == 2 or ff == 3:
        axs_hidrogeno[1,ff-2].plot(tempo/1000,enlaces_H,color='k')
        titulo_1=axs_hidrogeno[1,ff-2]
    else:
        axs_hidrogeno[2,ff-4].plot(tempo/1000,enlaces_H,color='k')
        titulo_1=axs_hidrogeno[2,ff-4]
    titulo_1.tick_params(labelsize=20)
    titulo_1=titulo_1.set_title(str(lista_de_simulaciones[ff]),y=0.6,x=0.5);titulo_1.set_position([.08,0.85]);titulo_1.set_bbox(dict(facecolor='white', edgecolor='black', pad=5.0))
figs_hidrogen.text(0.5, 0.04, r'Time (ns)', ha='center', fontsize=22)
figs_hidrogen.text(0.07, 0.5, r'# HYDROGEN BONDS', va='center', rotation='vertical', fontsize=22)


plt.show()
plt.close()


from matplotlib.ticker import MaxNLocator
####################RMSF
fig_rmsf, ax_rmsf = plt.subplots(3, 2,sharex=True, sharey=True)#son 6 simulaciones, por lo que 3*2
fig_rmsf.subplots_adjust(wspace=0, hspace=0)
ax_rmsf=ax_rmsf.ravel()
count=0
EJE_X=['0','ACE','LEU','ARG','ALA','LEU','ALA','ALA','LEU','ALA','ARG','ALA','ALA','ALA','ALA','ALA','ALA','ARG']
for elemento in range(len(lista_de_simulaciones)):
    u=mda.Universe(str(carpetas[elemento])+'/conc.pdb',str(carpetas[elemento])+'/conc.xtc',in_memory=True)
#########################################################
    protein=u.select_atoms('protein')#seleccionamos nuestro sistema
    protein=u.select_atoms('protein and name CA')
    R=MDAnalysis.analysis.rms.RMSF(protein)
    A=R.run()
    t=np.arange(1,len(A.rmsf)+1,1)
    t=np.arange(1,len(A.rmsf)+1,1)

# Creating plot
    ax_rmsf[elemento].plot(t,A.rmsf,marker='o',color='k')

    ax_rmsf[elemento].set_xticklabels(EJE_X)
    ax_rmsf[elemento].xaxis.set_major_locator(MaxNLocator(18))
    count+=1
    titulo_1=ax_rmsf[elemento]
    titulo_1.tick_params(labelsize=20)
    titulo_1=titulo_1.set_title(str(lista_de_simulaciones[elemento]),y=0.6,x=0.5);titulo_1.set_position([.08,0.85]);titulo_1.set_bbox(dict(facecolor='white', edgecolor='black', pad=5.0))
#ax.set_xlabel('X-axis') 
#ax.set_ylabel('Y-axis') 
 
fig_rmsf.text(0.5, 0.04, r'RMSF (${\AA}$)', ha='center', fontsize=22)
fig_rmsf.text(0.07, 0.5, r'C$\alpha$', va='center', rotation='vertical', fontsize=22)

plt.show()
plt.close()





