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
    return tiempos,helix_structure




globales=['AGUA','TFE','TFE_HELICE_HECHA']
carpetas=['GROMOS_QUAD','GROMOS_DEU','GROMOS_HMR','GROMOS_D','AMBER_D','CHARMM_D']
lista_de_simulaciones=['Quadium','Deuterium','HMR','GROMOS','AMBER','CHARMM']#lista de simulaciones a representar

figs, axsw = plt.subplots(3, 2,sharex=True, sharey=True)#son 6 simulaciones, por lo que 3*2
figs.subplots_adjust(wspace=0, hspace=0)







#figrama.subplots_adjust(wspace=0, hspace=0)



colores=['blue','k','orange']
color_indice=0
for directorio in globales:
    for ff in range(len(lista_de_simulaciones)):
        tiempos,helice=analisis_cluster(str(directorio)+'/'+str(carpetas[ff])+'/conc.pdb',str(directorio)+'/'+str(carpetas[ff])+'/conc.xtc',str(directorio)+'/'+str(carpetas[ff])+'/ref.pdb')
        # Plotear los puntos del RMSD con colores según el cluster
        tiempos=tiempos/1000
        if ff == 0 or  ff == 1:
            axsw[0,ff].plot(tiempos,helice,color=colores[color_indice])
            titulo_2=axsw[0,ff]
        elif ff == 2 or ff == 3:
            axsw[1,ff-2].plot(tiempos,helice,color=colores[color_indice])
            titulo_2=axsw[1,ff-2]
        else:
            axsw[2,ff-4].plot(tiempos,helice,color=colores[color_indice])
            titulo_2=axsw[2,ff-4]
        if color_indice==1:
            titulo_2.tick_params(labelsize=20)
            titulo_2=titulo_2.set_title(str(lista_de_simulaciones[ff]),y=0.6,x=0.5);titulo_2.set_position([.08,0.85]);titulo_2.set_bbox(dict(facecolor='white', edgecolor='black', pad=5.0))
    color_indice=color_indice+1


figs.text(0.5, 0.04, r'Time (ns)', ha='center', fontsize=22)
figs.text(0.07, 0.5, r'%$\mathbf{(\alpha)-helix}$', va='center', rotation='vertical', fontsize=22)
plt.minorticks_on()
custom_handles = [
    plt.Line2D([], [], color='blue', linestyle='-', label='Water'),
    plt.Line2D([], [], color='k', linestyle='-', label='TFE'),
    plt.Line2D([], [], color='orange', linestyle='-', label='TFE Initial Helix'),
]

# Crear la leyenda utilizando los handles personalizados
figs.legend(handles=custom_handles,loc='upper center',ncol=4,fontsize=22)#, bbox_to_anchor=(0.5, 1.15), ncol=3)

plt.show()



figs_hidrogen, axs_hidrogeno = plt.subplots(3, 2,sharex=True, sharey=True)#son 6 simulaciones, por lo que 3*2
figs_hidrogen.subplots_adjust(wspace=0, hspace=0)


import MDAnalysis
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import (
  HydrogenBondAnalysis as HBA)

def calcular_enlaces_hidrogeno(pdb_file,xtc_file):
    u = mda.Universe(pdb_file, xtc_file)
    hbonds = HBA(universe=u,donors_sel='name O or name N',hydrogens_sel='name H or name HN',acceptors_sel='name O or name N',d_a_cutoff=3.6,d_h_cutoff=1.8)
    hbonds.run()
    hbonds_smooth=explicit_heat_smooth(hbonds.count_by_time(),len(hbonds.count_by_time()))

    return hbonds.times,hbonds_smooth
color_indice=0
for directorio in globales:
# Calcular el número de enlaces de hidrógeno a lo largo de la trayectoria
    for ff in range(len(lista_de_simulaciones)):
        tempo,enlaces_H=calcular_enlaces_hidrogeno(str(directorio)+'/'+str(carpetas[ff])+'/conc.pdb',str(directorio)+'/'+str(carpetas[ff])+'/conc.xtc')
        if ff == 0 or  ff == 1:
            axs_hidrogeno[0,ff].plot(tempo/1000,enlaces_H,color=colores[color_indice])
            titulo_1=axs_hidrogeno[0,ff]
        elif ff == 2 or ff == 3:
            axs_hidrogeno[1,ff-2].plot(tempo/1000,enlaces_H,color=colores[color_indice])
            titulo_1=axs_hidrogeno[1,ff-2]
        else:
            axs_hidrogeno[2,ff-4].plot(tempo/1000,enlaces_H,color=colores[color_indice])
            titulo_1=axs_hidrogeno[2,ff-4]
        if color_indice==2:
            titulo_1.tick_params(labelsize=20)
            titulo_1=titulo_1.set_title(str(lista_de_simulaciones[ff]),y=0.6,x=0.5);titulo_1.set_position([.08,0.85]);titulo_1.set_bbox(dict(facecolor='white', edgecolor='black', pad=5.0))
    color_indice=color_indice+1

figs_hidrogen.text(0.5, 0.04, r'Time (ns)', ha='center', fontsize=22)
figs_hidrogen.text(0.07, 0.5, r'# HYDROGEN BONDS', va='center', rotation='vertical', fontsize=22)

figs.legend(handles=custom_handles,loc='upper center',ncol=4,fontsize=22)#, bbox_to_anchor=(0.5, 1.15), ncol=3)

plt.show()
plt.close()
