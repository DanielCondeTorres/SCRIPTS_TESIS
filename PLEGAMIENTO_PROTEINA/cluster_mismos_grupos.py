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

referencias = ['GROMOS_QUAD/ref.pdb', 'GROMOS_DEU/ref.pdb', 'GROMOS_HMR/ref.pdb', 'GROMOS_D/ref.pdb', 'AMBER_D/ref.pdb', 'CHARMM_D/ref.pdb']



import numpy as np
import MDAnalysis as mda
from sklearn.cluster import KMeans
from MDAnalysis.analysis import rms
import matplotlib.pyplot as plt

def analisis_cluster(archivos_gro, archivos_xtc, referencias):
    caracteristicas = []
    rmsd_totales = []
    etiquetas_cluster_totales = []
    tiempos_totales=[]
    for archivo_gro, archivo_xtc, referencia in zip(archivos_gro, archivos_xtc, referencias):
        u = mda.Universe(archivo_gro, archivo_xtc)
        ref = mda.Universe(referencia)
        i=0
        seleccion = u.select_atoms('backbone')
        seleccion_ref = ref.select_atoms('backbone')
        r=RAMACHANDRAM_PLOT(archivo_gro,archivo_xtc)
        tiempos=[]
        for ts in u.trajectory:
            phi=r.results.angles[i][:, 0]
            psi=r.results.angles[i][:, 1]
            i=i+1
            tiempos.append(ts.time)
            helices = np.logical_and(phi > -150, phi < -30) & np.logical_and(psi > -120, psi < 30)
            porcentaje_helix = np.sum(helices) / float(len(helices)) * 100
            radio_giro = seleccion.radius_of_gyration()

            caracteristicas.append([porcentaje_helix, radio_giro])

            rmsd_totales.append(rms.rmsd(seleccion.positions, seleccion_ref.positions, superposition=True))
        tiempos_totales.append(tiempos)		
    caracteristicas = np.array(caracteristicas)
    tiempos_totales=np.array(tiempos_totales)
    rmsd_totales = np.array(rmsd_totales)
	
    n_clusters = 6
    kmeans = KMeans(n_clusters=n_clusters)
    kmeans.fit(caracteristicas)
    etiquetas_cluster_totales = kmeans.labels_

    return rmsd_totales, etiquetas_cluster_totales,tiempos_totales


carpetas = ['GROMOS_QUAD', 'GROMOS_DEU', 'GROMOS_HMR', 'GROMOS_D', 'AMBER_D', 'CHARMM_D']
lista_de_simulaciones = ['Quadium', 'Deuterium', 'HMR', 'GROMOS', 'AMBER', 'CHARMM']

archivos_gro = ['GROMOS_QUAD/conc.pdb', 'GROMOS_DEU/conc.pdb', 'GROMOS_HMR/conc.pdb', 'GROMOS_D/conc.pdb', 'AMBER_D/conc.pdb', 'CHARMM_D/conc.pdb']
archivos_xtc = ['GROMOS_QUAD/conc.xtc', 'GROMOS_DEU/conc.xtc', 'GROMOS_HMR/conc.xtc', 'GROMOS_D/conc.xtc', 'AMBER_D/conc.xtc', 'CHARMM_D/conc.xtc']

fig, axs = plt.subplots(3, 2,sharex=True, sharey=True)#son 6 simulaciones, por lo que 3*2
fig.subplots_adjust(wspace=0, hspace=0)

for i, (gro, xtc, ref) in enumerate(zip(archivos_gro, archivos_xtc, referencias)):
    rmsd, etiquetas_cluster,tiempo = analisis_cluster([gro], [xtc], [ref])
    colores = ['r', 'g', 'b', 'c', 'm', 'y']
    print(i) 
    row = i % 3
    col = i // 3
    tiempo=tiempo.reshape(len(rmsd))/1000
    axs[row, col].set_title(lista_de_simulaciones[i])
    axs[row, col].plot(tiempo, rmsd,color='k')
    for cluster in range(6):
        indices_cluster = np.where(etiquetas_cluster == cluster)[0]
        porcentajes=len(indices_cluster)/len(rmsd)*100
        percentaje=np.round(porcentajes,2)
        axs[row, col].scatter(tiempo[indices_cluster], rmsd[indices_cluster], c=colores[cluster], label=f'Cluster {cluster+1}'+' '+str(percentaje)+'%')
        axs[row,col].legend(fontsize=10,ncol=3,markerscale=2.,loc='lower right');axs[row,col].tick_params(labelsize=18)
        
    print(str(lista_de_simulaciones[i]))
    titulo_1=axs[row, col]
    titulo_1=titulo_1.set_title(str(lista_de_simulaciones[i]),y=0.6,x=0.5);titulo_1.set_position([.08,0.85]);titulo_1.set_bbox(dict(facecolor='white', edgecolor='black', pad=5.0))

fig.text(0.5, 0.04, r'Time (ns)', ha='center', fontsize=22)
fig.text(0.07, 0.5, r'RMSD (${\AA}$)', va='center', rotation='vertical', fontsize=22)

plt.minorticks_on()
plt.show()
