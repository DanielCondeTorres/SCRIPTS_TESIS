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

u = mda.Universe('conc.gro', 'conc.xtc')
ref=mda.Universe('ref.pdb')

# Cargar el archivo xtc y pdb con MDAnalysis
# Seleccionar los átomos de interés (por ejemplo, átomos de carbono alfa)
seleccion = u.select_atoms('backbone')
seleccion_ref=ref.select_atoms('backbone')
# Crear una matriz de características basada en el porcentaje de hélice alfa y el radio de giro
caracteristicas = []
tiempos = []
for ts in u.trajectory:
    # Calcular los ángulos dihedrales para determinar la presencia de hélices alfa
    phi = np.degrees(calc_angles(seleccion.residues.atoms[1].position,
                                 seleccion.residues.atoms[0].position,
                                 seleccion.residues.atoms[-1].position))
    psi = np.degrees(calc_angles(seleccion.residues.atoms[-2].position,
                                 seleccion.residues.atoms[-1].position,
                                 seleccion.residues.atoms[1].position))
    
    # Calcular el porcentaje de hélice alfa
    helices = np.logical_and(phi > -150, phi < -30) & np.logical_and(psi > 30, psi < 150)  # Rangos de ángulos para considerar hélices alfa
    porcentaje_helix = np.sum(helices) / float(len(seleccion)) * 100

    # Calcular el radio de giro
    radio_giro = seleccion.radius_of_gyration()
    
    # Agregar las características a la matriz
    caracteristicas.append([porcentaje_helix, radio_giro])
    tiempos.append(ts.time)

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
plt.plot(tiempos,rmsd_valores,color='k')
# Plotear los puntos del RMSD con colores según el cluster
colores = ['r', 'g', 'b', 'c', 'm', 'y']  # Colores para cada cluster
for i in range(n_clusters):
    indices_cluster = np.where(etiquetas_cluster == i)[0]
    plt.scatter(tiempos[indices_cluster], rmsd_valores[indices_cluster], c=colores[i], label=f'Cluster {i+1}')

plt.xlabel('Tiempo (ps)')
plt.ylabel('RMSD (nm)')
plt.legend()
plt.show()
