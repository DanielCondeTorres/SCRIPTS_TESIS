#Importamos las librerías necesarias
import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from MDAnalysis.analysis import contacts
###################################
# LECTURA DE ARCHIVOS             #
###################################
import argparse

# Crear el objeto ArgumentParser
parser = argparse.ArgumentParser(description='Script que procesa archivos de entrada')

# Agregar argumentos para los archivos de entrada
parser.add_argument('-f', '--input1', required=True, help='Ruta del archivo .pdb')
parser.add_argument('-x', '--input2', required=True, help='Ruta del archivo .xtx')

# Leer los argumentos de la línea de comandos
args = parser.parse_args()

# Acceder a los archivos de entrada
archivo_pdb = args.input1
archivo_xtc = args.input2






###########################################################################################################
#                                           FUNCIONES                                                     #
###########################################################################################################
def marksize(my_list,factor):
    counts = {}
    for item in my_list:
        counts[item] = my_list.count(item)
    lista = [counts[x] for x in counts for _ in range(counts[x])]
    return np.array(lista)*factor

def media_movil(x,n):
    resultado=[]
    for i in range(n):
        if i==0:
            pass
        else:
            resultado.append(float((np.sum(x[0:i])/len(x[0:i]))))
    cumsum=np.cumsum(np.insert(x, 0, 0))
    res=(cumsum [n:]-cumsum[:-n])/float(n)
    for elemento in res:
        resultado.append(elemento)
    return resultado

#Funcion donde establecemos el numero de contactos con un cutoff de 6 Amstrings
def contacts_within_cutoff(u, group_a, group_b, radius=6.0,contactos_minimos=16):
    timeseries = [];time=[]
    for ts in u.trajectory[::50]:
        # calculate distances between group_a and group_b
        dist = contacts.distance_array(group_a.positions, group_b.positions)
        # determine which distances <= radius
        n_contacts = contacts.contact_matrix(dist, radius).sum()
        if n_contacts>contactos_minimos: #cutoff contactos minimos (!=0; criterio inicial) 
            variable=1
        else:
            variable=0
        if group_a==group_b:
            timeseries.append([ts.frame,0,0])
            time.append(u.trajectory.time)
        else:
            timeseries.append([ts.frame, n_contacts,variable])
            time.append(u.trajectory.time)
    return np.array(timeseries),time

#Función con la calculamos los clusters
#Funcion calculo de cluster:
def calculamos_clusters(contactos,numero_bolas):
    for bola in range(1, numero_bolas):
        cluster = set()
        visitados = set()
        cola = [bola]
        cluster.add(bola)#Tenemos en cuenta la actual
        while cola:
            actual = cola.pop(0)
            visitados.add(actual)
            for contacto in contactos[actual]:
                if contacto not in visitados:
                    cluster.add(contacto)
                    cola.append(contacto)
        #print(f"Bola {bola}: {cluster}")
        contactos[bola].update(cluster)
    return contactos
def contactos_a_calcular_peptido_membrana(numero_de_peptidos,peptide_atoms,lipidos,cutoff):
    residuos_por_peptido=int(len(peptide_atoms)/numero_peptidos_en_el_sistema)
    lista_peptidos=[]
    matriz_contactos=[]
    for buscar_peptidos in range(int(numero_de_peptidos)):
        lista_peptidos.append(peptide_atoms[residuos_por_peptido*buscar_peptidos:residuos_por_peptido*(buscar_peptidos+1)])
    for peptidos_i in ((lista_peptidos)):
        calculamos_contactos=[]
        print('LIPIDOS',lipidos,cutoff)
        calculamos_contactos,times=contacts_within_cutoff(u, peptidos_i.select_atoms('name BB'), lipidos,cutoff,0)
        matriz_contactos.append(np.array(calculamos_contactos[:,1]))
    print(matriz_contactos)
    return matriz_contactos

def contactos_totales(numero_de_peptidos,peptide_atoms,cutoff):
    residuos_por_peptido=int(len(peptide_atoms)/numero_peptidos_en_el_sistema)
    lista_peptidos=[]
    matriz_contactos=[]
    for buscar_peptidos in range(int(numero_de_peptidos)):
        lista_peptidos.append(peptide_atoms[residuos_por_peptido*buscar_peptidos:residuos_por_peptido*(buscar_peptidos+1)])
    contador_i=0
    for peptidos_i in ((lista_peptidos)):
        calculamos_contactos=[];
        contador_j=0
        for peptidos_j in ((lista_peptidos)):
            if contador_j>=contador_i:
                calculamos_contactos,times=contacts_within_cutoff(u, peptidos_i.select_atoms('name BB'), peptidos_j.select_atoms('name BB'),cutoff,0)
                contador_j+=1
            else:
                contador_j+=1
            try:
                matriz_contactos.append((calculamos_contactos[:,1]))
            except TypeError:
                continue

        contador_i+=1
    return np.sum(matriz_contactos,axis=0)


#Funcion que crea la matriz de contactos
def busqueda_peptidos(numero_de_peptidos,peptide_atoms,cutoff):
    residuos_por_peptido=int(len(peptide_atoms)/numero_peptidos_en_el_sistema)#Solo valido si los peptidos son iguales
    lista_peptidos=[]#lista que devuelve cada peptido del sistema
    matriz_contactos=[]#matriz que nos devuelve la matriz de contactos,1ºColumna: Peptido en el que estamos,2ºColumna peptido al que tocamos,3ºColumna frame (lo cambiaremos por tiempo)
    matriz_contactos_peptidos=[]
    for buscar_peptidos in range(int(numero_de_peptidos)):
        lista_peptidos.append(peptide_atoms[residuos_por_peptido*buscar_peptidos:residuos_por_peptido*(buscar_peptidos+1)])
    for peptidos_i in ((lista_peptidos)):
        alex_tonto=[];contactos_peptidos=[]
        for peptidos_j in ((lista_peptidos)):
            calculamos_contactos,times=contacts_within_cutoff(u, peptidos_i, peptidos_j,cutoff)
            calculamos_contactos_BB,times=contacts_within_cutoff(u, peptidos_i.select_atoms('name BB'), peptidos_j.select_atoms('name BB'),cutoff)
            alex_tonto.append(calculamos_contactos[:,2]);contactos_peptidos.append(calculamos_contactos_BB[:,1])
        matriz_contactos.append(alex_tonto);matriz_contactos_peptidos.append(sum(contactos_peptidos))
    return np.array(matriz_contactos),times,(matriz_contactos_peptidos)

def quien_toca_a_cada_peptido_por_cada_frame(peptido_situado,frame,matriz_contactos):
    #Calculamos los contactos
    numero_contactos_peptido_frame=np.zeros((peptido_situado,frame))
    quien_toca_al_peptido_por_frame=np.zeros((peptido_situado,frame),dtype=object)
    for elemento in range(peptido_situado):
        for tiempo in range(frame):
            con_quien_toca_por_tiempo=[]
            numero_contactos_peptido_frame[elemento][tiempo]=np.sum(matriz_contactos[elemento,:,tiempo])
            for con_quien_toca in range(peptido_con_el_que_contactamos):
                #quien_toca_al_peptido_por_frame[elemento][tiempo]
                toca=matriz_contactos[elemento,con_quien_toca,tiempo]
                if toca!=0:
                    con_quien_toca_por_tiempo.append(con_quien_toca+1)
                else: 
                    continue
            quien_toca_al_peptido_por_frame[elemento][tiempo]=con_quien_toca_por_tiempo    
    return quien_toca_al_peptido_por_frame






#Buscar valores equivalentes
def find_keys(dic):
    result = {}
    for key, value in dic.items():
        # Convertir los sets en tuplas
        value = tuple(value)
        # Usar setdefault con la tupla
        result.setdefault(value, []).append(key)
    return [keys for keys in result.values() if len(keys) > 1]

def unique_values(dic):
    counts = {}
    for value in dic.values():
        for item in value:
            counts[item] = counts.get(item, 0) + 1
    return [key for key, count in counts.items() if count == 1]



def almacenar_valores_peptido_contacto(dic,tiempo):
    llaves=[];valores=[]
    for key, value in dic.items():
        valores.append(len(value))
    return valores

def representacion(array,tiempos,numero_peptidos_en_el_sistema):
    plt.figure(figsize=(20, 10))
    f,c=np.shape(np.array(array))
    lstyle=['--','-.',':','-','--','-.',':','-']
    colors=['blue','red','k','green','orange','purple','gold','brown']
    for peptides in range(c):
        plt.plot(tiempos_append,list(np.array(array)[:,peptides]),label='Peptide '+str(peptides+1),color=colors[peptides],ls=lstyle[peptides])
    plt.legend(fontsize=20, bbox_to_anchor=(1.00,0.5), loc = "center left")
    #Tiempo ps, Distancias Amstrongs
    plt.title('Peptide clusters', fontsize=40)
    plt.xlabel('Time ('r'$\mu$s)',fontsize=30)
    plt.ylabel('Cluster size',fontsize=30)
    plt.xticks(fontsize=30);plt.yticks(np.arange(0,numero_peptidos_en_el_sistema+2,1),fontsize=30)
    plt.tight_layout()
    plt.savefig("cluster_lineas_dashed.png", dpi = 300)
    plt.show()

def representacion2(array,tiempos,numero_peptidos_en_el_sistema):
    plt.figure(figsize=(20, 10))
    f,c=np.shape(np.array(array))
    lstyle=['--','-.',':','-','--','-.',':','-']
    colors=['blue','red','k','green','orange','purple','gold','brown']
    for peptides in range(c):
        y = media_movil(list(np.array(array)[:,peptides]),1)
        plt.plot(tiempos_append,y,label=str(peptides+1) + ' Peptide(s) cluster', color=colors[peptides])#+str(peptides+1),color=colors[peptides])#,ls=lstyle[peptides],marker='o')
    plt.legend(fontsize=20, bbox_to_anchor=(1.00,0.5), loc="center left")
    #Tiempo ps, Distancias Amstrongs
    plt.title('Peptide clusters', fontsize=40)
    plt.xlabel('Time ('r'$\mu$s)',fontsize=30)
    plt.ylabel('# of peptides in each cluster',fontsize=30)
    plt.ylim(-0.5,numero_peptidos_en_el_sistema+1)
    plt.xticks(fontsize=30);plt.yticks(np.arange(0,numero_peptidos_en_el_sistema+2,1),fontsize=30)
    plt.tight_layout()
    plt.savefig("cluster_lineas.png", dpi = 300)
    plt.show()

import matplotlib.pyplot as plt
import numpy as np

def representacion(array, tiempos, numero_peptidos_en_el_sistema):
    f, c = np.shape(np.array(array))
    lstyle = ['--', '-.', ':', '-', '--', '-.', ':', '-']
    colors = ['blue', 'red', 'k', 'green', 'orange', 'purple', 'gold', 'brown']
    
    for peptides in range(c):
        plt.plot(tiempos, list(np.array(array)[:, peptides]),
                 label='Peptide ' + str(peptides + 1),
                 color=colors[peptides], ls=lstyle[peptides])
    
    plt.legend(fontsize=20, bbox_to_anchor=(1.0, 0.5), loc="center left")
    plt.title('Peptide clusters', fontsize=30)
    #plt.xlabel('Time (' r'$\mu$s)', fontsize=30)
    plt.ylabel('Cluster size', fontsize=30)
    plt.yticks(np.arange(0,numero_peptidos_en_el_sistema+2,1),fontsize=10)
    plt.tick_params(axis='x', labelsize=20)
    plt.tick_params(axis='y', labelsize=20)

def representacion2(array, tiempos, numero_peptidos_en_el_sistema):
    f, c = np.shape(np.array(array))
    lstyle = ['--', '-.', ':', '-', '--', '-.', ':', '-']
    colors = ['blue', 'red', 'k', 'green', 'orange', 'purple', 'gold', 'brown']
    
    for peptides in range(c):
        y = media_movil(list(np.array(array)[:, peptides]), 1)
        plt.plot(tiempos, y, label=str(peptides + 1) + ' Peptide(s) cluster',
                 color=colors[peptides])
    
    plt.legend(fontsize=20, bbox_to_anchor=(1.0, 0.5), loc="center left")
    #plt.title('Peptide clusters', fontsize=30)
    plt.xlabel('Time (' r'$\mu$s)', fontsize=30)
    plt.ylabel('# of peptides\nin each cluster', fontsize=30)
    plt.ylim(-0.5, numero_peptidos_en_el_sistema + 1)
    plt.yticks(np.arange(0,numero_peptidos_en_el_sistema+2,1),fontsize=10)
    plt.tick_params(axis='x', labelsize=20)
    plt.tick_params(axis='y', labelsize=20)


###########################################################################################################
#                                       FUNCIONES FIN                                                     #
###########################################################################################################


###########################################################################################################
#                                       MAIN PROGRAM                                                      #
###########################################################################################################
#Inputs:
# Carga el archivo PDB y XTC
factor_marker=100
u = mda.Universe(str(archivo_pdb), str(archivo_xtc))

numero_peptidos_en_el_sistema=int(input("Por favor, ingresa el número de Péptidos en el sistema: "))

cutoff=6

peptide_atoms = u.select_atoms('protein')
lipidos=u.select_atoms('not protein')
#No more inputs
print('Loading... Alex, you can take a coffe, Dani trabaja por ti <3')
# Selecciona todos los átomos de cada péptido
matriz_contactos,times,contactos_peptidos=busqueda_peptidos(numero_peptidos_en_el_sistema,peptide_atoms,cutoff)










peptido_situado,peptido_con_el_que_contactamos,frame=np.shape(matriz_contactos)
#Calculamos los contactos
quien_toca_al_peptido_por_frame=quien_toca_a_cada_peptido_por_cada_frame(peptido_situado,frame,matriz_contactos)
cluster_final=[];tiempos_final=[]
lista=list(range(1,numero_peptidos_en_el_sistema+1,1))

peptidos_por_cluster=[]
tiempos_append=[]
Peptido=[]
for elemento in range(frame):
    donde_estoy=0
    clusters=[];tiempos=[];zz=[]
    diccionario_contactos = {}
    for p in lista:
    	diccionario_contactos[p]=set()
    for te_toca in quien_toca_al_peptido_por_frame[:,elemento]:
        for otra_bola in te_toca:
            diccionario_contactos[donde_estoy+1].add(otra_bola)
        donde_estoy+=1
        clusters.append(len(te_toca))
    tiempos.append(times[elemento]*10**(-6))#poner tiempos
    dictado=calculamos_clusters(diccionario_contactos,numero_peptidos_en_el_sistema+1)
    tiempos_final.append(tiempos)
    cluster_final.append(clusters)
    cluster_lista=[]
    for lista_keys in find_keys(dictado):
        cluster_lista.append(len(lista_keys))
    for lista_keys in unique_values(dictado):
        cluster_lista.append(1)
    tamaño_marker=marksize(cluster_lista,factor_marker)
    for rr in range(1,numero_peptidos_en_el_sistema+1):#contamos los peptidos por cluster
        cluster_lista.count(rr);zz.append(cluster_lista.count(rr)*(rr))#por el valor minimo del cluster (rr)
    peptidos_por_cluster.append(zz)
    primero=len(cluster_lista);segundo=len(tamaño_marker);cociente=int(primero/segundo)
    size=tamaño_marker
    Peptido.append(almacenar_valores_peptido_contacto(dictado,float(tiempos[0])))
    tiempos_append.append(float(tiempos[0]))
print('Descargando virus ...')
# Create the subplots
fig, axs = plt.subplots(2, 1, figsize=(20, 20))#,shstex=True gridspec_kw={'hspace': 0})
#fig.subplots_adjust(hspace=0)
# Plot graph 1 in the first subplot
plt.sca(axs[0])
representacion(Peptido, tiempos_append, numero_peptidos_en_el_sistema)

# Plot graph 2 in the second subplot
plt.sca(axs[1])
representacion2(peptidos_por_cluster, tiempos_append, numero_peptidos_en_el_sistema)

# Adjust the spacing between subplots
plt.tight_layout(pad=10.0)

# Save and display the figure
plt.savefig("cluster_subplots.png", dpi=300)
plt.show()
plt.close()
###########################################################################################################
print('Creando figuras ...')




# Crear una figura con una cuadrícula de 4 filas y 2 columnas
fig, axs = plt.subplots(4, 2,sharex=True, sharey=True, figsize = (20,10))
fig.subplots_adjust(wspace=0, hspace=0)
fig.suptitle("Contacts Peptide (1to1)", fontsize=22)
for i in range(len(contactos_peptidos)):
    if i < int(len((contactos_peptidos))/2):
        axs[i,0].plot(tiempos_append,contactos_peptidos[i])
        axs[i,0].tick_params(labelsize=22)
    else:
        axs[i-4,1].plot(tiempos_append,contactos_peptidos[i])
        axs[i-4,1].tick_params(labelsize=22)
#Titulo de los ejes:
fig.text(0.5, 0.04, r'Time (' r'$\mu$s)', ha='center', fontsize=30)
fig.text(0.07, 0.5, r'# Contacts(BB)', va='center', rotation='vertical', fontsize=30)
plt.minorticks_on()
plt.savefig("Num_contacts_1toAll.png", dpi = 300)
#plt.show()


contactos_membrana=contactos_a_calcular_peptido_membrana(numero_peptidos_en_el_sistema,peptide_atoms,lipidos,cutoff)

# Crear una figura con una cuadrícula de 4 filas y 2 columnas
fig, axs = plt.subplots(4, 2,sharex=True, sharey=True, figsize = (20,10))
fig.subplots_adjust(wspace=0, hspace=0)
fig.suptitle("Contacts Peptide-Membrana", fontsize=22)
for i in range(len(contactos_membrana)):
    if i < int(len((contactos_membrana))/2):
        axs[i,0].plot(tiempos_append,contactos_membrana[i])
        axs[i,0].tick_params(labelsize=22)
    else:
        axs[i-4,1].plot(tiempos_append,contactos_membrana[i])
        axs[i-4,1].tick_params(labelsize=22)

#Titulo de los ejes:
fig.text(0.5, 0.04, r'Time (' r'$\mu$s)', ha='center', fontsize=30)
fig.text(0.07, 0.5, r'# Contacts(BB)', va='center', rotation='vertical', fontsize=30)
plt.minorticks_on()
plt.savefig("Num_contacts_Membrana_lipido.png", dpi = 300)


print('Descargando Counter Strike ...')
fig = plt.figure()
PRUEBA=contactos_totales(numero_peptidos_en_el_sistema,peptide_atoms,cutoff)
plt.plot(tiempos_append,PRUEBA)
plt.xlabel('Time (' r'$\mu$s)',fontsize=14)
plt.ylabel('Total contacts',fontsize=14)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.minorticks_on()
plt.savefig("CONTACTOS_TOTALES.png", dpi = 300)




