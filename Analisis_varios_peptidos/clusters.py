#Importamos las librerías necesarias
import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from MDAnalysis.analysis import contacts
###########################################################################################################
#                                           FUNCIONES                                                     #
###########################################################################################################
def marksize(my_list,factor):
    counts = {}
    for item in my_list:
        counts[item] = my_list.count(item)
    lista = [counts[x] for x in counts for _ in range(counts[x])]
    return np.array(lista)*factor

#Funcion donde establecemos el numero de contactos
def contacts_within_cutoff(u, group_a, group_b, radius=6.0):
    timeseries = [];time=[]
    for ts in u.trajectory:
        # calculate distances between group_a and group_b
        dist = contacts.distance_array(group_a.positions, group_b.positions)
        # determine which distances <= radius
        n_contacts = contacts.contact_matrix(dist, radius).sum()
        if n_contacts!=0:
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

#Funcion que crea la matriz de contactos
def busqueda_peptidos(numero_de_peptidos,peptide_atoms,cutoff):
    residuos_por_peptido=int(len(peptide_atoms)/numero_peptidos_en_el_sistema)#Solo valido si los peptidos son iguales
    lista_peptidos=[]#lista que devuelve cada peptido del sistema
    matriz_contactos=[]#matriz que nos devuelve la matriz de contactos,1ºColumna: Peptido en el que estamos,2ºColumna peptido al que tocamos,3ºColumna frame (lo cambiaremos por tiempo)
    for buscar_peptidos in range(int(numero_de_peptidos)):
        lista_peptidos.append(peptide_atoms[residuos_por_peptido*buscar_peptidos:residuos_por_peptido*(buscar_peptidos+1)])
    for peptidos_i in ((lista_peptidos)):
        alex_tonto=[]
        for peptidos_j in ((lista_peptidos)):
            calculamos_contactos,times=contacts_within_cutoff(u, peptidos_i, peptidos_j,cutoff)
            alex_tonto.append(calculamos_contactos[:,2])
        matriz_contactos.append(alex_tonto)
    return np.array(matriz_contactos),times

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
    f,c=np.shape(np.array(array))
    lstyle=['--','-.',':','-','--','-.',':','-']
    colors=['blue','red','k','green','orange','purple','yellow','brown']
    for peptides in range(c):
        plt.plot(tiempos_append,list(np.array(array)[:,peptides]),label='Peptide '+str(peptides+1),color=colors[peptides],ls=lstyle[peptides])
    plt.legend()
    #Tiempo ps, Distancias Amstrongs
    plt.xlabel('Time (microseconds)',fontsize=22)
    plt.ylabel('Cluster size',fontsize=22)
    plt.xticks(fontsize=20);plt.yticks(np.arange(0,numero_peptidos_en_el_sistema+2,1),fontsize=20)
    plt.show()

def representacion2(array,tiempos,numero_peptidos_en_el_sistema):
    f,c=np.shape(np.array(array))
    lstyle=['--','-.',':','-','--','-.',':','-']
    colors=['blue','red','k','green','orange','purple','yellow','brown']
    for peptides in range(c):
        plt.plot(tiempos_append,list(np.array(array)[:,peptides]),label='Cluster size '+str(peptides+1),color=colors[peptides],ls=lstyle[peptides],marker='o')
    plt.legend()
    #Tiempo ps, Distancias Amstrongs
    plt.xlabel('Time (microseconds)',fontsize=22)
    plt.ylabel('Peptides in the cluster size',fontsize=22)

    plt.ylim(0,9)
    plt.xticks(fontsize=20);plt.yticks(np.arange(0,numero_peptidos_en_el_sistema+2,1),fontsize=20)
    plt.show()


###########################################################################################################
#                                       FUNCIONES FIN                                                     #
###########################################################################################################


###########################################################################################################
#                                       MAIN PROGRAM                                                      #
###########################################################################################################
#Inputs:
# Carga el archivo PDB y XTC
factor_marker=100
u = mda.Universe('traj_0.pdb', 'borrar_dani.xtc')
numero_peptidos_en_el_sistema=8
cutoff=6
peptide_atoms = u.select_atoms('name BB')
#No more inputs
print('Loading... Alex, you can take a coffe, Dani trabaja por ti <3')
# Selecciona todos los átomos de cada péptido
matriz_contactos,times=busqueda_peptidos(numero_peptidos_en_el_sistema,peptide_atoms,cutoff)
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
    #plt.scatter([tiempos]*len(cluster_lista),cluster_lista,s=size,color='blue')
    Peptido.append(almacenar_valores_peptido_contacto(dictado,float(tiempos[0])))
    tiempos_append.append(float(tiempos[0]))
grafica=representacion(Peptido,tiempos_append,numero_peptidos_en_el_sistema)
grafica=representacion2(peptidos_por_cluster,tiempos_append,numero_peptidos_en_el_sistema)
###########################################################################################################
