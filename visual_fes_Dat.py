import numpy as np
from os.path import exists
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.interpolate import interp1d
mpl.use('Agg')
import argparse
# Crear el objeto ArgumentParser
parser = argparse.ArgumentParser(description='Script que procesa archivos de entrada')

# Agregar argumentos para los archivos de entrada
parser.add_argument('-f', '--input1', required=True, help='Ruta del archivo a representar')


from scipy.integrate import quad

def calcular_area_bajo_curva(x, y, y_incert_promedio,x_inicio, x_fin, R=8.314472/1000., T=300):
    # Encontrar los índices correspondientes a x=6 y x=7
    y_incert_promedio = np.mean(y_incert_promedio[:np.argmax(x >= x_fin)])
    indice_fin = np.argmax(x >= x_fin)
    indice_inicio = np.argmax(x >= x_inicio)
    

    y_promedio = np.mean(y[indice_inicio:indice_fin])
    #print('Promedio : ',y_promedio)
    # Restar ese valor de y a toda la gráfica
    y_corregido = y - y_promedio
    # Calcular la nueva función g
    g = np.exp(-y_corregido / (R*T))
    x_area = x[0:indice_fin]
    g_area = g[0:indice_fin]
    #print(g_area)
    # Calcular el área bajo la curva utilizando la regla del trapecio
    area = np.trapz(g_area, x_area)
    '''
    print('AREA:',area)
    # Graficar la función y sombrear el área bajo la curva
    plt.figure()
    plt.plot(x, g, label='g(x)')
    plt.fill_between(x, g, where=(x >= 0) & (x <= x_fin), color='skyblue', alpha=0.3)

    # Calcular el área bajo la curva utilizando la función quad de SciPy
    plt.xlabel('x')
    plt.ylabel('g(x)')
    plt.title('Gráfica de la función g y área sombreada')
    plt.legend()
    plt.grid(True)
    plt.show()
    '''
    incertidumbre_area = np.sqrt(np.sum((y_incert_promedio * np.gradient(g))[:indice_fin] ** 2))
    return area, incertidumbre_area
# Leer los argumentos de la línea de comandos
args = parser.parse_args()
T=298;R=8.314472/1000.
# Acceder a los archivos de entrada
carpeta = args.input1
print('FOLDER: ',str(carpeta))
#directorios = ['CANCER','BACTERIA','POPC']
directorios = ['CANCER','POPE_POPG_1_9','POPC']
nombres = ['Cancer','Bacteria','Mammal']
colors = ['orange','green','blue']
fig = plt.figure()
for archivo in range(len(directorios)):
    fes=open(str(carpeta)+'/'+str(directorios[archivo])+'/'+'fes.dat')
    xx=[]
    yy=[]
    zz=[]
    nfesfiles=0
    for i in range(1,1000): # Detect number of 1D fes files
        if exists("fes_"+str(i)+".dat"): nfesfiles=nfesfiles+1
        else: break
    n1Dplots=int(nfesfiles*10/100)

    correccion=1
    for line in fes:
        if '#' in line:#para saltar as liñas comentadas
            pass
        else:
            if round(float(line.split()[0]))==-6:
                correccion=-1
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
                print('a')
                if round(float(line.split()[0]))==-6:
                    correccion=-1
                x.append(float(line.split()[0]))
                y.append(float(line.split()[1]))
                z.append(float(line.split()[2])) 
        
        plt.xlim(0,30);plt.plot(x*correccion,y,linestyle='--')
        plt.ylim(0,150)
        plt.xlabel('Position (nm)',fontsize=38);print('AREA: ',calcular_area_bajo_curva(x*correccion, y))

    a=np.loadtxt(str(carpeta)+'/'+str(directorios[archivo])+'/'+'CV1_BA.dat')
    energia=a[:,1]
    pos=a[:,0]*correccion
    indices_ordenados = np.argsort(pos)
    pos = np.sort(pos)
    inc=a[:,2]
    #print('Correccion:',correccion)
    energia = [energia[i] for i in indices_ordenados]
    inc = [inc[i] for i in indices_ordenados]


    plt.plot(pos,energia,color=colors[archivo],label=nombres[archivo])
    area,incert_area = calcular_area_bajo_curva(pos, energia, inc, 5,6)
    incert_area = incert_area
    print(f'AREA en {nombres[archivo]}: {round(-R*T*np.log(area/166**(1/3)),3)}+-{abs(round(-R*T*np.log(incert_area/166**(1/3))/10,2)) } kJ/mol +-{round(max(inc),2)} ')#en nm estamos calculando el area
    plt.fill_between(pos, np.array(energia)-np.array(inc),np.array(energia)+np.array(inc),alpha=0.3,color=colors[archivo])
import os    
nombre_carpeta = (str(carpeta))
nombre_carpeta = nombre_carpeta.replace('/', '_')
plt.ylim(0,180)
plt.xlim(-0.5,7.5)
plt.xticks(fontsize=10)
plt.minorticks_on()
plt.yticks(fontsize=10)
plt.ylabel('Energy (kJ/mol)',fontsize=12)
plt.xlabel('Distance (nm)',fontsize=12)
plt.legend(fontsize=10)
#plt.savefig(str(carpeta)+".png", dpi = 300)
nombre_archivo = f"{nombre_carpeta}.png"
print(f"Guardado exitosamente en: {nombre_archivo}, {nombre_carpeta}")
plt.savefig(nombre_archivo, dpi=300)





