import os
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
mpl.use('Agg')


import numpy as np
import matplotlib.pyplot as plt
############################################################################################################
#                                 FUNCIÓN PARA O CÁLCULO DO BA                                             #
############################################################################################################
def blockAverage(datos,tamaño_bloque_maximo=False,grafica=False):
    "Método para calcular incertidumes con datos correlacionados temporalmente, block-avarage. Sería ideal engadir o método de bootstraping para comparar"
    Nobs = len(datos) # número total de datos, dar un array ou lista
    tamaño_bloque_minimo=1 # mínimo tamaño de bloque, corresponde a non facer block_avarage
    if tamaño_bloque_maximo==False:
        tamaño_bloque_maximo=int(Nobs/4)#Criterio por se non se selecciona un dato óptimo de bloque
    Numero_de_bloques = tamaño_bloque_maximo-tamaño_bloque_minimo # total number of block sizes
    Media_bloque = np.zeros(Numero_de_bloques) 
    co = np.zeros(Numero_de_bloques) # definese no seguinte papper, pero é a varianza
    #print('A información do método atópase en:  https://doi.org/10.1063/1.457480')
    bloqueCtr=0
    incertidume=np.zeros(Numero_de_bloques)
    barra_erros=np.zeros(Numero_de_bloques)
    for datos_bloque in range(tamaño_bloque_minimo,tamaño_bloque_maximo):
        Nbloque=int(Nobs/datos_bloque) # En que bloque estou según o numero de observacions
        almacenamento_temporal=np.zeros(Nbloque) # isto é para almacenar a media según en que bloque esté, é temporal
        # Loop to chop datastream into blocks
        # and take average
        for i in range(1,Nbloque+1):
            #Esta é a parte máis complicada do código pois (fago de 1 a N+1, para empregar logo as definicións matemáticas)
            #E como vou calculando as medias según o bloque no que me atope, e dicir o número de datos por bloque
            comeza = (i-1)*datos_bloque
            remata = comeza+datos_bloque
            almacenamento_temporal[i-1] = np.mean(datos[comeza:remata])
        Media_bloque[ bloqueCtr] = np.mean(almacenamento_temporal)
        co[bloqueCtr]  = np.var(almacenamento_temporal)
        incertidume[ bloqueCtr]=np.sqrt(np.var(almacenamento_temporal)/((Nbloque-1)))
        barra_erros[ bloqueCtr]=np.sqrt(np.var(almacenamento_temporal)/(Nbloque-1))*1/(np.sqrt(2*(Nbloque-1)))
        bloqueCtr=bloqueCtr+1
    tamaño_bloque= np.arange(tamaño_bloque_minimo,tamaño_bloque_maximo)
    if grafica == True:
        plt.errorbar(tamaño_bloque, incertidume,barra_erros, marker='o',ls='None',markersize=8, capsize=6)
        plt.xlabel('Tamaño do bloque (Número de datos por bloque, que deben ter todos o mesmo tamaño)')
        plt.ylabel('Desviación estándar')
        plt.show()
        plt.errorbar(tamaño_bloque, Media_bloque, incertidume,marker='o',ls='None',markersize=8, capsize=6,color='Orange')
        plt.ylabel('Media')
        plt.xlabel('Tamaño do bloque')
        plt.show()
        plt.tight_layout()
        plt.show()
        #print(f"Valor medio = {Media_bloque[-1]:.3f} +/- {incertidume[-1]:.3f}")
    return tamaño_bloque, Media_bloque, incertidume
#######################################################

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

# Lista de nombres de las carpetas
folder_names = ['700']#,'CA_300_EMPIEZA_MEMBRANA']

# Crear la figura y los ejes de subparcelas
fig, axs = plt.subplots(4, 2, sharex=True, sharey=True, figsize=(10, 10))
fig.subplots_adjust(wspace=0, hspace=0)
contador=0
# Leer los archivos y almacenar los datos en la lista
for folder_idx, folder in enumerate(folder_names):
    # Lista para almacenar los datos
    lista_valor_medio = []
    lista_incertidumbre = []
    ba_media = []
    ba_inc = []
    # Lista para almacenar los datos
    data = [[] for _ in range(8)]
    folder_path = os.path.join(folder)

    for i in range(8):
        print('Reading file:', f'num_bonds_{i+1}.xvg')
        filename = os.path.join(folder_path, f'num_bonds_{i+1}.xvg')
        with open(filename, 'r') as file:
            lines = file.readlines()
            # Ignorar la primera línea si contiene encabezados
            if lines[0].startswith('x'):
                lines = lines[1:]
            # Separar las columnas y convertirlas en listas de flotantes
            x_values = []
            y_values = []
            for line in lines:
                line = line.strip()
                if line.startswith('#') or line.startswith('@'):
                    continue  # Saltar líneas que empiezan con # o @
                values = line.split()
                if len(values) >= 2:
                    x_values.append(float(values[0]))
                    y_values.append(float(values[1]))
            if x_values and y_values:
                data[i].extend((x_values, y_values))
        print('Files readed')
        # Configurar los colores de las líneas para cada gráfica en subparcelas
        #colors = ['blue', 'orange', 'green','black']#,'purple']
        import matplotlib.colors as mcolors

        # Definir los códigos hexadecimales
        color_codes = ['purple']
                        
#['#feedde', '#fdd0a2', '#fdae6b', '#fd8d3c', '#e6550d', '#a63603']
        colors = [mcolors.hex2color(code) for code in color_codes]
        # Convertir a valores RGB
        # Plotear los datos en subparcelas separadas
        for j, values in enumerate(data):
            row = j // 2
            col = j % 2
            ax = axs[row, col]
            if values:
                x = np.array(values[0])*10**(-6)  # Convertir a microsegundos
                y = values[1]
                indices = [i for i, valor in enumerate(x) if valor > 5]
                # Filtrar los elementos de la lista1 cuyos índices coinciden con los valores mayores que 5 en la lista2
                medias = [y[i] for i in indices]

                lista_valor_medio.append(np.mean(medias))
                lista_incertidumbre.append(np.std(medias))

                #BA
                #tamaño_bloque, Media_bloque, incertidume = blockAverage(medias,tamaño_bloque_maximo=False,grafica=False)

                #ba_media.append(Media_bloque[-1])
                #ba_inc.append(incertidume[-1])
                #x=explicit_heat_smooth(np.array(values[0]),len(values[0]));y=explicit_heat_smooth(np.array(values[1]),len(values[1]))
                print('Media movil')
                n=800;x=media_movil(values[0],n);y=media_movil(values[1],n)
                x=np.array(x)*10**(-6) 
                if contador ==0 or contador ==3:
                    ax.plot(x, y, color=colors[contador],ls='-')
                elif contador!=4:
                    ax.plot(x, y, color=colors[contador],ls='-')  # Usar colors[i % len(colors)] para alternar los colores
                else:
                    ax.plot(x, y, color=colors[contador],ls='-')
                ax.tick_params(labelsize=14)
                ax.set_xlim(-0.3,10)
                ax.set_ylim(0,90)
                ax.set_xticks([1, 3, 5, 7, 9])
            titulo=ax.set_title( f'Peptide {j+1}',y=0.6,x=0.85,fontsize=12);titulo.set_position([.85,0.77]);titulo.set_bbox(dict(facecolor='white', edgecolor='black', pad=5.0))
    contador=contador+1
    print('Carpeta: ', folder)
    print('LEN: ',len(lista_valor_medio))
    print('Media: ', np.sum(lista_valor_medio))
    print('Inc: ', np.sqrt(sum(x**2 for x in lista_incertidumbre)))
    print(f"Valor medio = {np.sum(ba_media):.4f} +/- {1.96*np.sqrt(sum(x**2 for x in ba_inc)):.3f}")

# Ajustar los espacios entre las subparcelas

# Mostrar la gráfica
# Título de los ejes:
fig.text(0.5, 0.04, r'Time (' r'$\mu$s)', ha='center', fontsize=18)
fig.text(0.05, 0.5, r'# Contacts (BB-Water)', va='center', rotation='vertical', fontsize=18)
plt.minorticks_on()

##   f'# Ca$^{{2+}}$/#Total Ions: {ca_ions}/{total_ions}'
custom_handles = [
       # plt.Line2D([], [], color=color_codes[0] , linestyle='-', label=f' 0 M'),# ([Ca$^{{2+}}$])'),#'# Ca$^{{2+}}$/#Total Ions: 47/47'),
        #plt.Line2D([], [], color=color_codes[0] , linestyle='-', label=f'Reference: 0.86% (Na$^{{+}}$-Cl$^{{-}}$)'),#'# Na$^{{+}}$/#Total Ions: 94/94'),
    #plt.Line2D([], [], color=color_codes[1] , linestyle='-', label=f'0.43% (Ca$^{{2+}}$-Cl$^{{-}}$)'),#'# Ca$^{{2+}}$/#Total Ions: 47/47'),
   # plt.Line2D([], [], color=color_codes[1] , linestyle='-', label=f' 0.117 M '),#([Ca$^{{2+}}$])'),#'# Ca$^{{2+}}$/#Total Ions: 200/506'),
    #plt.Line2D([], [], color=color_codes[4] , linestyle='-', label=f'7.37% (Ca$^{{2+}}$-Cl$^{{-}}$)'),#'# Ca$^{{2+}}$/#Total Ions: 300/806'),
   # plt.Line2D([], [], color=color_codes[2] , linestyle='-', label=f' 0.263 M'),# ([Ca$^{{2+}}$])'),#'# Ca$^{{2+}}$/#Total Ions: 400/1106'),
    plt.Line2D([], [], color=color_codes[0], linestyle='-', label='0.485 M')# ([Ca$^{{2+}}$])')
]

# Crear la leyenda utilizando los handles personalizados
fig.legend(handles=custom_handles,loc='upper center',ncol=4,fontsize=12)#, bbox_to_anchor=(0.5, 1.15), ncol=3)
plt.show()
plt.savefig("contacts_water.png", dpi = 300)

