import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import griddata
import matplotlib.cm as cm
from os.path import exists
from matplotlib.colors import LogNorm
import matplotlib as mpl
from matplotlib.colors import Normalize
from scipy.spatial import ConvexHull #volumen calcular
from matplotlib.colors import DivergingNorm
from vedo import *
#import plotly.express as px
#import plotly.graph_objects as go
#ctes:
T=298;R=8.314472/1000.
################################################################################
#                   FUNCIÓN PARA EL CÁLCULO DEL ÁREA                           #
################################################################################
from sklearn.metrics import roc_curve,auc
def area(x,y,xmin,xmax,globalmin,T=298,R=8.314472/1000.):#listas de puntos en las coordenadas x e y, asi como el rango en x max y min
    x_corregida=[]#listas vacias para quedarnos con los elementos en x e y de interés
    y_corregida=[]
    y_corregida_exp=[]
    maxima=np.exp(max(-y)/(T*R))
    print('XMAX: ',xmax,'XMIN: ',xmin)
    if xmin<xmax:
        xmin2=xmin
        xmax2=xmax
    else:
        xmin2=xmax
        xmax2=xmin
    for i in range(len(x)):
        if (x[i])>(xmin2) and (x[i])<(xmax2):
            x_corregida.append(x[i])
            y_corregida.append(y[i])
    #return x_corregida,y_corregida#max(y_corregida)-y_corregida
    for elemento in y_corregida:
        y_corregida_exp.append(np.exp(-(elemento-globalmin)/(T*R)))
    norma=np.sum(y_corregida_exp)#menos
    areas=auc(x_corregida,y_corregida_exp)#maxima-y_corregida_exp
    return x_corregida,y_corregida_exp,areas#maxima-y_corregida_exp,areas
################################################################################
######################## INPUTS ################################################
################################################################################
#Aquí se seleccionan los valores máximos y mínimos de las dos CVS que se esten empleando CV1, CV2
cv1min=0
cv1max=10
cv2min=0.45
cv2max=2.5
#Aquí los valores máximos y mínimos del perfil energético
pmfmax=250;factor=pmfmax/abs(pmfmax)
pmfmin=0
# Read 2D fes file
cambio_de_ejes=1 #1 o -1







df=pd.DataFrame(pd.read_csv('fesd1d2.dat',delim_whitespace=True, header=None, skiprows=9))

#Variables x, y e z del archivo fesd1d2.dat
x=df[df.columns[0]].to_numpy()*cambio_de_ejes#Este -1, depende de la definición inicial, igual en tu caso se puede borrar
y=df[df.columns[1]].to_numpy()
z=df[df.columns[2]].to_numpy()
globalmin=0#;globalmin_p=108
z=z-globalmin
#Creamos un gradiente para a representación gráfica
xi = np.linspace(np.min(x), np.max(x), 1000)
yi = np.linspace(np.min(y), np.max(y), 1000)
zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method='linear')
#contamos el número de 1D fes files
nfesfiles=0
for i in range(1,1000): # Detect number of 1D fes files
    if exists("fesd1_"+str(i)+".dat"): nfesfiles=nfesfiles+1
    else: break
n1Dplots=int(nfesfiles*10/100)#Representamos el 30% de los últimos valores
xd1=[];yd1=[];xd2=[];yd2=[];yd1_trans=[];yd2_trans=[]
for i in range(nfesfiles,nfesfiles-n1Dplots,-1): # Read last 10 1D fes files
    df1Dd1=pd.DataFrame(pd.read_csv("fesd1_"+str(i)+".dat",delim_whitespace=True, header=None, skiprows=5))
    df1Dd2=pd.DataFrame(pd.read_csv("fesd2_"+str(i)+".dat",delim_whitespace=True, header=None, skiprows=5))#-df[df.columns[2]].to_numpy()/(R*T)
    xd1.append(df1Dd1[df1Dd1.columns[0]].to_numpy()); yd1.append(df1Dd1[df1Dd1.columns[1]].to_numpy()-globalmin);yd1_trans.append(np.exp(-(df1Dd1[df1Dd1.columns[1]].to_numpy()-globalmin)/(R*T))/np.sum(np.exp(-(df1Dd1[df1Dd1.columns[1]].to_numpy()-globalmin)/(R*T))))#OJO -
    xd2.append(df1Dd2[df1Dd2.columns[0]].to_numpy()); yd2.append(df1Dd2[df1Dd2.columns[1]].to_numpy()-globalmin);yd2_trans.append(np.exp(-(df1Dd2[df1Dd2.columns[1]].to_numpy()-globalmin)/(R*T))/np.sum(np.exp(-(df1Dd2[df1Dd2.columns[1]].to_numpy()-globalmin)/(R*T))))#OJO -


#Representación gráfica
fig, axs = plt.subplots(nrows=2, ncols=2, gridspec_kw={'width_ratios': [3, 1], 'height_ratios': [1, 3],'hspace':0.15,'wspace':0.1})
axs = axs.flatten()
for i in range(1,n1Dplots):
    axs[0].plot(xd1[i],yd1[i])#O -1 depende da tua definición de coordenada x.
    axs[3].plot(yd2[i], xd2[i])
axs[1].axis('off')
axs[3].yaxis.tick_right()
levels=list(range(0, pmfmax, 5))#Niveles a representarlevels=np.sort(levels)
#cpf = axs[2].contourf(xi, yi, zi,100,vmin=min(levels),vmax=max(levels),extend='max',cmap=cm.rainbow.reversed(),ticks=levels)
cpf = axs[2].contourf(xi, yi, zi,100,vmin=min(levels),vmax=max(levels),extend='max',cmap=cm.rainbow.reversed(),ticks=levels,levels=levels)
cpf.changed()
#Representamos la barra del gradiente de colores
cbar =fig.colorbar(cpf,ax=axs[3],ticks=levels)
cbar.ax.tick_params(labelsize=18) 
cpf.changed()
#cbar.ax.set_ylim(levels[0], levels[-1])
cp = axs[2].contour(xi, yi, zi, levels, colors='black', linestyles='dashed')
axs[2].clabel(cp, fmt='%d',fontsize=16, colors='black')
axs[0].axis([cv1min,cv1max,levels[0],levels[-1]])
axs[0].get_xaxis().set_visible(False)
axs[0].set_ylabel(r'$\mathbf{\Delta G\,\, (kJ/mol)}$', fontsize=20)
axs[2].axis([cv1min,cv1max,cv2min,cv2max])
axs[2].set_xlabel(r'$\mathbf{\alpha-helix\,\, (Residues)}$', fontsize=20)
axs[2].set_ylabel(r'$\mathbf{ rgyr\,\, (nm)}$', fontsize=20)
axs[2].set_ylim(cv2min,cv2max)
axs[0].tick_params(axis='both', which='major', labelsize=18)
axs[1].tick_params(axis='both', which='major', labelsize=18)
axs[2].tick_params(axis='both', which='major', labelsize=18)
axs[3].tick_params(axis='both', which='major', labelsize=18)

#Limites da gráfica
axs[2].set_xlim(cv1min,cv1max)
axs[0].set_xlim(cv1min,cv1max)
axs[3].axis([levels[0],levels[-1],cv2min,cv2max])
axs[3].get_yaxis().set_visible(False)
axs[3].set_xlabel(r'$\mathbf{\Delta G\,\, (kJ/mol)}$', fontsize=20)
axs[3].set_ylim(cv2min,cv2max)



##################################
#-----INCERTIDUMBRE----------#####
##################################
a=np.loadtxt('CV1_BA.dat')
energia=a[:,1]
pos=a[:,0]
inc=a[:,2]
axs[0].plot(pos,yd1[0],color='green')
y_incerteza=[]
for u in range(len(pos)):
    if abs(pos[u])<10 and abs(pos[u])>0.2:
        y_incerteza.append(yd1[0][u])
plato=np.mean(y_incerteza);plato_inc=np.std(y_incerteza)
print('PLATO: ',plato,'+-',plato_inc)
print('POSSS',pos)
axs[0].fill_between(pos, np.array(yd1[0])-np.array(inc),np.array(yd1[0])+np.array(inc),alpha=0.3,color='green')
a=np.loadtxt('CV2_BA.dat')
energia=a[:,1]
pos=a[:,0]
inc=a[:,2]
angulo_promedio_inc=np.mean(inc)
axs[3].plot(yd2[0],(pos),color='green')
axs[3].fill_betweenx((np.array(pos)), np.array(yd2[0])-np.array(inc),np.array(yd2[0])+np.array(inc),alpha=0.2,color='green')
axs[2].set_xlim(cv1min,cv1max)
axs[0].set_xlim(cv1min,cv1max)
from matplotlib.ticker import MaxNLocator
axs[0].yaxis.set_major_locator(MaxNLocator(5))
axs[0].set_ylim(0,pmfmax)
axs[3].set_xlim(0,pmfmax)
#Generamos linea posición fosforos
axs[0].plot
fig.tight_layout()
#plt.tight_layout(pad = 1, w_pad = -1.2, h_pad = -0.4)
plt.savefig("fes.png")
plt.show()
plt.close()
#3D#
import math as mt
import pyvista as pv
import matplotlib.ticker as mticker
from matplotlib.ticker import FuncFormatter
from scipy.spatial import ConvexHull, convex_hull_plot_2d
def volumen_3D(X,Y,Z,YD1,correccion=1.,volumen=False,grafica=False,grafica_convex=False,limite=3.67):
    print('Correccion: ',correccion)
    if grafica==True:
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.contour3D(X, Y, Z, 50,cmap=cm.rainbow.reversed())
        f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
        g = lambda x,pos : "${}$".format(f._formatSciNotation('%1.10e' % x))
        plt.gca().zaxis.set_major_formatter(mticker.FuncFormatter(g))
        cset=ax.contourf(X, Y, Z, zdir='z',offset=0,cmap=cm.rainbow.reversed())
        cset = ax.contourf(X, Y, Z, xdir='x',offset=np.max(X),cmap=cm.rainbow.reversed())
        cset = ax.contourf(X, Y, Z, ydir='y',offset=np.max(Y),cmap=cm.rainbow.reversed())
        ax.grid(False)
        #ax.set_xlim(-2.2,-1.3)
        #ax.set_ylim(-0.3,0.3)
        #ax.set_zlim(0,0.95)
        ax.set_xlabel(r'$\mathbf{ Distance\,\, (nm)}$', fontsize=22,labelpad=30)
        ax.set_ylabel(r'$\mathbf{ Cosine}$', fontsize=22,labelpad=30)
        #ax.set_zlabel(r'$\mathbf{\Delta G\,\, (kJ/mol)}$', fontsize=20,labelpad=20)
        ax.set_zlabel(r'$\mathbf{Probability }$',fontsize=22,labelpad=30)
        #ax.set_zlim(0,0.95)
        ax.tick_params(axis='both', which='major', labelsize=22)
        plt.show()
        plt.close()
    if volumen==True:
        #ampliacion start
        Zmesh=np.copy(Z)
        f,c=np.shape(Z)
        Zmesh[np.isnan(Zmesh)] = 0
        Xmesh,Ymesh=np.meshgrid(X,Y)
        for i in range(f):
            for j in range(c):
                if Z[i][j]>max(YD1):
                    Zmesh[i][j]=max(YD1)
                elif Z[i][j]<min(YD1):
                    Zmesh[i][j]=min(YD1)
                else:
                    Zmesh[i][j]=Z[i][j]
        ############################################VOL4 START
        delta_x=[];delta_y=[]
        for i in range(len(X)-1):
        	delta_x.append(abs(X[i+1]-X[i]))
        	delta_y.append(abs(Y[i+1]-Y[i]))
        matrix=np.zeros((len(delta_x),len(delta_y)))
        matrix_z=np.zeros((len(delta_x),len(delta_y)))
        for i in range(len(delta_x)):
            for j in range(len(delta_y)):
                matrix[i][j]=delta_x[i]*delta_y[j]
                if j==0:
                    matrix_z[i][j]=abs(Z[i+1][j]+Z[i][j])/2
                else:
                    matrix_z[i][j]=abs(Z[i][j+1]+Z[i][j])/2
        intermedio=np.multiply(matrix,matrix_z)
        intermedio[np.isnan(intermedio)] = 0
        volumen_total=np.sum(intermedio)
        print('VOLUMEN_TOTAL: ',volumen_total)
        print('-----------------------------------------------')
        print('---------CÁLCULO VOLUMEN OCUPADO/LIBRE---------')
        print('-----------------------------------------------')
        delta_x_libre=[];delta_y_libre=[]
        for i in range(len(X)-1):
            if abs(X[i])>=abs(limite):
                delta_x_libre.append(abs(X[i+1]-X[i]))
                delta_y_libre.append(abs(Y[i+1]-Y[i]))
        matrix_libre=np.zeros((len(delta_x_libre),len(delta_y_libre)))
        matrix_z_libre=np.zeros((len(delta_x_libre),len(delta_y_libre)))
        for i in range(len(delta_x_libre)):
            for j in range(len(delta_y_libre)):
                matrix_libre[i][j]=delta_x_libre[i]*delta_y_libre[j]
                if j==0:
                    matrix_z_libre[i][j]=abs(Z[i+1][j]+Z[i][j])/2
                else:
                    matrix_z_libre[i][j]=abs(Z[i][j+1]+Z[i][j])/2
        intermedio_libre=np.multiply(matrix_libre,matrix_z_libre)
        intermedio_libre[np.isnan(intermedio_libre)]=0
        volumen_libre=np.sum(intermedio_libre)
        volumen_ocupado=volumen_total-volumen_libre
        print('Volumen bonding: ',volumen_ocupado,'Volumen libre: ',volumen_libre)
        print('-----------------------------------------------')
        print('-------------NUEVO CÁLCULO DE LA K-------------')
        print('-----------------------------------------------')
       # print('K: ',volumen_ocupado*(8-limite)/(limite*volumen_libre))
       # print('K2: ',volumen_total/(2*1660**(1/3.)*np.exp(-globalmin_p/(R*T))/correccion))
       # print('CORRECCION =? Volumen libre',np.exp(-globalmin_p/(R*T))/correccion)
       # print('OKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK: ',volumen_ocupado)
        ############################################VOL4 END
        #grid = pv.StructuredGrid(Xmesh, Ymesh, Zmesh);grid.plot();mesh=pv.StructuredGrid(Xmesh, Ymesh, Zmesh);mesh.plot(show_edges=True, show_grid=True,cmap='turbo')
       # print('VOLUMEN_FORMA_2',grid.volume)
        '''Zmesh[np.isnan(Zmesh)] = 0
        intermediate = np.trapz(Zmesh, X)
        intermediate[np.isnan(intermediate)] = 0
        result = np.trapz(intermediate, Y)
        print('RESULTADO: ',result)
        with open('3D.dat','w') as files:
            for i in range (len(X)):
                for j in range(len(Y)):
                    if mt.isnan(Z[i][j])==True or Z[i][j]>max(YD1):
                        pass
                    else:
                        files.write(str(X[i])+'   '+ str(Y[j])+'   '+ str(Z[i][j])+"\n")
        if grafica_convex==True:
               archivo='3D.dat'#cargamos el archivo
               lectura=np.loadtxt(archivo,comments='%')#borramos comentarios
               hull = ConvexHull(lectura)
               fig = plt.figure()#figsize = (20,10),facecolor="w")
               ax = plt.axes(projection="3d")
               for simplex in hull.simplices:	
                      ax.plot3D(lectura[simplex, 0], lectura[simplex, 1],lectura[simplex, 2], 's-')
               plt.show()
               plt.close()
               print('Area: ',hull.area)
               print('Volumen: ',hull.volume)'''
        print('COMPARACION: ',volumen_total,'y ', volumen_ocupado)
        return volumen_ocupado#result#grid.volume
def transformada(YD1,T=298,R=8.314472/1000.):
    x=df[df.columns[0]].to_numpy()*cambio_de_ejes#Este -1, depende de la definición inicial, igual en tu caso se puede borrar
    y=df[df.columns[1]].to_numpy()
    z=np.exp(-(df[df.columns[2]].to_numpy()-globalmin)/(R*T))
    z=np.maximum(z,np.exp(-(plato-globalmin)/(R*T)));z[np.isnan(z)] = np.exp(-(plato-globalmin)/(R*T))
    normalizar_z=np.sum(z)
    z=z/normalizar_z
    print('ARGUMENTO: ',-(df[df.columns[2]].to_numpy()-globalmin)/(R*T))
    #print('MAX',np.exp(-max(YD1[-1])/(R*T)),'MIN: ',np.exp(-min(YD1[-1])/(R*T)))
    #Creamos un gradiente para a representación gráfica
    xi = np.linspace(np.min(x), np.max(x), 1000)
    yi = np.linspace(np.min(y), np.max(y), 1000)
    zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method='linear')
    zi[np.isnan(zi)] =0
    zi=np.maximum(zi,np.exp(-(plato-globalmin)/(R*T))/normalizar_z)
    correccion=np.sum(np.exp(-(df[df.columns[2]].to_numpy()-globalmin)/(R*T)))
    return zi,np.exp(-(YD1[-1]-globalmin)/(R*T))/np.sum(np.exp(-(YD1[-1]-globalmin)/(R*T))),normalizar_z#-

volumen_3D(xi,yi,zi,yd1[-1],1,False,True,False)

'''



#------------------------------------------------------------#
#-----------------------COLVAR_2D----------------------------#
#------------------------------------------------------------#

# Read 2D HILLS file
df=pd.DataFrame(pd.read_csv('HILLS2.dat',delim_whitespace=True, header=None, skiprows=9))

x=df[df.columns[1]].to_numpy()
y=df[df.columns[2]].to_numpy()
z=df[df.columns[0]].to_numpy()
import numpy as np
import matplotlib.pyplot as plt
import random
  
# Creating dataset and Creating bins
x_min = np.min(x)
x_max = np.max(x)
y_min = np.min(y)
y_max = np.max(y)
x_bins = np.linspace(x_min, x_max, 500)
y_bins = np.linspace(y_min, y_max, 200)
#Representación gráfica do HILLS file
fig, ax = plt.subplots(figsize =(10, 7))
from scipy.stats import gaussian_kde
import matplotlib.gridspec as gridspec
# fit an array of size [Ndim, Nsamples]
data = np.vstack([x, y])
kde = gaussian_kde(data)
# evaluate on a regular grid
xgrid = np.linspace(-8, 0,50)
ygrid = np.linspace(1, -1,50)
Xgrid, Ygrid = np.meshgrid(xgrid, ygrid)
Z = kde.evaluate(np.vstack([Xgrid.ravel(), Ygrid.ravel()]))
fig = plt.figure()
gs = gridspec.GridSpec(2, 2,hspace=0.1, wspace=0.1)
ax0 = plt.subplot(gs[1, 0])
axx = plt.subplot(gs[0, 0],sharex=ax)
axy = plt.subplot(gs[1, 1],sharey=ax)
nullfmt = plt.NullFormatter()
axx.xaxis.set_major_formatter(nullfmt)
axy.yaxis.set_major_formatter(nullfmt)
axx.tick_params(axis ='x', bottom=False)
axy.tick_params(axis ='y', left=False)




axx.tick_params(axis ='y',labelsize=18)
axy.tick_params(axis ='x', labelsize=18)
ax0.tick_params(axis ='both', labelsize=18)
ax0.set_yticks(np.arange(10, 170, 40))
ax0.set_ylim(10,160)
axx.set_ylabel('Probability',fontsize=22)
axy.set_xlabel('Probability',fontsize=22)
ax0.set_xlabel('Position (nm)',fontsize=22)
ax0.set_ylabel(r'Cosine',fontsize=22)
# Plot the result as an image
y=np.array(y)
y_min = np.min(y)
y_max = np.max(y)
x=np.array(x)*cambio_de_ejes# O -1 depende da tua definición inicial
x_min=np.min(x)
x_max=np.max(x)
x_bins = np.linspace(x_min, x_max, 200)
y_bins = np.linspace(y_min, y_max, 200)
c=ax0.hist2d(x,y,bins =[x_bins, y_bins], cmap = plt.cm.rainbow,density=True,vmin=0.0001,vmax=0.25)#norm=LogNorm()
axx.hist(x, density = True,bins=100)
axy.hist(y, density = True,orientation=u'horizontal',bins=180)
cb = fig.colorbar(c[3],ax=axy,spacing='uniform',ticks=[0,0.01,0.05,0.1,0.15,0.2,0.25,0.3])#ticks=[0.0001,0.001,0.01,0.1,1])
cb.ax.tick_params(labelsize=18)
cb.set_label("density",fontsize=22)
# Top plot
axx.set_xlim(ax0.get_xlim())
# Right plot
axy.set_ylim(ax0.get_ylim())
#Establecemos límites nas gráficas
ax0.set_yticks(np.arange(-1, 1., 0.5))
ax0.set_ylim(-1,1)
ax0.set_xlim(-8,0)

plt.tight_layout()
plt.show()
plt.close()'''

