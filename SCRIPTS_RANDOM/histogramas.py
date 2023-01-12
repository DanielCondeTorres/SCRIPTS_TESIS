import matplotlib.pyplot as plt
import numpy as np
# Valores medios
etiquetas=['Mag2','Mag2 Res','Mag2 AF','Mag2 Pol','Mag2 Z', '1Z64 Z', '2JMY Z','2K60 Z','6C41 Z']
#El orden es: Mamifero, Cancer, Bacteria
Mag2 = [101,117,116]
Mag2_inc=[2,5,2]
Mag2_res= [100,117,118]
Mag2_res_inc=[4,4,5]
Mag2_af=[115.42,132.3,127.55]
Mag2_af_inc=[0.58,3.1,3.3]
Mag2_pol=[86.6,122.13,118.8]
Mag2_pol_inc=[1.1,2.1,2.5]
mag2_z=[138.9,122.51,135.6]
mag2_z_inc=[1.2,0.7,2.7]
z64=[133.5,133.7,145.7]
z64_inc=[1.9,2.4,3.0]
jmy=[76.88,101.60,121.44]
jmy_inc=[0.84,0.36,0.88]
k60=[151.0,191,183.2]
k60_inc=[8.2,23,7.6]
c41=[132.8,157.6,170.1]
c41_inc=[2.3,6.4,7.8]

colores=['Green','Orange','Blue']



# Crear el gráfico de barras
bar_width = 0.2
etiequetas_posiciones=[]
x_pos = np.arange(0,bar_width*3-bar_width,bar_width);print(x_pos)
plt.bar(x_pos, Mag2, yerr=Mag2_inc, width=bar_width,color=colores)
etiequetas_posiciones.append(x_pos[1])
plt.legend(Mag2,['Mammal', 'Cancer','Bacteria'])
x_pos = np.linspace(0.2*2+x_pos[-1],bar_width*3+x_pos[-1]+0.2,3)
plt.bar(x_pos, Mag2_res, yerr=Mag2_res_inc, width=bar_width,color=colores)
etiequetas_posiciones.append(x_pos[1])

x_pos = np.linspace(0.2*2+x_pos[-1],bar_width*3+x_pos[-1]+0.2,3)
plt.bar(x_pos, Mag2_af, yerr=Mag2_af_inc, width=bar_width,color=colores)
etiequetas_posiciones.append(x_pos[1])

x_pos = np.linspace(0.2*2+x_pos[-1],bar_width*3+x_pos[-1]+0.2,3)
plt.bar(x_pos, Mag2_pol, yerr=Mag2_pol_inc, width=bar_width,color=colores)
etiequetas_posiciones.append(x_pos[1])

x_pos = np.linspace(0.2*2+x_pos[-1],bar_width*3+x_pos[-1]+0.2,3)
plt.bar(x_pos, mag2_z, yerr=mag2_z_inc, width=bar_width,color=colores)
etiequetas_posiciones.append(x_pos[1])

x_pos = np.linspace(0.2*2+x_pos[-1],bar_width*3+x_pos[-1]+0.2,3)
plt.bar(x_pos, z64, yerr=z64_inc, width=bar_width,color=colores)
etiequetas_posiciones.append(x_pos[1])


x_pos = np.linspace(0.2*2+x_pos[-1],bar_width*3+x_pos[-1]+0.2,3)
plt.bar(x_pos, jmy, yerr=jmy_inc, width=bar_width,color=colores)
etiequetas_posiciones.append(x_pos[1])

x_pos = np.linspace(0.2*2+x_pos[-1],bar_width*3+x_pos[-1]+0.2,3)
plt.bar(x_pos, k60, yerr=k60_inc, width=bar_width,color=colores)
etiequetas_posiciones.append(x_pos[1])

x_pos = np.linspace(0.2*2+x_pos[-1],bar_width*3+x_pos[-1]+0.2,3)
plt.bar(x_pos, c41, yerr=c41_inc, width=bar_width,color=colores)
etiequetas_posiciones.append(x_pos[1])
# Añadir etiquetas al eje x
plt.xticks(etiequetas_posiciones, etiquetas,fontsize=18)

plt.ylabel(r'$\mathbf{\Delta G\,\, (kJ/mol)}$', fontsize=28)
plt.yticks(fontsize=28)


colors = {'Mammal':'Green', 'Cancer':'Orange', 'Bacteria':'Blue'}         
labels = list(colors.keys())
handles = [plt.Rectangle((0,0),1,1, color=colors[label]) for label in labels]
plt.legend(handles, labels, fontsize=22)
# Mostrar el gráfico
plt.show()
