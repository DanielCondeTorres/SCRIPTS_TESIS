import cv2
import imageio
import os

# Ruta del directorio que contiene los archivos PNG
directorio = '.'

# Obtener la lista de archivos PNG en el directorio
archivos_png = [archivo for archivo in os.listdir(directorio) if archivo.endswith('.png')]

# Ordenar los archivos PNG en orden numérico
archivos_png.sort()

# Crear una lista para almacenar las imágenes
imagenes = []

# Iterar sobre los archivos PNG y agregar las imágenes a la lista
for archivo_png in archivos_png:
    ruta_archivo = os.path.join(directorio, archivo_png)
    imagen = imageio.imread(ruta_archivo)
    imagenes.append(imagen)

            # Guardar las imágenes como un archivo GIF
imageio.mimsave('animacion.gif', imagenes, duration=0.1)


