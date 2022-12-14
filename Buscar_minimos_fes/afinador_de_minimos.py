# abrimos el archivo en modo lectura
with open("fesd1d2.dat", "r") as f:
    # leemos todas las líneas del archivo y las almacenamos en una lista
    lineas = f.readlines()

    # asumimos que el valor mínimo de la primera columna es infinito
    minimo = float("inf")

    # iteramos sobre cada línea del archivo
    for linea in lineas:
        # dividimos cada línea en sus diferentes columnas y las convertimos a números
        col1, col2, col3 = map(float, linea.split())

        # si el valor de la primera columna en esta línea es menor que el valor mínimo actual,
        # actualizamos el valor mínimo y almacenamos los valores de esta línea
        if col3 < minimo:
            minimo = col3
            valores = (col1, col2, col3)

# imprimimos los valores de la fila que contiene el valor mínimo de la primera columna
print(valores)
