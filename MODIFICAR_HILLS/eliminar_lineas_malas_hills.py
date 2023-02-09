#Si el HILLS tiene alguna linea mal
new = open('HILLS_CORREGIDO2','w+')
with open( 'HILLS', 'r+' ) as file:
    for line in file:
        if len(line.split()) == 7:#numero lineas del HILLS
            try:
                # problema de ALEX
                if not round(float(line.split()[0]))==5 :#or float(line.split()[1])<0:#ver primer valor de hills, si es negativo así
                    new.write(line)
            except ValueError:
                pass
#new.close()
