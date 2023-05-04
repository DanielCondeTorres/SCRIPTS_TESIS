import numpy as np
import re
def modificamos_archivo(aminoacido,diccionario_beads,diccionario_carga,archivo_escritura='Protein_A_modificado.itp',pH=7):
    'Script que se encarga de modificar el .itp de interés en función del pH'
    # Formato de cada línea de átomo
    atom_line_format = "{0:>5} {1:>5} {2:>5} {3:>5} {4:>5} {5:>5} {6:>7.4f}; {7}\n"
    #atom_line_format = "{0:>6} {1:<5} {2:>5} {3:<5} {4:<5} {5:>5} {6:>10.4f} ; {7}\n"
    # Leer el archivo y almacenar su contenido en una lista de líneas
    with open(str(archivo_escritura), 'r') as f:
        lines = f.readlines()
    # Buscar la sección [atoms]
    for i, line in enumerate(lines):
        if line.startswith('[ atoms ]'):
            start = i
            break
    # Buscar el final de la sección [atoms]
    for i in range(start+1, len(lines)):
        if lines[i].startswith('['):
            end = i
            print('Numero de beads: ',end-start-2)#restamos 2 porque estamos cogiendo lineas de más
            aminoacido_inicial=lines[start+1].split()
            aminoacido_final=lines[end-2].split()
            print('Aminoacido inicial: ',aminoacido_inicial[3], 'aminoacido final: ',aminoacido_final[3])
            break
    # Reemplazar el residuo requerido por requerido0 en la sección [atoms]
    for i in range(start+1, end):
        # Modificamos extremos cargados de nuestro peptido
        if str(aminoacido_final[3]) in lines[i] and pH==0 and aminoacido_final[2]==lines[i].split()[2] and lines[i].split()[4]=='BB':#a) Por defecto los extremo van a estar cargados, por lo que a pH 0 debemos modificar el extremo Cterminal para que se encuentre neutro b) La carga la va a tener la bead del BB
            
            fields = lines[i].split();atom_index = int(fields[0]);residue_index = int(fields[2]);residue_name = fields[3];atom_name = fields[4]
            #Modificación deseada
            atom_type = 'Na' #Bead neutra del extremo Cterminal
            charge = float(0)
            #Guardamos los datos
            line_number = i - start - 1

            atom_line = atom_line_format.format(atom_index, atom_type, residue_index, residue_name, atom_name, residue_index, charge, line_number)
            lines[i] = atom_line
 
        if str(aminoacido_inicial[3]) in lines[i] and pH==14 and aminoacido_inicial[2]==lines[i].split()[2] and lines[i].split()[4]=='BB':
            fields = lines[i].split();atom_index = int(fields[0]);residue_index = int(fields[2]);residue_name = fields[3];atom_name = fields[4]
            #Modificación deseada
            atom_type = 'Na' #Bead neutra del extremo Cterminal
            charge = float(0)
            #Guardamos los datos
            line_number = i - start - 1
            
            atom_line = atom_line_format.format(atom_index, atom_type, residue_index, residue_name, atom_name, residue_index, charge, line_number)
            lines[i] = atom_line
        try:
            aminoacido_actual=lines[i].split()
            aminoacido_siguiente=lines[i+1].split()
        except IndexError:
            continue
      # Modificamos cadenas laterales de nuestro péptido
        if str(aminoacido) in lines[i]:
            if str(aminoacido)=='HIS':
                lines[i] = re.sub(str(aminoacido), str(aminoacido)+'H', lines[i])
            else:
                lines[i] = re.sub(str(aminoacido), str(aminoacido)+'0', lines[i])
            try:
                # Crear campos de la línea de átomo
                fields = lines[i].split()
                atom_index = int(fields[0])
                atom_type = fields[1]
                residue_index = int(fields[2])

                residue_name = fields[3]
                atom_name = fields[4]
                charge = float(fields[6])
                #Aplicamos el cambio en la última bead sólo, que es la que nos interesa, es decir cuando se cambia de aa
                if atom_type==diccionario_beads[str(aminoacido)][0] and aminoacido_actual[3]!=aminoacido_siguiente[3]:
                    atom_type=diccionario_beads[str(aminoacido)][1]
                    charge=float(diccionario_carga[str(aminoacido)][1])
                    print(charge)
                # Crear línea de átomo con el formato deseado
                line_number = i - start - 1
                atom_line = atom_line_format.format(atom_index, atom_type, residue_index, residue_name, atom_name, residue_index, charge, line_number)
                lines[i] = atom_line
            except IndexError:
                continue
    # Escribir el archivo modificado a un nuevo archivo llamado "topol_modificado.top"
    with open('Protein_A_modificado.itp', 'w') as f:
        f.writelines(lines)
    return 


############Diccionario para intercambiar los tipos de BEADS de interés#################################
diccionario_beads_aminoacidos = {'ASP': ['Qa', 'P3'],'GLU':['Qa','P1'],'ARG':['Qd','P4'],'HIS':['SP1','SQd'],'LYS':['Qd','P1']}
diccionario_cargas={'ASP': ['-1', '0'],'GLU':['-1','0'],'ARG':['1','0'],'HIS':['0','1'],'LYS':['1','0']}
#Si el pH es ácido, el péptido va a ganar H+, es decir tendra carga positova
#pI a tener en cuenta
#Aspartica: neutro  <----3.90-----> carga -
#Glutámica: neutro  <----4.25 ----> carga -
#Arginina:  carga + <----12.48----> neutro
#Cisteina:  neutro  <----8.3------> carga - #solo tienen el neutro
#Histidina: carga + <----6.0------> neutro #por defecto pone la neutra
#Lisina:    carga + <----10.0-----> neutro
#Tirosina:  neutro <----10.1------> carga - #solo tienen el neutro
#Treonina:  neutro  <----13.60-----> carga - #solo tienen el neutro
#Serina:    neutro  <----13.60------> carga - #solo tienen el neutro
def beads_segun_pH(pH):
    #Por defecto las cadenas laterales las va a dar cargadas
    if pH==0:
        modificamos_archivo('ASP',diccionario_beads_aminoacidos,diccionario_cargas,archivo_inicial,pH)#Solo lo hacemos una vez
        modificamos_archivo('GLU',diccionario_beads_aminoacidos,diccionario_cargas)
        modificamos_archivo('HIS',diccionario_beads_aminoacidos,diccionario_cargas)
        #completamente positivo; cambiar beads extremos;
        #Seria N terminal: Qd, Cterminal: Na
    elif pH < 3.90:
        modificamos_archivo('ASP',diccionario_beads_aminoacidos,diccionario_cargas,archivo_inicial)
        modificamos_archivo('GLU',diccionario_beads_aminoacidos,diccionario_cargas)
        modificamos_archivo('HIS',diccionario_beads_aminoacidos,diccionario_cargas)
    elif pH < 4.25:
        modificamos_archivo('GLU',diccionario_beads_aminoacidos,diccionario_cargas,archivo_inicial)
        modificamos_archivo('HIS',diccionario_beads_aminoacidos,diccionario_cargas)
    elif pH < 6:
        modificamos_archivo('HIS',diccionario_beads_aminoacidos,diccionario_cargas,archivo_inicial) #aqui no modificamos nada, todos cargados
    elif pH < 10.00:
        pass
        #modificamos_archivo('HIS')#la histidina por defecto la pone neutra en este script cambiamos la carga por defecto, entonces aqui la estaríamos cargando, lo cual no tiene sentido
    elif pH <12.48:
        #modificamos_archivo('HIS')
        modificamos_archivo('LYS',diccionario_beads_aminoacidos,diccionario_cargas,archivo_inicial)
    elif pH==14:
        #modificamos_archivo('HIS')
        modificamos_archivo('LYS',diccionario_beads_aminoacidos,diccionario_cargas,archivo_inicial,pH)#Solo lo hacemos una vez
        #completamente negativo; cambiar beads extremos
        #Seria N terminal: Nd, Cterminal: Qa
    return
print('Valores de interés: 0,2,4,5,7,11,13 y 14')
archivo_inicial='Protein_A.itp'
ph_del_medio = float(input("pH del medio: "))
beads_segun_pH(ph_del_medio)
