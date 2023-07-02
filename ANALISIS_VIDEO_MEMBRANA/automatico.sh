#!/bin/bash
#SBATCH -t 12:30:00 # execution time. Ex: 1 hour
#SBATCH -n 2 -c 4 # number of tasks, number of cores
#SBATCH --ntasks-per-node=2
#SBATCH --mem=64G



#CARGAMOS MODULOS NECESARIOS
module load gcc/system openmpi/4.0.5_ft3_cuda gromacs/2021.4-plumed-2.8.0
module load cesga/2020 mdanalysis/2.1.0

#Creamos la trayectoria a partir del archivo xtc original
#MEMBRANA Y PROTEINA 16, es el grupo que nos interesa, por lo que debemos crear el grupo proteina y lipidos
gmx trjconv -s ../prod.tpr -f ../traj_comp.part0001.xtc -o new.xtc -n ../membrana_indice.ndx -e 2225000 <<EOF
16
q
EOF 
gmx trjconv -s ../prod.tpr -f ../traj_comp.part0001.xtc -o new.pdb -n ../membrana_indice.ndx -b 0 -e  0 <<EOF
16
q
EOF
#Empleamos el escript de filtrado
./filter.sh new.xtc new.pdb 10
python creador_leaflet.py new.pdb "name PO4" l0 l1
gmx trjconv -f new.xtc -s new.pdb -n l0.ndx -o l0.xtc 
gmx trjconv -f new.xtc -s new.pdb -n l1.ndx -o l1.xtc

#Seleccion de la proteina
gmx trjconv -s new.pdb -f filter_10.xtc -o protein.xtc <<EOF
1
q
EOF
 
gmx trjconv -f new.xtc -s new.pdb -o protein.gro -b 0 -e 0 <<EOF
1
q
EOF


#Extraemos x y z de la membrana

python tananho_membrana.py -f new.pdb
output=$(python tananho_membrana.py -f new.pdb)
numbers=$(echo "$output" | awk '{print $2}')
zmax=$(echo "$output" | awk '/zmax/ {print $2}' | bc -l | xargs printf "%.0f")
zmin=$(echo "$output" | awk '/zmin/ {print $2}' | bc -l | xargs printf "%.0f")
xmax=$(echo "$output" | awk '/xmax/ {print $2}' | bc -l | xargs printf "%.0f")
xmin=$(echo "$output" | awk '/xmin/ {print $2}' | bc -l | xargs printf "%.0f")
ymax=$(echo "$output" | awk '/ymax/ {print $2}' | bc -l | xargs printf "%.0f")
ymin=$(echo "$output" | awk '/ymin/ {print $2}' | bc -l | xargs printf "%.0f")

# Utilizar las variables individuales como desees
echo "zmax: $zmax"
echo "zmin: $zmin"
echo "xmax: $xmax"
echo "xmin: $xmin"
echo "ymax: $ymax"
echo "ymin: $ymin"

################ejecutamos los scripts y los guardamos en sus carpetas

python profundidad.py protein.gro protein.xtc l0.gro l0.xtc "name PO4" 8  $xmin $xmax $ymin $ymax $zmin $zmax 0 6005 2 test .png true 4
mkdir imagenes_profundidad_l0
mv *png imagenes_profundidad_l0

python profundidad.py protein.gro protein.xtc l1.gro l1.xtc "name PO4" 8  $xmin $xmax $ymin $ymax $zmin $zmax 0 6005 2 test .png true 4
mkdir imagenes_profundidad_l1
mv *png imagenes_profundidad_l1

python difusion_lipidos.py protein.gro protein.xtc l0.gro l0.xtc l1.gro l1.xtc "name PO4" 8  $xmin $xmax $ymin $ymax  0 6005 2 test png true 4
mkdir imagenes_difusion
mv *png imagenes_difusion
#Entramos en cada carpeta
cp hacer_video.py imagenes_difusion
cp hacer_video.py imagenes_profundidad_l0
cp hacer_video.py imagenes_profundidad_l1

cd imagenes_difusion
python hacer_video.py
cd ..
cd imagenes_profundidad_l0
python hacer_video.py
cd ..
cd imagenes_profundidad_l1
python hacer_video.py

