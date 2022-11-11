##################################################################
#  Creacion de los archivos necesarios para obtener los PMFs     #
##################################################################
#bin/bash
module load gcc/system openmpi/4.0.5_ft3_cuda gromacs/2021.4-plumed-2.8.0
# compute distances
plumed sum_hills --hills HILLS --mintozero --outfile fesd1d2d2.dat
plumed sum_hills --hills HILLS --mintozero --idw angle --kt 2.5 --outfile fesd22_ --stride 10000
plumed sum_hills --hills HILLS --mintozero --idw D.z --kt 2.5 --outfile fesd11_ --stride 10000
egrep -v "\#" fesd1d2d2.dat | awk '{if (NF==5) print $1, $2, $3}' > fesd1d2.dat
lineas=`ls fesd11_* | wc -l`
echo $lineas
for ((i=0; i<$lineas; i=i+1))
do
        echo "jay fesd22_$i.dat"
        egrep -v "\#" fesd22_$i.dat > fesd2_$i.dat
        egrep -v "\#" fesd11_$i.dat > fesd1_$i.dat
done
rm fesd22_* fesd11*
exit;
