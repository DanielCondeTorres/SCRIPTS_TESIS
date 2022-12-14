#!/bin/bash
# añade grupos al index.ndx
# r 1-7 & a O4 \n 1 | 7 \n name 7 O4CDA \n name 8 LIG_O4CDA \n q

#module load gcc/system openmpi/4.0.5 gromacs/2021-PLUMED-2.7.1
module load gcc/system openmpi/4.0.5_ft3_cuda gromacs/2021.4-plumed-2.8.0
# Los valores de cv1 y cv2 salen del PMF calculados con fes
#En este caso:
#$1 = cv1
#$2 = cv2
#$3 =energia
cv1=`awk 'NR == 1 {line = $1; min = $3}
          NR > 1 && $1 > 1.2 && $3 < min {line = $1; min = $3}
          END{print line}' fesd1d2.dat `
cv2=`awk 'NR == 1 {line = $2; min = $3}
          NR > 1 && $2 > 1.2 && $3 < min {line = $2; min = $3}
          END{print line}' fesd1d2.dat `
echo $cv1
echo $cv2 


declare -i n
n=0

echo
echo
#Esto debo modificarlo para que busque en todas las cvs, e igual que busque más minimos:
#En mi caso en el COLVAR
# CV1 Se encuentra en la columna 3 ($3)
# CV2 se encuentra en la columna 4 ($4)
    for i in {1..100}
    do
        times=`awk -v cv1=$cv1 -v cv2=$cv2 -v margin=$i '{if (($3>cv1-cv1*margin/100) && ($3<cv1+cv1*margin/100) && ($4>cv1-cv1*margin/100) && ($4<cv1+cv1*margin/100) && ($1%25==0)) printf "%f ", $1}' ../w_1/COLVAR* `   # Solo coge los tiempos que al dividirlos por 25 el resto es 0
        nfields=`echo $times|awk -F';' '{print NF}'`
        if [ $nfields -gt "0" ]
        then
            echo "Number of structures= " $nfields > out.log
            echo "Times= " $times >> out.log
            echo "Percentual deviation= " $i >> out.log
            break
        fi
    done

    echo
    echo


    for time in $times
    do
        gmx trjcat -f ../w_1/traj_comp*.xtc -o $time.xtc -b $time -e $time
        #echo -e "0 \n0 \n0" | gmx trjconv -s ../prod_1.tpr -f $time.xtc -o $time.xtc -center -pbc cluster -n ../index.ndx 
        echo -e "0 \n0" | gmx trjconv -s ../prod_1.tpr -f $time.xtc -o $time.pdb -fit rot+trans -n ../index.ndx
        rm  \#*
        rm *xtc
        n=$n+1
        if [ $n -gt 20 ]; then
            break
        fi
    done
