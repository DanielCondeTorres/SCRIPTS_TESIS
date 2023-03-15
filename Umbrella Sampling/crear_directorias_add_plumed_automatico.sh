export LC_NUMERIC=en_US.UTF-8
for AT in $(seq 0.1 0.1 7.0)

do
cat >plumed.dat << EOF
RESTART
MOLINFO STRUCTURE=HELICE_GROMOS.pdb
#Potentials to avoid the membrane and/or PBC crossing
ALFA: ALPHARMSD RESIDUES=4-15 LESS_THAN={RATIONAL R_0=0.08 NN=8 MM=12} #RESIDUES=1-357 TYPE=OPTIMAL 
ALFA_RESTRAINT: RESTRAINT ARG=ALFA.lessthan KAPPA=1000 AT=$AT
PRINT ARG=ALFA.lessthan,ALFA_RESTRAINT.bias  STRIDE=1000  FILE=COLVAR #,metad.bias
EOF
result=$(echo "$AT * 10" | bc)
carpeta=$(printf "%.0f" $result)
echo $result
echo $carpeta
mkdir "$carpeta"
cp plumed.dat $carpeta
cp *tpr $carpeta
cp *pdb $carpeta
done
