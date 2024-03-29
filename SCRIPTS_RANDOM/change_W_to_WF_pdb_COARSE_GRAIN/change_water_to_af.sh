#/bin/bash
python change_W_to_AF_percentage.py -f equilibrado.pdb

grep -wv "WF" new.pdb > new2.pdb
grep -wv 'TER' new2.pdb > new22.pdb
grep -wv 'ENDMDL' new22.pdb > new222.pdb
grep -w "WF" new.pdb > new3.pdb

grep -w "TER" new.pdb > new4.pdb
grep -w "ENDMDL" new.pdb > new5.pdb
cat new222.pdb new3.pdb new4.pdb new5.pdb > output.pdb
rm new.pdb new?.pdb
