for tar in *.gz;  do tar xvzf $tar; done
egrep -v "\#" HILLS > HILLS2.dat
