#bin/bash

squeue --format="%.18i %.9P %.17j %.12u %.8T %.10M %.10l %.6D %R" --me --sort=+j | awk '{if($5=="RUNNING"){printf "\033[32;1m" "%-12s %-18s %-25s %-15s %-15s %-15s %-15s \n" ,$1,$2,$3,$4,$5,$6,$7}else if($5=="PENDING"){printf "\033[31;1m"  "%-12s %-18s %-25s %-15s %-15s %-15s %-15s \n" ,$1,$2,$3,$4,$5,$6,$7}else if($5=="STATE"){printf "\033[33;1m"  "%-12s %-18s %-25s %-15s %-15s %-15s %-15s \n" ,$1,$2,$3,$4,$5,$6,$7}}'


echo " No more jobs bitch..." | awk '{print"\033[37m" $1" "$2" "$3 " "$4}'
                                                                             
