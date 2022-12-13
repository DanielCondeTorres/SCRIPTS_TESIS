#bin/bash

squeue --format="%.18i %.9P %.22j %.12u %.8T %.10M %.10l %.6D %R" --me --sort=+j | awk '{if($5=="RUNNING"){printf "\033[32;1m" "%-12s %-20s %-25s %-15s %-15s %-15s %-15s \n" ,$1,$2,$3,$4,$5,$6,$7}else if($5=="PENDING"){printf "\033[31;1m"  "%-12s %-20s %-25s %-15s %-15s %-15s %-15s%-45s\n" ,$1,$2,$3,$4,$5,$6,$7,$9}else if($5=="STATE"){printf "\033[33;1m"  "%-12s %-20s %-25s %-15s %-15s %-15s %-15s \n" ,$1,$2,$3,$4,$5,$6,$7}}'

echo -e " \033[5mNo more jobs bitch..." | awk '{print"\033[37m" $1" "$2" "$3 " "$4}'
echo -e "\033[0m"
