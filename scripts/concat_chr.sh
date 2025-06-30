#!/bin/bash 

sw=({2..23}) 
rv=(1 2 3 5 6 7 8 10  11  12  13  14  15  16  17  18  19  20  21  22  23  24) 
rv_asm=( RVCU0 RVND1 RVND2 RVUO0 RVAZ1 RVAZ2 RVNR1 RVNR2 ) 
sw_asm=( SWCU0 SWPC0 SWWA0 ) 

for ((idx=0; idx<${#sw[@]}; idx++)); do 
    rv_out="" 
    sw_out="" 
    for ((i=0; i<${#rv_asm[@]}; i++)); do 
        rv_out+=" ${rv_asm[i]}${rv[idx]}_ungapped.fa" 
    done 

    for ((i=0; i<${#sw_asm[@]}; i++)); do 
    	  sw_out+=" ${sw_asm[i]}${sw[idx]}_ungapped.fa" 
    done 

  seqs=$(echo "$rv_out $sw_out") 
  echo "$seqs > chr_"${sw[idx]}".fa" 
  cat $seqs > chr_"${sw[idx]}".fa 
done 