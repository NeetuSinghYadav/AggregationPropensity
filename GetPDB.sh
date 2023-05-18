#!/bin/bash

filename1='Tuttle_pdbseq.csv'
filename2='Tuttle_pdbseq_NoHyphen.csv'

while IFS= read -r line1 && read -r line2  <&3;
do
   cp leaprc_org.in leaprc.in
   seq="$line1"
   file="$line2"
   echo -e "$line\n"
   sed -i "s/SSEEQ/$file/g" leaprc.in
   sed -i "s/SEQ/$seq/g" leaprc.in 
   tleap -f leaprc.in
   echo $seq, $file
   
done  < $filename2 3< $filename1 
rm leaprc.in




#tr "-" " "< xx.csv > xxx.csv
