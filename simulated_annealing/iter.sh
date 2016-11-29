#!/bin/bash

rm output5*

for run in {1..10}
do
python spin_glass_SA.py input5.txt >> output5_SA.txt
done

sort -n output5_SA.txt -o output5_SA.txt 
head -n 1 output5_SA.txt 


