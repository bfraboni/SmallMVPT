#!/bin/bash

make 

# compute reference image : 4h / scene for 9 images
# for i in $(seq 0 3) 
# do 
# 	./smallmvpt -a mvpt -s $i -t 432000; 
# done
# for ((exp=1; exp<129; exp*=2)) 
# do 
#   ./smallmvpt -s 4 -t 432000; 
# done

# compare algorithm with scenes from smallmvpt
echo "" > render.txt
for ((scene=0; scene<4; scene++)) 
do 
	for ((iter=1; iter<513; iter*=2)) 
	do 
		time1=$(./smallmvpt -r -s $scene -i $iter)  
		time2=$(./smallmvpt -s $scene -t $time1)
    	echo "$scene $iter $time1 $time2" >> render.txt
	done
done

# glossy scene with constant phong exponent
echo "" > glossy.txt
for ((exp=1; exp<=128; exp*=2)) 
do 
    time1=$(./smallmvpt -r -s 4 $exp -i 64)  
    time2=$(./smallmvpt -s 4 $exp -t $time1)
	echo "4 $exp $time1 $time2" >> glossy.txt
done
