#!/bin/bash

echo "################################################################################"
echo "#             WARN: MAY STUCK AT SOME STEP. STOP AND RUN AGAIN                 #"
echo "################################################################################"

i=0
N=100
while (( i < $N ))
do
    echo "i=$i"
    ./streambox/Release/test-yahoo.bin --input_file /home/zxchen/inputdata.bin --window_size 2 --records 12820 --cores 8 &> ./results/latency/ysb/$i.txt
    (( $? == 0 )) && ((i++))
    sleep 3
done
