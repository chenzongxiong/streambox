#!/bin/bash

echo "################################################################################"
echo "#             WARN: MAY STUCK AT SOME STEP. STOP AND RUN AGAIN                 #"
echo "################################################################################"
INPUT_FILE=/home/zxchen/inputdata.bin
CORES=8
i=0
N=100
while (( i < $N ))
do
    echo "i=$i"
    ./streambox/Release/test-yahoo.bin --input_file $INPUT_FILE --window_size 2 --records 12820 --cores $CORES &> ./results/latency/ysb/$i.txt
    (( $? == 0 )) && ((i++))
    sleep 3
done
