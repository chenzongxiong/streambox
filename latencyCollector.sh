#!/bin/bash

i=0
N=100

LATENCY_FILE=./results/latency/ysb/latency.txt
test -f "$LATENCY_FILE" && mv "$LATENCY_FILE" "$LATENCY_FILE.bk"

while (( i < $N ))
do
    latency=`cat ./results/latency/ysb/$i.txt |grep "Latency:" | awk '{print $2}'`
    echo $latency >> "$LATENCY_FILE"
    ((i++))
done

cat "$LATENCY_FILE"
