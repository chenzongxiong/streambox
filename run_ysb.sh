#!/bin/bash

TIMEOUT=30                      # kill program after 30 seconds
echo "Measure thread 1"
timeout $TIMEOUT ./streambox/Release/test-yahoo.bin --records 1000000 --window_size 1 --input_file /home/zxchen/inputdata.bin --cores 1 &> throughputs/throughput-1.txt
echo "Measure thread 2"
timeout $TIMEOUT ./streambox/Release/test-yahoo.bin --records 1000000 --window_size 1 --input_file /home/zxchen/inputdata.bin --cores 2 &> throughputs/throughput-2.txt
echo "Measure thread 4"
timeout $TIMEOUT ./streambox/Release/test-yahoo.bin --records 1000000 --window_size 1 --input_file /home/zxchen/inputdata.bin --cores 4 &> throughputs/throughput-4.txt
echo "Measure thread 8"
timeout $TIMEOUT ./streambox/Release/test-yahoo.bin --records 1000000 --window_size 1 --input_file /home/zxchen/inputdata.bin --cores 8 &> throughputs/throughput-8.txt
echo "Measure thread 12"
timeout $TIMEOUT ./streambox/Release/test-yahoo.bin --records 1000000 --window_size 1 --input_file /home/zxchen/inputdata.bin --cores 12 &> throughputs/throughput-12.txt
echo "Measure thread 18"
timeout $TIMEOUT ./streambox/Release/test-yahoo.bin --records 1000000 --window_size 1 --input_file /home/zxchen/inputdata.bin --cores 18 &> throughputs/throughput-18.txt
echo "Measure thread 24"
timeout $TIMEOUT ./streambox/Release/test-yahoo.bin --records 1000000 --window_size 1 --input_file /home/zxchen/inputdata.bin --cores 24 &> throughputs/throughput-24.txt
echo "Measure thread 36"
timeout $TIMEOUT ./streambox/Release/test-yahoo.bin --records 1000000 --window_size 1 --input_file /home/zxchen/inputdata.bin --cores 36 &> throughputs/throughput-36.txt
echo "Measure thread 42"
timeout $TIMEOUT ./streambox/Release/test-yahoo.bin --records 1000000 --window_size 1 --input_file /home/zxchen/inputdata.bin --cores 42 &> throughputs/throughput-42.txt
