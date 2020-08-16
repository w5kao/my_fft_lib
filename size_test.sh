#!/bin/bash

# If you don't have root rights, usually you cannot execute this command:
#cpupower frequency-set -g performance

for i in {4..12}; do
	size=$(echo $i*1024| bc)
	./my_fft_test $size $size 10 8 size_test.dat
done

