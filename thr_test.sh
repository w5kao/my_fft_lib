#!/bin/bash

cpupower frequency-set -g performance

for i in {0..4}; do
	./my_fft_test 12288 12288 5 $i 12k_thr_test.dat
done

for i in {5..8}; do
	./my_fft_test 12288 12288 10 $i 12k_thr_test.dat
done
