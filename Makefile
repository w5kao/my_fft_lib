
#CXXFLAGS=-O3 -Wall -W -Wcast-qual -Wpointer-arith -Wcast-align -pedantic -fno-schedule-insns -fschedule-insns2 -fstrict-aliasing -funroll-loops -fprefetch-loop-arrays -fomit-frame-pointer -I .
#LDFLAGS=-L ../../fftw3lib/lib -lfftw3 -lm -static
CXXFLAGS=-march=native -Ofast -flto -std=gnu99 -W -pedantic -Wcast-qual -Wpointer-arith -Wcast-align -fno-schedule-insns -fschedule-insns2 -fstrict-aliasing -funroll-loops -fprefetch-loop-arrays -fomit-frame-pointer -I ../../fftw3lib/include
#CXXFLAGS=-O0 -g -I ../../fftw3lib/include

LDFLAGS=-march=native -Ofast -flto -L ../../fftw3lib/lib -lfftw3_threads -lfftw3 -lm -lpthread -flto -static

cCC=/usr/bin/gcc
CXX=/usr/bin/gcc

all: my_fft_test

%.o: %.c
	$(cCC) -c $(CXXFLAGS) -o $@ $<

my_fft_test: my_fft_test.o my_fft_lib.o my_thr_lib.o
	$(cCC) -Wall -o $@ my_fft_test.o my_fft_lib.o my_thr_lib.o $(LDFLAGS)

clean:
	rm -f *.o
distclean: clean 
	rm my_fft_test; rm fftw.wisdom;
