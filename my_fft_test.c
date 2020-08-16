/* This is SMP-parallel library for 2D DFTs */
/* v 0.1 Copyright (c) Summer 2019, Summer 2020 Alexander O. Korotkevich, w5kao@yahoo.com */
/***********************************************************************
 * This file is part of my_fft_lib.

    fft_lib is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License or Lesser General Public License
    as published by the Free Software Foundation, either version 3 of the License, or
    any later version.

    fft_lib is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with fft_lib.  If not, see <https://www.gnu.org/licenses/>. 
*************************************************************************/

/* Program for benchmarking of my_fft library. */

/*** N and M have to be multiple of BLOCK_SIDE_SIZE ***/
#include <stdio.h>
#include <stdlib.h>

#include <complex.h>
#include <math.h>

#include <fftw3.h>

#include <time.h>

#include "my_fft_lib.h"

#define VERBOSE
#define _SUPER_VERBOSE

#define ACCURACY_CONTROL

#define RUNS_NUMBER 10
#define _FORGET_WISDOM_BEFORE_FFTW
#define _READ_THR_WISDOM_FOR_FFTW
#define _SAVE_THR_WISDOM_FOR_FFTW
#define _FFTW_PATIENT_PLAN

fftw_plan *plans_fwd_X, *plans_fwd_Y;
double complex *in, *out, *out2;

unsigned long int NX=0, NY=0, num_reps=0, threads_number=0;

void complex_array_print (double complex *input, unsigned long int NX_size, unsigned long int NY_size);

double complex_array_control_sum (double complex *input, unsigned long int num_els) {
	unsigned long int i;
	double sum = 0.0;

	for (i = 0; i < num_els; ++i) {
		sum += cabs(input[i]);
	}

	return sum;
}

void complex_array_print (double complex *input, unsigned long int rows_N, unsigned long int cols_M) {
	unsigned long int i, j;
	double complex *temp_ptr;

	printf ("\n");
	temp_ptr = input;
	for (i = 0; i < rows_N; ++i) {
		for (j = 0; j < cols_M; ++j) {
			printf ("%.2e+i*%.2e ", creal(*temp_ptr), cimag(*temp_ptr));
			++temp_ptr;
		}
		printf ("\n");
	}
}

void complex_array_indexes (double complex *input, unsigned long int rows_N, unsigned long int cols_M) {
	unsigned long int i, j;
	double complex *temp_ptr;

	temp_ptr = input;
	for (i = 0; i < rows_N; ++i) {
		for (j = 0; j < cols_M; ++j) {
			(*temp_ptr) = (i+1) + I*(j+1);
			++temp_ptr;
		}
		printf ("\n");
	}
}

void complex_array_harmonics (double complex *input, unsigned long int NX_size, unsigned long int NY_size, double k_x, double k_y) {
	unsigned long int i, j;
	double complex *temp_ptr;
	double x, y, delta_x, delta_y;

	delta_x = 2.0*(M_PI)/NX_size;
	delta_y = 2.0*(M_PI)/NY_size;
	temp_ptr = input;
	for (i = 0; i < NX_size; ++i) {
		for (j = 0; j < NY_size; ++j) {
			x = delta_x*i;
			y = delta_y*j;
			(*temp_ptr) = cexp(I*(k_x*x + k_y*y));
			++temp_ptr;
		}
	}
}

void complex_array_element_print (double complex *input, unsigned long int rows_N, unsigned long int cols_M, unsigned long int row, unsigned long int col) {
	double complex element = input[row*cols_M + col];
	printf ("array[%lu,%lu] = %.15e + i*%.15e\n", row, col, creal(element), cimag(element));

}

int main (int argc, char** argv) {
	double mismatch=0.0, seconds, my_seconds_min, FFTW_seconds_min, ops, MFLOPS, my_MFLOPS_max, FFTW_MFLOPS_max, infty_norm, abs_mismatch;
	unsigned long int i, runs;
	FILE *f_p;
	fftw_plan normal_plan;
	my_fft_plan plan_fwd;
	struct timespec start, finish;

	/*Parsing command line.*/
	if (argc < 5)
	{
		printf ("Usage is the following:\n");
		printf("%s NX NY num_of_repeats number_of_threads [output_file]\n", argv[0]);
		exit (1);
	}
	
	sscanf(argv[1], "%lu", &NX);
	sscanf(argv[2], "%lu", &NY);
	sscanf(argv[3], "%lu", &num_reps);
	sscanf(argv[4], "%lu", &threads_number);
#ifdef VERBOSE
	printf ("NX = %lu, NY = %lu, num_reps = %lu, threads_number=%lu\n", NX, NY, num_reps, threads_number);
#endif /* VERBOSE */
/*** initialization and calculation of arrays ***/

	in = (double complex*) fftw_malloc ((NX)*(NY)*sizeof(double complex));
	out = (double complex*) fftw_malloc ((NX)*(NY)*sizeof(double complex));
	out2 = (double complex*) fftw_malloc ((NX)*(NY)*sizeof(double complex));
	my_fft_init (threads_number);
/*  Check that transpose works
	complex_array_indexes (in, NX, NY);
	complex_array_print (in, NX, NY);
	arrays_transpose (in, out, NX, NY);
	complex_array_print (out, NY, NX);
	exit(0);
*****************************/

	plan_fwd = my_fft_plan_dft_2d (in, out, out2, NX, NY, +1, FFTW_EXHAUSTIVE, threads_number);
	/*
	for (i = 0; i < NX*NY; ++i) {
		in[i] = sin(i) + I*cos(i);
	}
	*/
	my_MFLOPS_max = 0.0;
	for (runs=0; runs < (RUNS_NUMBER); ++runs) {
		complex_array_harmonics (in, NX, NY, -1.0, -2.0);

		clock_gettime(CLOCK_MONOTONIC, &start);
		for (i = 0; i < num_reps; ++i) {
			my_fft_execute (plan_fwd);
		}
		clock_gettime(CLOCK_MONOTONIC, &finish);

		seconds = (finish.tv_sec - start.tv_sec);
		seconds += (finish.tv_nsec - start.tv_nsec) / (1.0e+9);
		printf ("Time for %lu Fourier transforms using my_fft is %.11e seconds.\n", num_reps, seconds);
		ops = 5*NX*NY*log2(NX*NY);
		MFLOPS = (1.0e-6)*num_reps*ops/seconds;
		if (MFLOPS > my_MFLOPS_max) {
			my_MFLOPS_max = MFLOPS;
			my_seconds_min = seconds;
		}
	}
	printf ("MFLOPS = %lf\n", my_MFLOPS_max);
#ifdef SUPER_VERBOSE
	printf ("Control sum of cabs(out-array) is: %.15e\n", complex_array_control_sum(out, NX*NY));
	complex_array_element_print (out, NX, NY, 1, 2);
	complex_array_element_print (out, NX, NY, 2, 1);
#endif /* SUPER_VERBOSE */

	my_fft_destroy_plan (plan_fwd);
#ifdef FORGET_WISDOM_BEFORE_FFTW
	fftw_forget_wisdom();
#endif /* FORGET_WISDOM_BEFORE_FFTW */
	fflush(stdout);

#ifdef ACCURACY_CONTROL
	if (threads_number > 0) {
#ifdef READ_THR_WISDOM_FOR_FFTW
		if ((f_p = fopen("fftw.thr.wisdom", "rb")) == NULL) {
			printf ("Can't read wisdom file!\n");
			printf ("Trying to calculate plans without wisdom...\n");
		} else {
			if (fftw_import_wisdom_from_file (f_p)) {
				printf ("Wisdom successfully imported.\n");
			} else {
				printf ("File with wisdom is corrupted.\n");
				printf ("Trying to calculate plans without wisdom...\n");
			}
			fclose (f_p);
		}
#endif /* READ_THR_WISDOM_FOR_FFTW */
		fftw_init_threads();
		fftw_plan_with_nthreads(threads_number);

#ifdef FFTW_PATIENT_PLAN
		printf ("Using FFTW_PATIENT option for FFTW plan. Usually takes very long time.\n");
		normal_plan = fftw_plan_dft_2d(NX, NY, in, out2, +1, FFTW_PATIENT);
#else /* FFTW_PATIENT_PLAN */
		printf ("Using FFTW_MEASURE option for FFTW plan. Usually quite fast.\n");
		normal_plan = fftw_plan_dft_2d(NX, NY, in, out2, +1, FFTW_MEASURE);
#endif /* FFTW_PATIENT_PLAN */

#ifdef SAVE_THR_WISDOM_FOR_FFTW
		if ((f_p = fopen("fftw.thr.wisdom", "w+")) == NULL) {
			printf ("Can't open wisdom file for writing!\n");
			printf ("Wisdom has been lost!\n");
		} else {
			printf ("Exporting wisdom to file...\n");
			fftw_export_wisdom_to_file (f_p);
			fclose (f_p);
		}
#endif /* SAVE_THR_WISDOM_FOR_FFTW */
	} else {
		if ((f_p = fopen("fftw.wisdom", "rb")) == NULL) {
			printf ("Can't read wisdom file!\n");
			printf ("Trying to calculate plans without wisdom...\n");
		} else {
			if (fftw_import_wisdom_from_file (f_p)) {
				printf ("Wisdom successfully imported.\n");
			} else {
				printf ("File with wisdom is corrupted.\n");
				printf ("Trying to calculate plans without wisdom...\n");
			}
			fclose (f_p);
		}

		normal_plan = fftw_plan_dft_2d(NX, NY, in, out2, +1, FFTW_PATIENT);

		if ((f_p = fopen("fftw.wisdom", "w+")) == NULL) {
			printf ("Can't open wisdom file for writing!\n");
			printf ("Wisdom has been lost!\n");
		} else {
			printf ("Exporting wisdom to file...\n");
			fftw_export_wisdom_to_file (f_p);
			fclose (f_p);
		}
	}

	FFTW_MFLOPS_max = 0.0;
	for (runs=0; runs < (RUNS_NUMBER); ++runs) {
		complex_array_harmonics (in, NX, NY, -1.0, -2.0);

		clock_gettime(CLOCK_MONOTONIC, &start);
		for (i = 0; i < num_reps; ++i) {
			fftw_execute (normal_plan);
		}
		clock_gettime(CLOCK_MONOTONIC, &finish);

		seconds = (finish.tv_sec - start.tv_sec);
		seconds += (finish.tv_nsec - start.tv_nsec) / (1.0e+9);
		printf ("Time for %lu Fourier transforms using FFTW is %.11e seconds.\n", num_reps, seconds);
		ops = 5*NX*NY*log2(NX*NY);
		MFLOPS = (1.0e-6)*num_reps*ops/seconds;
		if (MFLOPS > FFTW_MFLOPS_max) {
			FFTW_MFLOPS_max = MFLOPS;
			FFTW_seconds_min = seconds;
		}
	}
	printf ("MFLOPS = %lf\n", FFTW_MFLOPS_max);
#ifdef SUPER_VERBOSE
	printf ("Control sum of cabs(out2-array) is: %.15e\n", complex_array_control_sum(out2, NX*NY));
	complex_array_element_print (out2, NX, NY, 1, 2);
	complex_array_element_print (out2, NX, NY, 2, 1);
#endif /* SUPER_VERBOSE */

	fftw_destroy_plan (normal_plan);
	fflush(stdout);

#ifdef VERBOSE
	infty_norm = 0.0;
	for (i = 0; i < NX*NY; ++i) {
		out[i] -= out2[i];
		abs_mismatch = cabs(out[i]);
		if (infty_norm < abs_mismatch) {
			infty_norm = abs_mismatch;
		}
	}
	printf ("Normed L_infty norm: ||Mismatch between my_fft and fftw||_infty/(NX*NY) = %.15e.\n", infty_norm/(NX*NY));
#endif /* VERBOSE */

	if (argc == 6) {
		if ((f_p = fopen(argv[5], "a+")) == NULL) {
			printf ("Can't open output file %s for writing!\n", argv[5]);
		} else {
			fprintf (f_p, "%lu\t%lu\t%lu\t%lf\t%lf\t%lf\t%lf\n", NX, NY, threads_number, my_MFLOPS_max, FFTW_MFLOPS_max, my_seconds_min, FFTW_seconds_min);
			fclose (f_p);
		}
	}
#endif /* ACCURACY_CONTROL */

	return 0;
}
