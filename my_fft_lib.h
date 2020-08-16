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
#ifndef _MY_FFT_H
#define _MY_FFT_H

#include <stdio.h>
#include <stdlib.h>

#include <complex.h>
#include <math.h>

#include <fftw3.h>

#include "my_thr_lib.h"

#define _VERBOSE

#define _SUPER_VERBOSE

/* For Intel processors optimal block size is 64, for AMD Opteron - 16. */
#define BLOCK_SIDE_SIZE 64

double complex *scratch_array_my_fft;

typedef struct {
	fftw_plan *plans_X, *plans_Y;
	unsigned long int NX_size, NY_size, number_of_threads;
	fftw_complex *input, *output, *scratch_array;
} my_plan, *my_fft_plan;

void my_fft_init (unsigned long int number_of_threads);
my_fft_plan my_fft_plan_dft_2d (double complex *input, double complex *output, double complex *scratch_array, unsigned long int NX_size, unsigned long int NY_size, int DIR, unsigned FLAGS, unsigned long int number_of_threads);
void my_fft_destroy_plan (my_fft_plan plan);
void my_fft_execute (my_fft_plan plan);

void arrays_transpose (double complex *input, double complex *output, unsigned long int row_elements_N, unsigned long int column_elements_M);
void arrays_transpose_with_threads (double complex *input, double complex *output, unsigned long int column_elements_M, unsigned long int row_elements_N, unsigned long int number_of_threads);
void unit_transpose (double complex *input, double complex *output, unsigned long int current_block, unsigned long int row_elements_N, unsigned long int column_elements_M, unsigned long int row_blocks_number, unsigned long int column_blocks_number, unsigned long int row_block_skip, unsigned long int column_block_skip);

struct DFT_thread_input {
	fftw_plan *plan;
	unsigned long int elements_number;
};

void DFT_for_arrays_thr (void *input_thr);

struct arrays_transpose_thread_input {
	double complex *input, *output;
	unsigned long int index, elements_number, row_elements_N, column_elements_M, row_blocks_number, column_blocks_number, row_block_skip, column_block_skip;
};

void arrays_transpose_thr (void *input_thr);

struct plans_init_thread_input {
	unsigned long int index, elements_number;
};

void plans_init_thr (void *input_thr);

#endif
