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
#include "my_fft_lib.h"

unsigned long int each_thread_elements, last_thread_elements;

struct DFT_thread_input *thr_array_DFT;
struct arrays_transpose_thread_input *thr_arrays_transpose;

void my_fft_init (unsigned long int number_of_threads) {
	if (number_of_threads > 0) {
		my_thr_pool_init (number_of_threads);
		thr_array_DFT = (struct DFT_thread_input *) malloc (number_of_threads*sizeof(struct DFT_thread_input));
		thr_arrays_transpose = (struct arrays_transpose_thread_input *) malloc (number_of_threads*sizeof(struct arrays_transpose_thread_input));
	}
}


/* FLAGS = FFTW_EXHAUSTIVE | FFTW_DESTROY_INPUT, DIR = +1 OR -1 */
my_fft_plan my_fft_plan_dft_2d (double complex *input, double complex *output, double complex *scratch_array, unsigned long int NX_size, unsigned long int NY_size, int DIR, unsigned FLAGS, unsigned long int number_of_threads) {
	FILE *f_p;
	unsigned long int i;
	double complex *temp_u, *temp_u_omega;
	fftw_plan *temp_plan;
	my_plan *return_plan;

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

	return_plan = (my_plan *) malloc (sizeof(my_plan));
	return_plan->NX_size = NX_size;
	return_plan->NY_size = NY_size;
	return_plan->number_of_threads = number_of_threads;
	return_plan->plans_Y = malloc(NX_size*sizeof(fftw_plan));
	return_plan->plans_X = malloc(NY_size*sizeof(fftw_plan));
	return_plan->input = input;
	return_plan->output = output;
	return_plan->scratch_array = scratch_array;

	temp_plan = return_plan->plans_Y;
	temp_u_omega = input;
	temp_u = scratch_array;
	for (i=0; i < NX_size; ++i) {
		(*temp_plan) = fftw_plan_dft_1d (NY_size, temp_u_omega, temp_u, DIR, FLAGS);

		++temp_plan;
		temp_u_omega += NY_size;
		temp_u += NY_size;
	}

	temp_plan = return_plan->plans_X;
	temp_u = output;
	temp_u_omega = scratch_array;
	for (i=0; i < NY_size; ++i) {
		(*temp_plan) = fftw_plan_dft_1d (NX_size, temp_u, temp_u_omega, DIR, FLAGS);

		++temp_plan;
		temp_u += NX_size;
		temp_u_omega += NX_size;
	}

	if ((f_p = fopen("fftw.wisdom", "w+")) == NULL) {
		printf ("Can't open wisdom file for writing!\n");
		printf ("Wisdom has been lost!\n");
	} else {
		printf ("Exporting wisdom to file...\n");
		fftw_export_wisdom_to_file (f_p);
		fclose (f_p);
	}

	return (my_fft_plan) return_plan;
}

void my_fft_destroy_plan (my_fft_plan plan) {
	unsigned long int i;
	fftw_plan *temp_plan;

	temp_plan = plan->plans_X;
	for (i=0; i < plan->NY_size; ++i) {
		fftw_destroy_plan((*temp_plan));

		++temp_plan;		
	}

	temp_plan = plan->plans_Y;
	for (i=0; i < plan->NX_size; ++i) {
		fftw_destroy_plan((*temp_plan));

		++temp_plan;		
	}
	free (plan);
}

void my_fft_execute (my_fft_plan plan) {
	unsigned long int i;
	fftw_plan *temp_plan;

	if (plan->number_of_threads == 0) {
		temp_plan = plan->plans_Y;
		for (i=0; i < plan->NX_size; ++i) {
			fftw_execute((*temp_plan));

			++temp_plan;		
		}

		arrays_transpose (plan->scratch_array, plan->output, plan->NX_size, plan->NY_size);

		temp_plan = plan->plans_X;
		for (i=0; i < plan->NY_size; ++i) {
			fftw_execute((*temp_plan));

			++temp_plan;		
		}

		arrays_transpose (plan->scratch_array, plan->output, plan->NY_size, plan->NX_size);
	} else {
		each_thread_elements = (unsigned long int) (0.5 + (1.0*(plan->NX_size))/(plan->number_of_threads));
		last_thread_elements = plan->NX_size - each_thread_elements*((plan->number_of_threads)-1);

		temp_plan = plan->plans_Y;
		for (i = 0; i < ((plan->number_of_threads)-1); ++i) {
			thr_array_DFT [i].plan = temp_plan;
			thr_array_DFT [i].elements_number = each_thread_elements;
			my_thr_data_assign (i,  (void *) &thr_array_DFT[i]);

			temp_plan += each_thread_elements;
		}

		thr_array_DFT [((plan->number_of_threads)-1)].plan = temp_plan;
		thr_array_DFT [((plan->number_of_threads)-1)].elements_number = last_thread_elements;

		my_thr_data_assign (((plan->number_of_threads)-1),  (void *) &thr_array_DFT[((plan->number_of_threads)-1)]);

		my_thr_manager (DFT_for_arrays_thr);

		arrays_transpose_with_threads (plan->scratch_array, plan->output, plan->NX_size, plan->NY_size, plan->number_of_threads);

		each_thread_elements = (unsigned long int) (0.5 + (1.0*(plan->NY_size))/(plan->number_of_threads));
		last_thread_elements = plan->NY_size - each_thread_elements*((plan->number_of_threads)-1);

		temp_plan = plan->plans_X;
		for (i = 0; i < ((plan->number_of_threads)-1); ++i) {
			thr_array_DFT [i].plan = temp_plan;
			thr_array_DFT [i].elements_number = each_thread_elements;
			my_thr_data_assign (i,  (void *) &thr_array_DFT[i]);

			temp_plan += each_thread_elements;
		}

		thr_array_DFT [((plan->number_of_threads)-1)].plan = temp_plan;
		thr_array_DFT [((plan->number_of_threads)-1)].elements_number = last_thread_elements;

		my_thr_data_assign (((plan->number_of_threads)-1),  (void *) &thr_array_DFT[((plan->number_of_threads)-1)]);

		my_thr_manager (DFT_for_arrays_thr);

		arrays_transpose_with_threads (plan->scratch_array, plan->output, plan->NY_size, plan->NX_size, plan->number_of_threads);
	}
}

void DFT_for_arrays_thr (void * input_thr) {
	struct DFT_thread_input *input;
	fftw_plan *temp_plan;
	unsigned long int i;

	input = (struct DFT_thread_input *) input_thr;

	temp_plan = input->plan;
	for (i = 0; i < input->elements_number; ++i) {
		fftw_execute(*temp_plan);

		++temp_plan;
	}
}

void arrays_transpose (double complex *input, double complex *output, unsigned long int column_elements_M, unsigned long int row_elements_N) {
	unsigned long int row_blocks_number, column_blocks_number, total_number_of_blocks, i,j, row_block_skip, column_block_skip;
	double complex *block_ptr_in, *block_ptr_out;

/*** We suppose that row_elements_N and column_elements_M are multiples of BLOCK_SIDE_SIZE, it can be generalized though. ***/
	row_blocks_number = row_elements_N/(BLOCK_SIDE_SIZE);
	column_blocks_number = column_elements_M/(BLOCK_SIDE_SIZE);

	total_number_of_blocks = row_blocks_number*column_blocks_number;
	row_block_skip = row_elements_N*(BLOCK_SIDE_SIZE);
	column_block_skip = column_elements_M*(BLOCK_SIDE_SIZE);

	for (i = 0; i < total_number_of_blocks; ++i) {
		unit_transpose (input, output, i, row_elements_N, column_elements_M, row_blocks_number, column_blocks_number, row_block_skip, column_block_skip);
	}
}

void arrays_transpose_with_threads (double complex *input, double complex *output, unsigned long int column_elements_M, unsigned long int row_elements_N, unsigned long int number_of_threads) {
	unsigned long int row_blocks_number, column_blocks_number, total_number_of_blocks, i,j, row_block_skip, column_block_skip;
	double complex *block_ptr_in, *block_ptr_out;

/*** We suppose that row_elements_N and column_elements_M are multiples of BLOCK_SIDE_SIZE, it can be generalized though. ***/
	row_blocks_number = row_elements_N/(BLOCK_SIDE_SIZE);
	column_blocks_number = column_elements_M/(BLOCK_SIDE_SIZE);

	total_number_of_blocks = row_blocks_number*column_blocks_number;
	row_block_skip = row_elements_N*(BLOCK_SIDE_SIZE);
	column_block_skip = column_elements_M*(BLOCK_SIDE_SIZE);

	each_thread_elements = (unsigned long int) (0.5 + (1.0*total_number_of_blocks)/(number_of_threads));
	last_thread_elements = total_number_of_blocks - each_thread_elements*(number_of_threads-1);

	for (i = 0; i < (number_of_threads-1); ++i) {
		thr_arrays_transpose [i].input = input;
		thr_arrays_transpose [i].output = output;
		thr_arrays_transpose [i].index = i*each_thread_elements;
		thr_arrays_transpose [i].elements_number = each_thread_elements;
		thr_arrays_transpose [i].row_elements_N = row_elements_N;
		thr_arrays_transpose [i].column_elements_M = column_elements_M;
		thr_arrays_transpose [i].row_blocks_number = row_blocks_number;
		thr_arrays_transpose [i].column_blocks_number = column_blocks_number;
		thr_arrays_transpose [i].row_block_skip = row_block_skip;
		thr_arrays_transpose [i].column_block_skip = column_block_skip;

		my_thr_data_assign (i,  (void *) &thr_arrays_transpose[i]);
	}

	thr_arrays_transpose [(number_of_threads-1)].input = input;
	thr_arrays_transpose [(number_of_threads-1)].output = output;
	thr_arrays_transpose [(number_of_threads-1)].index = (number_of_threads-1)*each_thread_elements;
	thr_arrays_transpose [(number_of_threads-1)].elements_number = last_thread_elements;
	thr_arrays_transpose [(number_of_threads-1)].row_elements_N = row_elements_N;
	thr_arrays_transpose [(number_of_threads-1)].column_elements_M = column_elements_M;
	thr_arrays_transpose [(number_of_threads-1)].row_blocks_number = row_blocks_number;
	thr_arrays_transpose [(number_of_threads-1)].column_blocks_number = column_blocks_number;
	thr_arrays_transpose [(number_of_threads-1)].row_block_skip = row_block_skip;
	thr_arrays_transpose [(number_of_threads-1)].column_block_skip = column_block_skip;

	my_thr_data_assign ((number_of_threads-1),  (void *) &thr_arrays_transpose[(number_of_threads-1)]);

	my_thr_manager (arrays_transpose_thr);
}

void arrays_transpose_thr (void * input_thr) {
	struct arrays_transpose_thread_input *input;
	unsigned long int i, final_index;

	input = (struct arrays_transpose_thread_input *) input_thr;

	final_index = input->index + input->elements_number;

	for (i = input->index; i < final_index; ++i) {
		unit_transpose (input->input, input->output, i, input->row_elements_N, input->column_elements_M, input->row_blocks_number, input->column_blocks_number, input->row_block_skip, input->column_block_skip);
	}
}

void unit_transpose (double complex *input, double complex *output, unsigned long int current_block, unsigned long int row_elements_N, unsigned long int column_elements_M, unsigned long int row_blocks_number, unsigned long int column_blocks_number, unsigned long int row_block_skip, unsigned long int column_block_skip) {
	double complex *block_input, *block_output, *temp_ptr_in, *temp_ptr_out;
	unsigned long int i, j;

	j = (unsigned long int) (current_block/row_blocks_number);
	i = current_block - (row_blocks_number*j);

	block_input = input + i*(BLOCK_SIDE_SIZE) + j*row_block_skip;
	block_output = output + i*column_block_skip + j*(BLOCK_SIDE_SIZE);

	for (j = 0; j < (BLOCK_SIDE_SIZE); ++j) {
		temp_ptr_in = block_input + row_elements_N*j;
		temp_ptr_out = block_output + j;
		for (i = 0; i < (BLOCK_SIDE_SIZE); ++i) {
			(*temp_ptr_out) = (*temp_ptr_in);

			++temp_ptr_in;
			temp_ptr_out += column_elements_M;
		}
	}
}
