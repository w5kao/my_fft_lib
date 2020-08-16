/* This is library for parallelizing of simple linear cycles */
/* v 0.1 Copyright (c) 2009 Alexander O. Korotkevich, w5kao@yahoo.com */
/***********************************************************************
 * This file is part of my_thr_lib.

    my_thr_lib is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License or Lesser General Public License
    as published by the Free Software Foundation, either version 3 of the License, or
    any later version.

    my_thr_lib is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with my_thr_lib.  If not, see <https://www.gnu.org/licenses/>. 
*************************************************************************/
#include "my_thr_lib.h"

unsigned long int thr_number, finished_thr_number, initialized_thr_number;

void (*ptr_thr_function) (void *);

struct thread_input *thread_data;

pthread_mutex_t finished_thr_number_mutex, initialized_thr_number_mutex;

pthread_cond_t	thread_manager_wait;

unsigned long int thr_i;

void my_thr_pool_init (unsigned long int thread_number) {
	thr_number = thread_number;
	
	thread_data = malloc (thread_number * sizeof(struct thread_input));

	pthread_mutex_init (&finished_thr_number_mutex, NULL);
	pthread_mutex_init (&initialized_thr_number_mutex, NULL);
	pthread_cond_init (&thread_manager_wait, NULL);

	initialized_thr_number = 0;
	for (thr_i = 0; thr_i < thread_number; ++thr_i) {
#ifdef DEBUG
		printf ("Initializing thread # %lu\n", thr_i); fflush(stdout);
#endif
		pthread_mutex_init (&(thread_data [thr_i].run), NULL);
		pthread_cond_init (&(thread_data [thr_i].wait), NULL);

		(thread_data [thr_i].thr_index) = thr_i;
		pthread_create (&(thread_data [thr_i].th), NULL,
			elementary_thread, (void*) &(thread_data [thr_i]));
	}
	
	pthread_mutex_lock (&initialized_thr_number_mutex);
	while (initialized_thr_number < thr_number) 
		pthread_cond_wait(&thread_manager_wait, &initialized_thr_number_mutex);
	pthread_mutex_unlock (&initialized_thr_number_mutex);
	printf ("%lu Threads initialized.\n", thread_number); fflush(stdout);
}

void my_thr_data_assign (unsigned long int thr_index, void * data_input) {
	
	(thread_data [thr_index].data_in) = data_input;
#ifdef DEBUG
	printf ("Data for thread # %lu assigned.\n", thr_index); fflush(stdout);
#endif
}

void my_thr_pool_clear (void) {
	for (thr_i = 0; thr_i < thr_number; ++thr_i) {
#ifdef DEBUG
		printf ("Clearing thread # %lu\n", thr_i); fflush(stdout);
#endif
		pthread_cancel (thread_data [thr_i].th);
		pthread_mutex_destroy (&(thread_data [thr_i].run));
		pthread_cond_destroy (&(thread_data [thr_i].wait));
	}
	free (thread_data);
	pthread_mutex_destroy (&finished_thr_number_mutex);
	pthread_mutex_destroy (&initialized_thr_number_mutex);
	pthread_cond_destroy (&thread_manager_wait);

	printf ("Threads cleared.\n");
}

void* elementary_thread (void *thr_input) {
	struct thread_input *input;
	
	input = (struct thread_input *) thr_input;
#ifdef DEBUG
	printf ("Entered thread # %lu.\n", input->thr_index); fflush(stdout);
#endif
	pthread_mutex_lock (&initialized_thr_number_mutex);
	++initialized_thr_number;
#ifdef DEBUG
	printf ("Thread # %lu initialized. Number of initialized threads = %lu\n",
			input->thr_index, initialized_thr_number); fflush(stdout);
#endif
	pthread_cond_signal(&thread_manager_wait);
	pthread_mutex_lock (&(input->run));
	pthread_mutex_unlock (&initialized_thr_number_mutex);

	do {
#ifdef DEBUG
		printf ("Waiting for signal for thread # %lu.\n", input->thr_index); fflush(stdout);
#endif
		pthread_cond_wait (&(input->wait), &(input->run));
		pthread_mutex_unlock (&(input->run));
#ifdef DEBUG
		printf ("Got signal for thread # %lu.\n", input->thr_index); fflush(stdout);
#endif

		ptr_thr_function (input->data_in);

		pthread_mutex_lock (&finished_thr_number_mutex);
		++finished_thr_number;
		pthread_cond_signal (&thread_manager_wait);
		pthread_mutex_lock (&(input->run));
		pthread_mutex_unlock (&finished_thr_number_mutex);

	} while (1);
	return  (void*) 0;
}

void my_thr_manager (void (*thr_function_requested) (void*)) {

	finished_thr_number = 0;
#ifdef DEBUG
	printf ("Entered manager.\n"); fflush(stdout);
#endif
	ptr_thr_function = thr_function_requested;

	for (thr_i = 0; thr_i < thr_number; ++thr_i) {
#ifdef DEBUG
		printf ("Signaling to thread # %lu.\n", thr_i); fflush(stdout);
#endif
		pthread_mutex_lock (&(thread_data[thr_i].run));
		pthread_cond_signal(&(thread_data[thr_i].wait));
		pthread_mutex_unlock (&(thread_data[thr_i].run));
	}
#ifdef DEBUG
	printf ("Threads started. Waiting for finish.\n"); fflush(stdout);
#endif
	pthread_mutex_lock (&finished_thr_number_mutex);
	while (finished_thr_number < thr_number) 
		pthread_cond_wait(&thread_manager_wait, &finished_thr_number_mutex);
	pthread_mutex_unlock (&finished_thr_number_mutex);

#ifdef DEBUG
	printf ("Threads finished. Exiting.\n"); fflush(stdout);
#endif
}
