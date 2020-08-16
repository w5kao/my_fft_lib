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
/* This is library for parallelizing of simple linear cycles */

#ifndef _MY_THR_LIB_H
#define _MY_THR_LIB_H

#include <stdio.h>
#include <stdlib.h>

#include <pthread.h>
#include <signal.h>

#define _DEBUG

struct thread_input {
	pthread_t th;
	pthread_mutex_t run;
	pthread_cond_t wait;
	unsigned long int thr_index;
	void *data_in;
	void (*f_ptr) (void*);
};

void my_thr_pool_init (unsigned long int thread_number);
void my_thr_data_assign (unsigned long int thr_index, void * data_input);
void my_thr_pool_clear (void);

void* elementary_thread (void* thr_input);
void my_thr_manager (void (*thr_function_requested) (void*));

#endif /* _MY_THR_LIB_H */
