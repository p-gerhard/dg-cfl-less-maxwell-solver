/*
 * Copyright (c) 
 * 2019-2021 Philippe Helluy <helluy@math.unistra.fr>
 * 2020-2021 Pierre Gerhard <pierre.gerhard@gmail.com>
 *
 * gdon3d is free software; you can redistribute it and/or modify
 * it under the terms of the GPLv3 license. See LICENSE for details.
 */
#ifndef SPARSE_H
#define SPARSE_H

#include <stdbool.h>

#include <gdon3d.h>

#include "cs.h" /* CSparse in Suitesparse lib*/
#include "klu.h" /* KLU in Suitesparse lib    */

#ifndef _DOUBLE_PRECISION
#define cs_di cs
#endif

typedef struct gdn_sparse {
	int neq;
	bool is_init;

	bool is_coo_alloc;
	bool is_csr_alloc;
	bool is_lu_alloc;

	cs_di *coo; /* CSparse Triplet format        */
	cs_di *csr; /* CSparse Compressed-row format */

	/* Handler for KLU solver */
	bool is_solver_init;
	klu_symbolic *symbolic;
	klu_numeric *numeric; /* KLU LU factorisation           */
	klu_common common;

} gdn_sparse;

gdn_sparse *gdn_sparse_malloc();

gdn_sparse *sparse_coo_init(int n);

gdn_sparse *sparse_coo_init_identity(int n);

gdn_sparse *sparse_init_copy(gdn_sparse *M);

void sparse_free(gdn_sparse *M);

void gdn_sparse_coo_switch_on(gdn_sparse *M, int i, int j);

void sparse_coo_insert(gdn_sparse *M, int i, int j, gdn_real val);

void sparse_csr_allocate(gdn_sparse *M, bool free_coo);

void sparse_csr_add_entry(gdn_sparse *M, int i, int j, gdn_real val);

void sparse_csr_set_entry(gdn_sparse *M, int i, int j, gdn_real val);

gdn_real sparse_csr_get_entry(gdn_sparse *M, int i, int j);

gdn_sparse *sparse_csr_matrix_multiplication(const gdn_sparse *A,
											 const gdn_sparse *B);

gdn_sparse *sparse_csr_matrix_addition(const gdn_sparse *A, const gdn_sparse *B,
									   const gdn_real alpha,
									   const gdn_real beta);

void sparse_spmdv(const gdn_sparse *A, const gdn_real *x, gdn_real *res);

void sparse_solver_init(gdn_sparse *A);

void sparse_solver_lu_factor(gdn_sparse *M);

void sparse_solve(gdn_sparse *A, gdn_real *rhs);

#endif