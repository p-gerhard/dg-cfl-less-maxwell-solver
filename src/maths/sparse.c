/*
 * Copyright (c) 
 * 2019-2021 Philippe Helluy <helluy@math.unistra.fr>
 * 2020-2021 Pierre Gerhard <pierre.gerhard@gmail.com>
 *
 * gdon3d is free software; you can redistribute it and/or modify
 * it under the terms of the GPLv3 license. See LICENSE for details.
 */
#undef NDEBUG
#include <assert.h>
#include <stdio.h>
#include <string.h>

#include <gdon3d.h>
#include "sparse.h"

#include "cs.h" /* Csparse in Suitesparse lib*/
#include "klu.h" /* KLU in lib Suitesparse lib*/

/**
 * \fn sparse_coo_init
 * \brief Allocate and initialize a new gdn_sparse struct.
 *        By default this struct contains a cs_di struct which represents 
 *        in triplet format, a unit matrix of dimension 1x1. 
 *        This matrix will grow up later when we will insert entries.
 *        The maximal growth of the matrix is controlled by the input parameter 
 *        n defines which defines the maximal number of row and column 
 *        of the matrix.
 * \param[in] n  Maximal number of row and column of the matrix 
 * \return       Pointer on the new allocated gdn_struct
 */
gdn_sparse *sparse_coo_init(int n)
{
	gdn_sparse *M = (gdn_sparse *)malloc(sizeof(gdn_sparse));

	M->neq = n;
	M->coo = cs_di_spalloc(0, 0, 1, 1, 1);

	M->is_coo_alloc = true;
	M->is_csr_alloc = false;
	M->csr = NULL;

	M->is_solver_init = false;
	M->is_lu_alloc = false;
	M->symbolic = NULL;
	M->numeric = NULL;

	klu_defaults(&M->common);
	M->common.tol = 0.1;
	M->is_init = true;

	return M;
}

/**
 * \fn sparse_coo_init_identity
 * \brief Allocate and initialize a new gdn_sparse struct.
 *        By default this struct contains a cs_di struct which represents 
 *        in triplet format, a identity matrix of dimension n x n. 
 * \param[in] n  Size of the identity matrix
 * \return       Pointer on the new allocated gdn_struct
 */
gdn_sparse *sparse_coo_init_identity(int n)
{
	gdn_sparse *Id = sparse_coo_init(n);

	if (Id->coo) {
		cs_di_spfree(Id->coo);
		Id->coo = NULL;
		Id->is_coo_alloc = false;
	}

	Id->coo = cs_di_spalloc(n, n, n, 1., 1);
	Id->is_coo_alloc = true;

	for (int i = 0; i < n; i++) {
		cs_di_entry(Id->coo, i, i, 1);
	}

	return Id;
}

/**
 * \fn sparse_init_copy
 * \brief Allocate and initialize a new gdn_sparse struct as a full copy of 
 *        the one passed in argument.
 * \param[in] M  Pointer on the gdn_struct to copy 
 * \return       Pointer on the new copy of the one passed in argument
 */
gdn_sparse *sparse_init_copy(gdn_sparse *M)
{
	assert(M != NULL && M->is_init);
	gdn_sparse *C = sparse_coo_init(M->neq);

	if (M->is_coo_alloc) {
		for (int k = 0; k < M->coo->nzmax; k++) {
			sparse_coo_insert(C, M->coo->i[k], M->coo->p[k], M->coo->x[k]);
		}
	}

	if (M->is_csr_alloc) {
		C->csr = cs_di_add(M->csr, M->csr, 1, 0);
		C->is_csr_alloc = true;
	}

	if (M->is_solver_init) {
		sparse_solver_init(C);
	}

	if (M->is_lu_alloc) {
		sparse_solver_lu_factor(C);
	}

	return C;
}

/**
 * \fn sparse_free
 * \brief Free completely a gdn_struct.
 * \param[in] M Pointer on a gdn_struct
 */
void sparse_free(gdn_sparse *M)
{
	assert(M->is_init);

	if (M->coo) {
		cs_di_spfree(M->coo);
		M->coo = NULL;
		M->is_coo_alloc = false;
	}

	if (M->csr) {
		cs_di_spfree(M->csr);
		M->csr = NULL;
		M->is_csr_alloc = false;
	}

	if (M->is_solver_init) {
		klu_free_symbolic(&M->symbolic, &M->common);
		M->symbolic = NULL;
		M->is_solver_init = false;
	}

	if (M->is_lu_alloc) {
		klu_free_numeric(&M->numeric, &M->common);
		M->numeric = NULL;
		M->is_lu_alloc = false;
	}

	M->is_init = false;
	free(M);
	M = NULL;
}

/**
 * \fn sparse_coo_insert
 * \brief Add an entry of value 1 at position (i, j) in the triplet 
 *        representation (coo) of the matrix.
 * \param[inout] M Pointer on a gdn_struct
 * \param[in] i    Row index 
 * \param[in] j    Column index
 */
void sparse_coo_insert(gdn_sparse *M, int i, int j, gdn_real val)
{
	assert(M->is_init && M->is_coo_alloc);
	cs_di_entry(M->coo, i, j, val);
}

/**
 * \fn sparse_csr_allocate
 * \brief Allocate and create the compressed-row (CSR) format of the triplet 
 *        matrix stored in the gdn_sparse struct passed in argument.
 * \param[inout] M     Pointer on a gdn_struct
 * \param[in] free_coo Boolean to delete the triplet (coo) representation
 */
void sparse_csr_allocate(gdn_sparse *M, bool free_coo)
{
	assert(M->is_init && M->is_coo_alloc);

	/* Transpose COO */
	int *tmp = M->coo->i;
	M->coo->i = M->coo->p;
	M->coo->p = tmp;

	/* Destroy the current CSR */
	if (M->is_csr_alloc) {
		cs_di_spfree(M->csr);
		M->csr = NULL;
		M->is_csr_alloc = false;
	}

	/* CSR = compressed-column (CSC) format of the transpose of the COO. */

	/* Build CSC */
	M->csr = cs_di_compress(M->coo);
	assert(M->csr);
	/* Remove duplicate entries by reduction */
	cs_di_dupl(M->csr);
	M->is_csr_alloc = true;

	assert((M->csr->m == M->neq) && (M->csr->n == M->neq));

	/* Transpose back the COO */
	tmp = M->coo->i;
	M->coo->i = M->coo->p;
	M->coo->p = tmp;

	if (free_coo) {
		cs_di_spfree(M->coo);
		M->coo = NULL;
		M->is_coo_alloc = false;
	}
}

/**
 * \fn sparse_csr_add_entry
 * \brief Add a value (NOT overwrite) to an existing entry at position (i,j) 
 *        in the CSR format of the matrix stored in the gdn_sparse struct 
 *        passed in argument. 
 * 
 *        Warning : - This function does no create a new non zero entry in the 
 *                    CSR format.
 *                  - There is not control to check if insertion was done.
 *                  - This function does no touch the COO format. 
 *                  - Using this function leads to an incoherent state between 
 *                    CSR and COO formats. 
 *        
 * \param[inout] M  Pointer on a gdn_struct
 * \param[in] i     Row index of the existing entry
 * \param[in] j     Column index of the existing entry
 * \param[in] val   Value to add
 */
void sparse_csr_add_entry(gdn_sparse *M, int i, int j, gdn_real val)
{
	assert(M->is_init && M->is_csr_alloc);
	assert(i < M->neq && j < M->neq);

	for (int k = M->csr->p[i]; k < M->csr->p[i + 1]; k++) {
		if (M->csr->i[k] == j) {
			M->csr->x[k] += val;
			return;
		}
	}
}

/**
 * \fn sparse_csr_set_entry
 * \brief Set a value (overwrite) to an existing entry at position (i,j) 
 *        in the CSR format of the matrix stored in the gdn_sparse struct 
 *        passed in argument. 
 * 
 *        Warning : - This function does no create a new non zero entry in the 
 *                    CSR format.
 *                  - There is not control to check if insertion was done.
 *                  - This function does no touch the COO format. 
 *                  - Using this function leads to an incoherent state between 
 *                    CSR and COO formats. 
 *        
 * \param[inout] M  Pointer on a gdn_struct
 * \param[in] i     Row index of the existing entry
 * \param[in] j     Column index of the existing entry
 * \param[in] val   Value to set
 */
void sparse_csr_set_entry(gdn_sparse *M, int i, int j, gdn_real val)
{
	assert(M->is_init && M->is_csr_alloc);
	assert(i < M->neq && j < M->neq);

	for (int k = M->csr->p[i]; k < M->csr->p[i + 1]; k++) {
		if (M->csr->i[k] == j) {
			M->csr->x[k] = val;
			return;
		}
	}
}

/**
 * \fn sparse_csr_get_entry
 * \brief Get the value of an existing entry at position (i,j) 
 *        in the CSR format of the matrix stored in the gdn_sparse struct 
 *        passed in argument. 
 * 
 *        Warning : If the entry doesn't exist it returns 0 anyway.

 * \param[inout] M  Pointer on a gdn_struct
 * \param[in] i     Row index of the existing entry
 * \param[in] j     Column index of the existing entry
 * \param[in] val   Value to set
 */
gdn_real sparse_csr_get_entry(gdn_sparse *M, int i, int j)
{
	assert(M->is_init && M->is_csr_alloc);
	assert(i < M->neq && j < M->neq);

	for (int k = M->csr->p[i]; k < M->csr->p[i + 1]; k++) {
		if (M->csr->i[k] == j) {
			return M->csr->x[k];
		}
	}
	return 0.;
}

/**
 * \fn sparse_csr_matrix_multiplication
 * \brief Matrix multiplication (C = A.B) between the two CSR matrices stored 
 *        in gdn_struct passed in argument.
 * 
 *        Warning : - This function allocates a new gdn_sparse struct which 
 *                    contains the results of the multiplication.
 *                  - This function does not create the COO format in the 
 *                    returned struct.
 *    
 * \param[in] A  Pointer on a gdn_struct 
 * \param[in] B  Pointer on a gdn_struct
 * \return       Pointer on a gdn_struct containing the result of the CSR matrix
 *               product
 */
gdn_sparse *sparse_csr_matrix_multiplication(const gdn_sparse *A,
											 const gdn_sparse *B)
{
	assert(A->is_init && A->is_csr_alloc);
	assert(B->is_init && B->is_csr_alloc);

	gdn_sparse *M = sparse_coo_init(A->neq);

	/* Back to CSC (Transpose of CSR) format */
	cs_di *AT = cs_di_transpose(A->csr, 1);
	cs_di *BT = cs_di_transpose(B->csr, 1);

	cs_di *C = cs_di_multiply(AT, BT);

	cs_di_spfree(AT);
	cs_di_spfree(BT);

	/* Back to CSR (Transpose of CSC) format */
	M->csr = cs_di_transpose(C, 1);
	cs_di_spfree(C);
	cs_di_dupl(M->csr);

	M->is_csr_alloc = true;
	return M;
}

/**
 * \fn sparse_csr_matrix_addition
 * \brief Matrix addition (C = alpha*A + beta*B) between the two CSR matrices
 *        stored in gdn_struct passed in argument.
 * 
 *        Warning : - This function allocates a new gdn_sparse struct which 
 *                    contains the results of the addition.
 *                  - This function does not create the COO format in the 
 *                    returned struct.
 *    
 * \param[in] A      Pointer on a gdn_struct 
 * \param[in] B      Pointer on a gdn_struct
 * \param[in] alpha  Coefficient associated to A
 * \param[in] beta   Coefficient associated to B
 * \return           Pointer on a gdn_struct containing the result of the CSR 
 *                   matrix addiction
 */
gdn_sparse *sparse_csr_matrix_addition(const gdn_sparse *A, const gdn_sparse *B,
									   const gdn_real alpha,
									   const gdn_real beta)
{
	assert(A->is_init && A->is_csr_alloc);
	assert(B->is_init && B->is_csr_alloc);

	gdn_sparse *C = sparse_coo_init(A->neq);

	C->csr = cs_di_add(A->csr, B->csr, alpha, beta);
	cs_di_dupl(C->csr);
	C->is_csr_alloc = true;

	return C;
}

/**
 * \fn sparse_spmdv
 * \brief Sparse Matrix Dense Vector (SPMDV) product. Computes the product
 *        between the CSR matrix passed in argument and the dense vector
 *        x, (ie. res = A.x). This function uses OpenMP.
 * 
 * \param[in] A       Pointer on a gdn_struct 
 * \param[in] x       Vector used in the prodcut
 * \param[inout] res  Vector containing the result
 */
void __attribute__((optimize("Ofast")))
sparse_spmdv(const gdn_sparse *A, const gdn_real *x, gdn_real *res)
{
	gdn_real res_i = 0;

#pragma omp parallel for
	for (int i = 0; i < A->neq; i++) {
		res_i = 0.;

		for (int j = A->csr->p[i]; j < A->csr->p[i + 1]; j++) {
			res_i += A->csr->x[j] * x[A->csr->i[j]];
		}

		res[i] = res_i;
	}
}

/**
 * \fn sparse_solver_init
 * \brief Initialise the solver (KLU) struct on the CSR matrix stored in
 *        the gdn_struct passed in argument.
 * 
 * \param[in] M Pointer on a gdn_struct 
 */
void sparse_solver_init(gdn_sparse *M)
{
	assert(M->is_init);

	const bool free_coo = true;

	if (!M->is_csr_alloc) {
		sparse_csr_allocate(M, free_coo);
	}
	assert(M->is_csr_alloc);

	M->symbolic = klu_analyze(M->neq, M->csr->p, M->csr->i, &(M->common));
	assert(M->symbolic);
	M->is_solver_init = true;
}

/**
 * \fn sparse_solver_lu_factor
 * \brief Compute the LU factorisation (using KLU) of the CSR matrix stored in
 *        the gdn_struct passed in argument.
 * 
 * \param[in] M Pointer on a gdn_struct 
 */
void sparse_solver_lu_factor(gdn_sparse *M)
{
	assert(M->is_init && M->is_solver_init);

	if (!M->is_solver_init) {
		sparse_solver_init(M);
	}

	assert(M->is_solver_init);

	M->numeric =
		klu_factor(M->csr->p, M->csr->i, M->csr->x, M->symbolic, &(M->common));

	// klu_sort(M->symbolic, M->numeric, &(M->common));

	assert(M->numeric);
	M->is_lu_alloc = true;
}

/**
 * \fn sparse_solve
 * \brief Solve the linear system A.X = rhs using LU factorisation of the CSR 
 *        matrix stored in the gdn_struct passed in argument.
 * 
 * \param[in] A    Pointer on a gdn_struct 
 * \param[in] rhs  Right hand side vector 
 * \param[inout] x Solution vector 
 */
void sparse_solve(gdn_sparse *A, gdn_real *rhs)
{
	/* klu_solve deals with CSC format, for CSR format we use klu_tsolve */
	klu_tsolve(A->symbolic, A->numeric, A->neq, 1, rhs, &(A->common));
}