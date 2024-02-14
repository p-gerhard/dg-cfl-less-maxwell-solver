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
#include <math.h>
#include <sys/time.h>

#include <gdon3d.h>
#include <mesh/t10.h>
#include <models/model.h>
#include <maths/sparse.h>
#include <simulation/simulation.h>
#include <utils/timing.h>

#include "dg.h"
#include "relaxation.h"

static void assemble_theta_scheme(gdn_simulation *simu, const gdn_real theta,
								  const gdn_real dt, gdn_sparse **M1,
								  gdn_sparse **M2)
{
	int mat_size;
	int mat_nb;

	if (simu->use_kinetic_scheme) {
		mat_size = simu->neq;
		mat_nb = simu->mdl->nb_v;
	} else {
		mat_size = simu->neq * simu->m;
		mat_nb = 1;
	}

	gdn_sparse *dg_mass = sparse_coo_init(mat_size);

	/* WARNING : Allocation on the stack using VLA */
	gdn_sparse *dg_flux[mat_nb];

	for (int k = 0; k < mat_nb; k++) {
		dg_flux[k] = sparse_coo_init(mat_size);
	}

	/* DG Assembly */
	const bool is_inverse = true;
	simulation_assemble_mass(simu, is_inverse, dg_mass);
	simulation_assemble_internal(simu, dg_flux);

	simulation_assemble_flux_kinetic(simu, dg_flux);
	simulation_extend_boundaries(simu, dg_mass, dg_flux);

	/* Allocate compressed-row format and destroy the triplet format */
	bool free_coo = true;
	sparse_csr_allocate(dg_mass, free_coo);

	for (int k = 0; k < mat_nb; k++) {
		sparse_csr_allocate(dg_flux[k], free_coo);
	}

	/* Compute theta-scheme matrices :
   *  - Step 1 :
   *    Id := Identity matrix
   *    M  := Mass matrix
   *    F  := Flux matrix
   *
   *    M1 := (Id + (     theta * dt ) * (M^-1) . F)
   *    M2 := (Id - ( 1 - theta * dt ) * (M^-1) . F)
   *
   *  - Step 2 :
   *    LU Factorisation of M1
   */
	gdn_sparse *Id = sparse_coo_init_identity(mat_size);
	sparse_csr_allocate(Id, free_coo);

	for (int k = 0; k < mat_nb; k++) {
		gdn_sparse *A = sparse_csr_matrix_multiplication(dg_mass, dg_flux[k]);
		M1[k] = sparse_csr_matrix_addition(Id, A, 1, theta * dt);
		M2[k] = sparse_csr_matrix_addition(Id, A, 1, -(1 - theta) * dt);
		sparse_free(A);
		sparse_free(dg_flux[k]);
		sparse_solver_init(M1[k]);
		sparse_solver_lu_factor(M1[k]);
	}

	sparse_free(dg_mass);
	sparse_free(Id);
}

static void apply_transport_step_theta_scheme(gdn_simulation *simu,
											  gdn_real *sol, gdn_sparse **M1,
											  gdn_sparse **M2)
{
	const int neq = simu->neq;
	const int nb_w = simu->mdl->nb_w;
	const int nb_v = simu->mdl->nb_v;
	const int m = simu->m;

	/* Method 1 */
	for (int iw = 0; iw < nb_w; iw++) {
#pragma omp parallel for
		for (int iv = 0; iv < nb_v; iv++) {
			int offset = iw * (nb_v * neq) + iv * neq;
			sparse_spmdv(M2[iv], simu->wn + offset, sol + offset);
			sparse_solve(M1[iv], sol + offset);
		}
	}

	/* Methode 2 */
	// 	for (int iv = 0; iv < nb_v; iv++) {
	// #pragma omp parallel for
	// 		for (int iw = 0; iw < nb_w; iw++) {
	// 			int offset = iv * nb_w * neq + iw * neq;
	// 			sparse_spmdv(M2[iv], simu->wn + offset, sol + offset);
	// 		}
	// 	}

	// for (int iv = 0; iv < nb_v; iv++) {
	// 	for (int iw = 0; iw < nb_w; iw++) {
	// 		int offset = iv * nb_w * neq + iw * neq;
	// 		sparse_solve(M1[iv], sol + offset);
	// 	}
	// }
}

void dg_solve_relaxation(gdn_simulation *simu)
{
	const gdn_real dt = simu->dt;
	const gdn_real tmax = simu->tmax;
	const gdn_real theta = 0.5;
	gdn_real t = 0;
	gdn_real *tmp_ptr = NULL;

	const int nb_v = simu->mdl->nb_v;
	const int iter_feq = 1;
	int iter = 0;

	/* Warning allocation on the stack using VLA */
	gdn_sparse *M1[nb_v];
	gdn_sparse *M2[nb_v];

	/* Temporary buffer */
	gdn_real *sol = (gdn_real *)malloc(simu->wlen * sizeof(gdn_real));

	struct timeval start;
	tic(&start);
	assemble_theta_scheme(simu, theta, dt, M1, M2);
	simu->cpu_assembly_t = toc(&start);

	dg_init_sol_macro(simu);

	tic(&start);
	while (t < tmax - 1e-12) {
		dg_update_boundary(simu, t + dt / 2);
		apply_transport_step_theta_scheme(simu, sol, M1, M2);

		/* Shallow copy */
		tmp_ptr = simu->wn;
		simu->wn = sol;
		sol = tmp_ptr;

		dg_apply_relaxation(simu);

		if (iter % iter_feq == 0) {
			printf("\r[Info] t = %f on %f (s)", t, tmax);
			fflush(stdout);
		}
		t += dt;
		simu->t = t;
		iter += 1;
	}
	dg_update_boundary(simu, t);
	simu->cpu_running_t = toc(&start);
	assert(iter == simu->itermax);
	assert(fabs(t - tmax) < 1e-12);

	for (int k = 0; k < nb_v; k++) {
		sparse_free(M1[k]);
		sparse_free(M2[k]);
	}
	free(sol);
}
