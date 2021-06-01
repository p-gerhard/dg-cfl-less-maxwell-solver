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

#include <mesh/t10.h>
#include <models/model.h>
#include <maths/sparse.h>

#include "dg.h"
#include <gdon3d.h>
#include "implicit.h"
#include <simulation/simulation.h>

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

void dg_implicit_solve_theta_scheme_relax_v2(gdn_simulation *simu)
{
	const gdn_real dt = simu->dt;
	const gdn_real tmax = simu->tmax;
	const gdn_real theta = 0.5;
	const int iter_feq = 1;
	gdn_real t = simu->t;
	int iter = 0;
	gdn_real *tmp_ptr = NULL;
	int mat_nb;

	if (simu->use_kinetic_scheme) {
		mat_nb = simu->mdl->nb_v;
	} else {
		mat_nb = 1;
	}

	/* DEBUG */
	// char filename[1024] = "debug_";
	// char mesh_filename[1024] = {0};
	// char data_filename[1024] = {0};
	// char xmf_filename[1024]  = {0};
	// set_name_mesh(mesh_filename, filename);
	// io_hdf5_dump_mesh(simu, mesh_filename);

	/* Warning allocation on the stack using VLA */
	gdn_sparse *M1[mat_nb];
	gdn_sparse *M2[mat_nb];

	/* Temporary buffer */
	gdn_real *sol = (gdn_real *)malloc(simu->wlen * sizeof(gdn_real));
	assemble_theta_scheme(simu, theta, dt, M1, M2);
	dg_init_sol_macro(simu);
	while (t < tmax - 1e-12) {
		/* DEBUG */
		// set_name_data(data_filename, filename, iter);
		// set_name_xmf(xmf_filename,   filename, iter);
		// io_hdf5_dump_data(simu, data_filename);
		// io_hdf5_write_xdmf(simu, xmf_filename, mesh_filename,
		// data_filename);

		dg_update_boundary(simu, t + dt / 2);
		apply_transport_step_theta_scheme(simu, sol, M1, M2);

		/* Adress swap in order to avoid a buffer recopy in sparse_solve */
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
	printf("\n");
	dg_update_boundary(simu, t);

	assert(fabs(t - tmax) < 1e-12);
	assert(iter == simu->itermax);

	for (int k = 0; k < mat_nb; k++) {
		sparse_free(M1[k]);
		sparse_free(M2[k]);
	}
	free(sol);
}

void dg_implicit_solve_theta_scheme_relax_v1(gdn_simulation *simu)
{
	const gdn_real dt = simu->dt;
	const gdn_real tmax = simu->tmax;
	const gdn_real theta = 0.5;
	const int iter_feq = 1;
	gdn_real t = simu->t;
	gdn_real *tmp_ptr = NULL;
	int iter = 0;

	int mat_nb;
	if (simu->use_kinetic_scheme) {
		mat_nb = simu->mdl->nb_v;
	} else {
		mat_nb = 1;
	}

	/* Warning allocation on the stack using VLA */
	gdn_sparse *M1[mat_nb];
	gdn_sparse *M2[mat_nb];

	/* Temporary buffer */
	gdn_real *sol = (gdn_real *)malloc(simu->wlen * sizeof(gdn_real));

	assemble_theta_scheme(simu, theta, dt / 2, M1, M2);
	dg_init_sol_macro(simu);
	dg_update_boundary(simu, t + dt / 2.);

	while (t < tmax - 1e-12) {
		apply_transport_step_theta_scheme(simu, sol, M1, M2);
		/* Adress swap in order to avoid a buffer recopy in sparse_solve */
		tmp_ptr = simu->wn;
		simu->wn = sol;
		sol = tmp_ptr;

		dg_apply_relaxation(simu);
		apply_transport_step_theta_scheme(simu, sol, M1, M2);
		/* Adress swap in order to avoid a buffer recopy in sparse_solve */
		tmp_ptr = simu->wn;
		simu->wn = sol;
		sol = tmp_ptr;

		dg_update_boundary(simu, t + dt);

		if (iter % iter_feq == 0) {
			printf("\r[Info] t = %f on %f (s)", t, tmax);
			fflush(stdout);
		}

		t += dt;
		simu->t = t;
		iter += 1;
	}

	dg_update_boundary(simu, t + dt / 2.);
	assert(fabs(t - tmax) < 1e-12);
	assert(iter == simu->itermax);

	for (int k = 0; k < mat_nb; k++) {
		sparse_free(M1[k]);
		sparse_free(M2[k]);
	}
	free(sol);
}