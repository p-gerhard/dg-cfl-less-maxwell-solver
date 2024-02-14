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
#include "explicit.h"


static void assemble_rk3(gdn_simulation *simu, gdn_sparse *dg_mass,
						 gdn_sparse *dg_flux[])
{
	const bool is_inverse = true;
	simulation_assemble_mass(simu, is_inverse, dg_mass);
	simulation_assemble_internal(simu, dg_flux);
	simulation_assemble_flux(simu, dg_flux[0]);
	simulation_extend_boundaries(simu, dg_mass, dg_flux);

	/* Allocate compressed-row format and destroy the triplet format */
	sparse_csr_allocate(dg_mass, true);
	sparse_csr_allocate(dg_flux[0], true);
}

void dg_solve_rk3(gdn_simulation *simu)
{
	const gdn_real A[3] = { 0., -5. / 9., -153. / 128. };
	const gdn_real B[3] = { 1. / 3., 15. / 16., 8. / 15. };
	const gdn_real C[3] = { 0., 1. / 3., 3. / 4. };
	const int size = simu->wlen;
	const int iter_feq = 100;
	const int order = 3;
	const gdn_real tmax = simu->tmax;
	const gdn_real dt = simu->dt;

	gdn_real t = 0;
	int iter = 0;

	gdn_real *k2 = calloc(size, sizeof(gdn_real));
	gdn_real *dw = calloc(size, sizeof(gdn_real));
	gdn_real *temp = calloc(size, sizeof(gdn_real));
	gdn_real *k1 = simu->wn;

	gdn_sparse *dg_mass = sparse_coo_init(size);

	int nb_mat = 1;
	gdn_sparse *dg_flux[nb_mat];
	dg_flux[0] = sparse_coo_init(size);

	struct timeval start;
	tic(&start);
	assemble_rk3(simu, dg_mass, dg_flux);
	simu->cpu_assembly_t = toc(&start);
	
	dg_init_sol_macro(simu);

	tic(&start);
	while (t < tmax - 1e-12) {
		simu->t = t;
		/* RK3 STEP */
		for (int i = 0; i < order; i++) {
			gdn_real t = simu->t + C[i] * dt;
			dg_update_boundary(simu, t);
			sparse_spmdv(dg_flux[0], k1, temp);
			sparse_spmdv(dg_mass, temp, dw);

			if (i == 0) {
#pragma omp parallel for
				for (int j = 0; j < size; j++) {
					k2[j] = -dt * dw[j];
				}
			} else {
#pragma omp parallel for
				for (int j = 0; j < size; j++) {
					k2[j] = A[i] * k2[j] - dt * dw[j];
				}
			}
#pragma omp parallel for
			for (int j = 0; j < size; j++) {
				k1[j] = k1[j] + B[i] * k2[j];
			}
		}

		if (iter % iter_feq == 0) {
			printf("\r[Info] t = %f on %f (s)", t, tmax);
			fflush(stdout);
		}
		iter += 1;
		t += dt;
	}
	simu->t = t;
	simu->cpu_running_t = toc(&start);

	assert(fabs(t - tmax) < 1e-12);
	assert(iter == simu->itermax);

	sparse_free(dg_flux[0]);
	sparse_free(dg_mass);

	free(k2);
	free(dw);
	free(temp);
}