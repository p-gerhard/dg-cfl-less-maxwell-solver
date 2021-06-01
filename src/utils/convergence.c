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
#include <stdio.h>
#include <stdbool.h>

#include <gdon3d.h>
#include <dg/dg.h>
#include <io/io_hdf5.h>
#include <mesh/mesh.h>
#include <models/model.h>
#include <simulation/simulation.h>

#include "convergence.h"
#include "timing.h"

static gdn_real get_order(const gdn_real y_1, const gdn_real y_2,
						  const gdn_real x_1, const gdn_real x_2)
{
	const gdn_real tol = 1e-12;

	assert(fabs(y_1) > tol);
	assert(fabs(x_1) > tol);

	const gdn_real r1 = y_2 / y_1;
	const gdn_real r2 = x_2 / x_1;

	assert(r1 > 0);
	assert(r2 > 0);

	return log(r1) / log(r2);
}

void utils_l2_time_convergence(gdn_model *model, const int nb_run,
							   const char *msh_filename, const gdn_real init_dt,
							   const gdn_real tmax, const bool export_xdmf)
{
	gdn_real error[nb_run];
	gdn_real order[nb_run - 1];
	gdn_real dt[nb_run];
	gdn_real iter[nb_run];

	gdn_simulation simu = { 0 };
	gdn_mesh mesh = { 0 };

	/* Initial time settings */
	const int init_iter_max = (int)floor(tmax / init_dt) - 1;
	const gdn_real init_tmax = init_iter_max * init_dt;

	order[0] = -1;
	dt[0] = init_dt;
	iter[0] = init_iter_max;

	for (int k = 1; k < nb_run; k++) {
		dt[k] = 0.5 * dt[k - 1];
		iter[k] = 2 * iter[k - 1];
	}

	mesh_read_from_msh22_file(&mesh, msh_filename);
	mesh_build_connectivity(&mesh);
	simulation_init(&mesh, model, &simu);

	for (int k = 0; k < nb_run; k++) {
		simulation_set_time_parameters(&simu, init_tmax, dt[k], -1);

		assert(simu.tmax == init_tmax);
		assert(simu.itermax == iter[k]);

		dg_solve(&simu);
		error[k] = simulation_error_l2(&simu);

		if (k > 0) {
			order[k - 1] = get_order(error[k], error[k - 1], dt[k], dt[k - 1]);
		}

		simulation_dump_info(k, &simu, error[k], order[k]);

		// if (export_xdmf){
		// 	char filename[1024];
		// 	snprintf(filename, 1024, "wave_run_%d", k);
		// 	// snprintf(filename, 1024, "maxwell_cos5_d3q4_cfl_mul_%d", coef[k]);
		// 	io_save_all(&simu, filename, simu.itermax);
		// }
	}
	simulation_free(&simu);
}

void utils_l2_space_convergence(gdn_model *model, const int nb_run,
							 const char *msh_filename[], const gdn_real init_dt,
							 const gdn_real tmax,
							 const bool use_time_refinement,
							 const bool export_xdmf)
{
	gdn_real error[nb_run];
	gdn_real order;
	gdn_real dt[nb_run];
	gdn_real dh[nb_run];
	gdn_real iter[nb_run];

	gdn_simulation simu = { 0 };
	gdn_mesh mesh = { 0 };

	/* Initial time settings */
	const int init_iter_max = (int)floor(tmax / init_dt) - 1;
	const gdn_real init_tmax = init_iter_max * init_dt;

	dt[0] = init_dt;
	iter[0] = init_iter_max;

	if (use_time_refinement) {
		for (int k = 1; k < nb_run; k++) {
			dt[k] = 0.5 * dt[k - 1];
			iter[k] = 2 * iter[k - 1];
		}
	} else {
		for (int k = 1; k < nb_run; k++) {
			dt[k] = init_dt;
			iter[k] = init_iter_max;
		}
	}

	for (int k = 0; k < nb_run; k++) {
		mesh_read_from_msh22_file(&mesh, msh_filename[k]);
		mesh_build_connectivity(&mesh);
		simulation_init(&mesh, model, &simu);
		simulation_set_time_parameters(&simu, init_tmax, dt[k], -1);

		assert(simu.tmax == init_tmax);
		assert(simu.itermax == iter[k]);

		dh[k] = simu.hmin;
		dg_solve(&simu);

		error[k] = simulation_error_l2(&simu);

		if (k > 0) {
			assert(fabs((dh[k - 1] / dh[k]) - 2) < 1e-6);
			order = get_order(error[k], error[k - 1], dh[k], dh[k - 1]);
		}
		simulation_dump_info(k, &simu, error[k], order);

		// if (export_xdmf) {
		// 	char filename[1024];
		// 	snprintf(filename, 1024, "test_order_%s_run_%d", model->name, k);
		// 	io_save_all(&simu, filename, simu.itermax);
		// }

		simulation_free(&simu);
	}
}

// void utils_benchmark_run(gdn_model *model, const int nb_run,
// 						 const char *mesh_name[], const int itermax,
// 						 const gdn_real dt)
// {
// 	gdn_simulation simu = { 0 };
// 	gdn_mesh mesh = { 0 };

// 	/* Initial time settings */
// 	const gdn_real init_dt = dt;

// 	/* Adjust initial tmax as a multiple of dt */
// 	gdn_real init_tmax = itermax * dt;

// 	/* Warning allocation on the stack using VLA */
// 	gdn_real cpu_list[nb_run];
// 	int nb_elem_list[nb_run];
// 	int nb_inc[nb_run];

// 	struct timeval start;
// 	for (int k = 0; k < nb_run; k++) {
// 		mesh_read_from_msh22_file(&mesh, mesh_name[k]);
// 		mesh_build_connectivity(&mesh);
// 		simulation_init(&mesh, model, &simu);

// 		simulation_set_time_parameters(&simu, init_tmax, dt, -1);
// 		simulation_display_info(&simu);

// 		tic(&start);
// 		dg_solve(&simu);
// 		cpu_list[k] = toc(&start);
// 		nb_elem_list[k] = simu.tt->nbelems;
// 		nb_inc[k] = simu.m * simu.tt->nbelems * 10;
// 		simulation_free(&simu);
// 		printf("\n");
// 	}

// 	/* Display order */
// 	printf("Sumary : \n");
// 	printf("%-30s %-12s %-12s %-12s\n", "Filename", "nb. elem.", "nb. inc.",
// 		   "CPU (s)");

// 	for (int k = 0; k < nb_run; k++) {
// 		printf("%-30s %-12d %-12d %-12f\n", mesh_name[k], nb_elem_list[k],
// 			   nb_inc[k], cpu_list[k]);
// 	}
// }
