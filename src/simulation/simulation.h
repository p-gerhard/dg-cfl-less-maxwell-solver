/*
 * Copyright (c) 
 * 2019-2021 Philippe Helluy <helluy@math.unistra.fr>
 * 2020-2021 Pierre Gerhard <pierre.gerhard@gmail.com>
 *
 * gdon3d is free software; you can redistribute it and/or modify
 * it under the terms of the GPLv3 license. See LICENSE for details.
 */
#ifndef SIMULATION_H
#define SIMULATION_H

#include <gdon3d.h>
#include <maths/sparse.h>
#include <mesh/mesh.h>
#include <models/model.h>

typedef struct gdn_simulation {
	int m;
	int neq;
	int wlen;
	int itermax;

	bool use_kinetic_scheme;
	gdn_real t;
	gdn_real tmax;
	gdn_real dt;
	gdn_real cfl;
	gdn_real hmin;
	
	gdn_real cpu_running_t;
	gdn_real cpu_assembly_t;

	gdn_real *wn;

	gdn_model *mdl;
	gdn_mesh *tt;

} gdn_simulation;

void simulation_get_f(gdn_simulation *simu, const int id_pg, gdn_real *f);

void simulation_set_f(gdn_simulation *simu, const int id_pg, const gdn_real *f);

void simulation_relax_f(gdn_simulation *simu, const int id_pg,
						const gdn_real *f_old, const gdn_real *f_new);

int simulation_get_varindex(int m, int ipg, int iv);

void simulation_init(gdn_mesh *mesh, gdn_model *model, gdn_simulation *simu);

void simulation_free(gdn_simulation *simu);

void simulation_set_time_parameters(gdn_simulation *simu, const gdn_real tmax,
									const gdn_real dt, const gdn_real cfl);

void simulation_display_info(gdn_simulation *simu);

/* Assembly routine for non kinetic solver */

void simulation_assemble_mass(gdn_simulation *simu, const bool is_inverse,
							  gdn_sparse *dg_mass);

void simulation_assemble_inverse_mass(gdn_simulation *f, gdn_sparse *dg_mass);

void simulation_assemble_flux(gdn_simulation *f, gdn_sparse *dg_flux);

void simulation_extend_boundaries(gdn_simulation *simu, gdn_sparse *dg_mass,
								  gdn_sparse *dg_flux[]);

/* Assembly routine for kinetic relaxation solver */
void simulation_assemble_inverse_mass_kinetic(gdn_simulation *simu,
											  gdn_sparse *dg_mass);

void simulation_assemble_internal(gdn_simulation *simu, gdn_sparse *dg_flux[]);

void simulation_assemble_flux_kinetic(gdn_simulation *simu,
									  gdn_sparse *dg_flux[]);

void simulation_extend_for_boundaries_kinetic(gdn_simulation *simu,
											  gdn_sparse *dg_mass,
											  gdn_sparse *dg_flux[]);

/* Common routine */
gdn_real simulation_error_l2(gdn_simulation *simu);

void simulation_dump_info(const int num_run, gdn_simulation *simu,
						  const gdn_real error, const gdn_real order);

#endif