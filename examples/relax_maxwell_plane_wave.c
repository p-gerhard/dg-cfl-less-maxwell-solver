/*
 * Copyright (c) 
 * 2019-2021 Philippe Helluy <helluy@math.unistra.fr>
 * 2020-2021 Pierre Gerhard <pierre.gerhard@gmail.com>
 *
 * gdon3d is free software; you can redistribute it and/or modify
 * it under the terms of the GPLv3 license. See LICENSE for details.
 */
#undef NDEBUG
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <time.h>

#include <dg/dg.h>
#include <models/model.h>
#include <models/maxwell.h>
#include <models/d3q4.h>
#include <io/io_hdf5.h>
#include <simulation/simulation.h>
#include <utils/timing.h>
#include <utils/convergence.h>


static char base_name[256] = "maxwell_rk3_plane_wave";

/**
 * \fn imposed_macro_cos
 * \brief Imposed oscillatory solution of frequence nu for Maxwell's equations.
 *        This solution is constructed using the method of characteristics, 
 *        therby it is an exact solution for all (x, t) of Maxwell's equations.
 * \param[in] x    Position vector
 * \param[in] t    Time
 * \param[inout] w Macro state
 */
void static imposed_macro_cos(const gdn_real *x, const gdn_real t, gdn_real *w)
{
	// const gdn_real nu = 0.5;
	const gdn_real nu = 1;
	const gdn_real s = cos(nu * M_PI * (x[0] - t));
	// const gdn_real s = (x[0] - t) * (x[0] - t);
	w[0] = 0;
	w[1] = 0;
	w[2] = s;
	w[3] = 0;
	w[4] = -s;
	w[5] = 0;
}

int relax_maxwell_plane_wave(void)
{
	bool export_xdmf = true;

	gdn_mesh mesh = { 0 };
	gdn_simulation simu = { 0 };
	gdn_model model = { 0 };

	model.nb_w = 6;
	model.get_num_flux = maxwell_num_flux_upwind;
	model.get_num_flux_boundary = maxwell_num_flux_upwind;
	model.get_imposed_data = imposed_macro_cos;
    model.get_proj_pos     = maxwell_get_proj_pos_null_prod_w;
    model.get_proj_neg     = maxwell_get_proj_neg_prod_w;
    model.use_relax_scheme = true;
    d3q4_set_model(&model);

	char *mesh_name = "../data/mesh/t4_cube_8.msh";

	/* Time settings */
	int mesh_raf = 9;
	gdn_real tmax = 1.0;
	gdn_real dt = 0.032 * (1. / mesh_raf);
	
	export_xdmf = true;

	utils_l2_time_convergence(&model, 4, mesh_name, dt, tmax, export_xdmf);
	
	return 1;
}

int main(void)
{
	int resu = relax_maxwell_plane_wave();
	if (resu)
		printf("[Info] - relax_maxwell_plane_wave : OK\n");
	else
		printf("[Info] -  relax_maxwell_plane_wave : FAILED\n");
	return !resu;
}