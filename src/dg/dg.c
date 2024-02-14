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
#include <stdbool.h>

#include <gdon3d.h>
#include <simulation/simulation.h>
#include <mesh/mesh.h>
#include <mesh/t10.h>
#include <models/model.h>

#include "explicit.h"
#include "relaxation.h"

void dg_solve(gdn_simulation *simu)
{
	if (simu->use_kinetic_scheme) {
		dg_solve_relaxation(simu);
	} else {
		dg_solve_rk3(simu);
	}
}

void dg_init_sol_macro(gdn_simulation *simu)
{
	const gdn_mesh *mesh = simu->tt;
	const gdn_model *model = simu->mdl;
	const int m = simu->m;
	const int nb_w = model->nb_w;

	/* Warning allocation on the stack using VLA */
	gdn_real f_eq[m];
	gdn_real wI[nb_w];

	for (int id_elem = 0; id_elem < mesh->nbelems; id_elem++) {
		for (int id_loc = 0; id_loc < NQN; id_loc++) {
			int id_pg = id_elem * NQN + id_loc;
			int id_node = mesh->elem2qnode[id_pg];
			gdn_real *x = mesh->node + 3 * id_node;

			if (simu->use_kinetic_scheme) {
				model->get_imposed_data(x, 0., wI);
				model->get_feq(model, wI, f_eq);
				simulation_set_f(simu, id_pg, f_eq);
			} else {
				gdn_real *w = simu->wn + id_pg * m;
				model->get_imposed_data(x, 0., w);
			}
		}
	}
}

void dg_update_boundary(gdn_simulation *simu, const gdn_real t)
{
	const gdn_mesh *mesh = simu->tt;
	const gdn_model *model = simu->mdl;
	const int offset = mesh->nbelems * NQN; /* TO CHANGE */

	const int m = simu->m;
	const int nb_w = model->nb_w;

	/* Warning allocation on the stack using VLA */
	gdn_real fL[m];
	gdn_real fR[m];
	gdn_real wI[nb_w];
	gdn_real vn[3];

	for (int id_face = 0; id_face < mesh->nboundaryfaces; id_face++) {
		/* Indirection using boundaryface[id_face] in order to get the 
     * global face id and compute the outgoing normal vector vn 
     */
		int id_face_glob = mesh->boundaryface[id_face];
		mesh_get_face_vn(mesh, id_face_glob, vn);
		assert(mesh_get_bool_is_boundary_face(mesh, id_face_glob));

		/* Loop over the six node of the T10 face */
		for (int id_loc = 0; id_loc < NQV; id_loc++) {
			/* Gauss point in domain */
			int id_pgL =
				mesh_get_ipg_elem_left_from_face(mesh, id_face_glob, id_loc);
			/* Gauss point on boundary */
			int id_pgR = offset + id_face * NQV + id_loc;

			int id_nodeL = mesh->elem2qnode[id_pgL];
			gdn_real *xL = mesh->node + 3 * id_nodeL;

			if (simu->use_kinetic_scheme) {
				model->get_imposed_data(xL, t, wI);
				simulation_get_f(simu, id_pgL, fL);

				model->get_fR(model, fL, vn, wI, fR);
				simulation_set_f(simu, id_pgR, fR);
			} else {
				gdn_real *w = simu->wn + id_pgR * m;
				model->get_imposed_data(xL, t, w);
			}
		}
	}
}

void dg_apply_relaxation(gdn_simulation *simu)
{
	const bool relax_boundary = true; /* USELESS ? */
	const gdn_mesh *mesh = simu->tt;
	const gdn_model *model = simu->mdl;
	const int m = simu->m;
	const int nb_w = model->nb_w;

	/* Warning allocation on the stack using VLA */
	gdn_real f_old[m];
	gdn_real f_new[m];
	gdn_real w_loc[nb_w];

	/* Relax on domain */
	// #pragma omp parallel for
	for (int id_elem = 0; id_elem < mesh->nbelems; id_elem++) {
		for (int id_loc = 0; id_loc < NQN; id_loc++) {
			int id_pg = id_elem * NQN + id_loc;
			simulation_get_f(simu, id_pg, f_old);
			model->get_w(model, f_old, w_loc);
			model->get_feq(model, w_loc, f_new);
			simulation_relax_f(simu, id_pg, f_old, f_new);
		}
	}
	/* Relax on boundaries */
	if (relax_boundary) {
		const int offset = mesh->nbelems * NQN;
		// #pragma omp parallel for
		for (int id_face = 0; id_face < mesh->nboundaryfaces; id_face++) {
			for (int id_loc = 0; id_loc < NQV; id_loc++) {
				int id_pg = offset + id_face * NQV + id_loc;
				simulation_get_f(simu, id_pg, f_old);
				model->get_w(model, f_old, w_loc);
				model->get_feq(model, w_loc, f_new);
				simulation_relax_f(simu, id_pg, f_old, f_new);
			}
		}
	}
}