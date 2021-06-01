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
#include <stdarg.h>

#include <gdon3d.h>
#include "simulation.h"

#include <maths/sparse.h>
#include <mesh/mesh.h>
#include <mesh/t10.h>
#include <mesh/tetrageometry.h>
#include <models/model.h>

void static get_f_v1(gdn_simulation *simu, const int id_pg, gdn_real *f)
{
	const int m = simu->m;
	const int neq = simu->neq;

	for (int im = 0; im < m; im++) {
		f[im] = simu->wn[im * neq + id_pg];
	}
}

void static set_f_v1(gdn_simulation *simu, const int id_pg, const gdn_real *f)
{
	const int m = simu->m;
	const int neq = simu->neq;

	for (int im = 0; im < m; im++) {
		simu->wn[im * neq + id_pg] = f[im];
	}
}

void static relax_f_v1(gdn_simulation *simu, const int id_pg,
					   const gdn_real *f_old, const gdn_real *f_new)
{
	const gdn_model *model = simu->mdl;
	const int m = simu->m;
	const int neq = simu->neq;
	const gdn_real omega = model->omega;

	for (int im = 0; im < m; im++) {
		simu->wn[im * neq + id_pg] =
			omega * f_new[im] + (1. - omega) * f_old[im];
	}
}

void static get_f_v2(gdn_simulation *simu, const int id_pg, gdn_real *f)
{
	const gdn_model *model = simu->mdl;
	// const int m = simu->m;
	const int neq = simu->neq;
	const int nb_w = model->nb_w;
	const int nb_v = model->nb_v;

	for (int iv = 0; iv < nb_v; iv++) {
		for (int iw = 0; iw < nb_w; iw++) {
			int offset_wn = iv * nb_w * neq + iw * neq;
			f[iw * nb_v + iv] = simu->wn[offset_wn + id_pg];
		}
	}
}

void static set_f_v2(gdn_simulation *simu, const int id_pg, const gdn_real *f)
{
	const gdn_model *model = simu->mdl;
	// const int m = simu->m;
	const int neq = simu->neq;
	const int nb_w = model->nb_w;
	const int nb_v = model->nb_v;

	for (int iv = 0; iv < nb_v; iv++) {
		for (int iw = 0; iw < nb_w; iw++) {
			int offset_wn = iv * nb_w * neq + iw * neq;
			simu->wn[offset_wn + id_pg] = f[iw * nb_v + iv];
		}
	}
}

void static relax_f_v2(gdn_simulation *simu, const int id_pg,
					   const gdn_real *f_old, const gdn_real *f_new)
{
	const gdn_model *model = simu->mdl;
	// const int m = simu->m;
	const int neq = simu->neq;
	const int nb_w = model->nb_w;
	const int nb_v = model->nb_v;
	const gdn_real omega = model->omega;

	for (int iv = 0; iv < nb_v; iv++) {
		for (int iw = 0; iw < nb_w; iw++) {
			int offset_wn = iv * nb_w * neq + iw * neq;
			simu->wn[offset_wn + id_pg] = omega * f_new[iw * nb_v + iv] +
										  (1. - omega) * f_old[iw * nb_v + iv];
		}
	}
}

void inline simulation_get_f(gdn_simulation *simu, const int id_pg, gdn_real *f)
{
	get_f_v1(simu, id_pg, f);
	// get_f_v2(simu, id_pg, f);
}

void inline simulation_set_f(gdn_simulation *simu, const int id_pg,
							 const gdn_real *f)
{
	set_f_v1(simu, id_pg, f);
	// set_f_v2(simu, id_pg, f);
}

void simulation_relax_f(gdn_simulation *simu, const int id_pg,
						const gdn_real *f_old, const gdn_real *f_new)
{
	relax_f_v1(simu, id_pg, f_old, f_new);
	// relax_f_v2(simu, id_pg, f_old, f_new);
}

static void kinetic_num_flux_upwind(gdn_real wL, gdn_real wR, gdn_real *vnorm,
									gdn_real *vi, gdn_real *flux)
{
	gdn_real vn = vi[0] * vnorm[0] + vi[1] * vnorm[1] + vi[2] * vnorm[2];

	gdn_real vnp = vn > 0 ? vn : 0;
	gdn_real vnm = vn - vnp;

	*flux = vnp * wL + vnm * wR;
}

static void kinetic_num_flux_boundary(gdn_real wL, gdn_real wR, gdn_real *vnorm,
									  gdn_real *vi, gdn_real *flux)
{
	kinetic_num_flux_upwind(wL, wR, vnorm, vi, flux);
}

void simulation_set_time_parameters(gdn_simulation *simu, const gdn_real tmax,
									const gdn_real dt, const gdn_real cfl)
{
	const gdn_model *model = simu->mdl;
	simu->tmax = tmax;

	if (cfl != -1) {
		simu->cfl = cfl;
		simu->dt = (simu->hmin * cfl) / model->lambda;
	} else {
		simu->dt = dt;
		if (simu->mdl->use_relax_scheme)
			simu->cfl = (model->lambda * simu->dt) / (simu->hmin);
		else {
			simu->cfl = (model->lambda * simu->dt) / (simu->hmin);
		}
	}

	simu->itermax = (int)floor(tmax / simu->dt);
	simu->tmax = simu->itermax * simu->dt;
}


void simulation_init(gdn_mesh *msh, gdn_model *mdl, gdn_simulation *f)
{
	assert(msh);
	/* TODO : Check model */

	f->tt = msh;
	f->mdl = mdl;
	f->use_kinetic_scheme = mdl->use_relax_scheme;

	if (f->use_kinetic_scheme) {
		f->m = mdl->nb_v * mdl->nb_w;
	} else {
		f->m = mdl->nb_w;
	}

	f->neq = msh->nbelems * NQN + msh->nboundaryfaces * NQV;
	f->wlen = f->m * f->neq;
	f->wn = calloc(f->wlen, sizeof(gdn_real)); /* Solution buffer */
	assert(f->wn);

	f->t = 0;
	f->hmin = mesh_get_t4_hmin(msh);
	f->cpu_running_t = 0;
	f->cpu_assembly_t = 0;
}

void simulation_free(gdn_simulation *f)
{
	f->m = 0;
	f->wlen = 0;
	f->neq = 0;
	f->t = 0;
	f->hmin = 0;

	if (f->tt) {
		gdn_gdn_mesh_free(f->tt);
	}

	f->mdl = NULL;

	if (f->wn)
		free(f->wn);
}

int simulation_get_varindex(int m, int ipg, int iv)
{
	return iv + m * ipg;
}

/* Assembly routine for non kinetic solver */

void simulation_assemble_mass(gdn_simulation *f, const bool is_inverse,
							  gdn_sparse *dg_mass)
{
	const gdn_mesh *msh = f->tt;
	const int nb_w = f->mdl->nb_w;
	const bool use_kinetic_scheme = f->use_kinetic_scheme;

	gdn_real det;
	gdn_real val;

	const int nb_loc_qnode = mesh_get_nb_loc_qnode(msh);
	const int nb_elem = mesh_get_nb_elem(msh);
	/* Loop over elements */
	for (int id_elem = 0; id_elem < nb_elem; id_elem++) {
		mesh_get_elem_det(msh, id_elem, &det);
		/* Loop over local qnodes */
		for (int iloc = 0; iloc < nb_loc_qnode; iloc++) {
			int ipg = mesh_get_gauss_point_id(msh, id_elem, iloc);

			/* Loop over local qnodes */
			for (int jloc = 0; jloc < nb_loc_qnode; jloc++) {
				int jpg = mesh_get_gauss_point_id(msh, id_elem, jloc);

				if (is_inverse) {
					val = t10_get_inv_mass_ref_entry(iloc, jloc) / det;
				} else {
					val = t10_get_mass_ref_entry(iloc, jloc) * det;
				}

				if (use_kinetic_scheme) {
					sparse_coo_insert(dg_mass, ipg, jpg, val);
				} else {
					/* Loop over conservative variables */
					for (int iw = 0; iw < nb_w; iw++) {
						/* Get corresponding location in the unknown vector */
						int i = simulation_get_varindex(nb_w, ipg, iw);
						int j = simulation_get_varindex(nb_w, jpg, iw);
						sparse_coo_insert(dg_mass, i, j, val);
					}
				}
			}
		}
	}
}

void simulation_assemble_internal(gdn_simulation *simu, gdn_sparse *dg_flux[])
{
	gdn_mesh *mesh = simu->tt;
	gdn_model *model = simu->mdl;

	const int nb_v = model->nb_v;
	const int nb_w = model->nb_w;

	/* Only used if we use non kinetic scheme */
	gdn_real w[nb_w];
	gdn_real flux[nb_w];

	if (!(simu->use_kinetic_scheme)) {
		for (int iw = 0; iw < nb_w; iw++) {
			w[nb_w] = 0;
		}
	}

	for (int id_elem = 0; id_elem < mesh->nbelems; id_elem++) {
		gdn_real codtau[3][3];
		mesh_get_elem_codtau(mesh, id_elem, codtau);

		for (int ii = 0; ii < NQN; ii++) {
			int ipg = mesh_get_gauss_point_id(mesh, id_elem, ii);
			for (int jj = 0; jj < NQN; jj++) {
				int jpg = mesh_get_gauss_point_id(mesh, id_elem, jj);

				/* @TODO : Encapsulate c_grad computation */
				gdn_real c_grad[3] = { 0, 0, 0 };
				for (int l = 0; l < 3; l++) {
					for (int k = 0; k < 3; k++) {
						c_grad[l] += codtau[l][k] *
									 t10_get_internal_ref_entry(ii, jj, k);
					}
				}
				if (simu->use_kinetic_scheme) {
					gdn_real flux_transport;
					for (int iv = 0; iv < nb_v; iv++) {
						gdn_real *v = &(model->vi[3 * iv]);

						kinetic_num_flux_upwind(1., 1., c_grad, v,
												&flux_transport);
						int i = simulation_get_varindex(1, ipg, 0);
						int j = simulation_get_varindex(1, jpg, 0);
						sparse_coo_insert(dg_flux[iv], i, j, -flux_transport);
					}
				} else {
					for (int k = 0; k < nb_w; k++) {
						w[k] = 1;
						int i = simulation_get_varindex(nb_w, ipg, k);
						for (int l = 0; l < nb_w; l++) {
							model->get_num_flux(w, w, c_grad, flux);
							int j = simulation_get_varindex(nb_w, jpg, l);
							sparse_coo_insert(dg_flux[0], i, j, -flux[l]);
						}
						w[k] = 0;
					}
				}
			}
		}
	}
}

void simulation_extend_boundaries(gdn_simulation *simu, gdn_sparse *dg_mass,
								  gdn_sparse *dg_flux[])
{
	const gdn_mesh *mesh = simu->tt;
	const gdn_model *model = simu->mdl;
	const int nb_v = model->nb_v;
	const int nb_w = model->nb_w;

	/* Extend flux matrix to known boundary values */
	int start = mesh->nbelems * NQN;
	for (int ifa = 0; ifa < mesh->nboundaryfaces; ifa++) {
		for (int ino = 0; ino < NQV; ino++) {
			int ipg = start + ifa * NQV + ino;
			if (simu->use_kinetic_scheme) {
				int i = simulation_get_varindex(1, ipg, 0);
				sparse_coo_insert(dg_mass, i, i, 1);
				for (int iv = 0; iv < nb_v; iv++) {
					sparse_coo_insert(dg_flux[iv], i, i, 0);
				}
			} else {
				for (int iw = 0; iw < nb_w; iw++) {
					int i = simulation_get_varindex(nb_w, ipg, iw);
					sparse_coo_insert(dg_mass, i, i, 1);
					sparse_coo_insert(dg_flux[0], i, i, 0);
				}
			}
		}
	}
}

void simulation_assemble_flux(gdn_simulation *simu, gdn_sparse *dg_flux)
{
	const gdn_mesh *mesh = simu->tt;
	const gdn_model *model = simu->mdl;
	const int nb_w = model->nb_w;

	int facount = -1;
	int start = mesh->nbelems * NQN;

	gdn_real flux[nb_w];
	gdn_real *wL = calloc(nb_w, sizeof(gdn_real));
	gdn_real *wR = calloc(nb_w, sizeof(gdn_real));

	for (int ifa = 0; ifa < mesh->nbfaces; ifa++) {
		int ieL = mesh->face2elem[E2F * ifa + 0];
		int ifL = mesh->face2elem[E2F * ifa + 1];
		int ieR = mesh->face2elem[E2F * ifa + 2];

		if (ieR < 0)
			facount++;

		gdn_real cr = 1;
		if (ifL == 3)
			cr = sqrt(3);

		gdn_real physnode[NNO][3];

		for (int iloc = 0; iloc < NNO; iloc++) {
			int ino = mesh->elem2qnode[ieL * NQN + iloc];
			physnode[iloc][0] = mesh->node[ino * 3 + 0];
			physnode[iloc][1] = mesh->node[ino * 3 + 1];
			physnode[iloc][2] = mesh->node[ino * 3 + 2];
		}

		gdn_real dtau[3][3], codtau[3][3], xref[3], vn[3];
		gdn_tetra_ref2phy(physnode, xref, NULL, ifL, NULL, dtau, codtau, NULL,
						  vn);

		for (int ii = 0; ii < 3; ii++)
			vn[ii] *= cr;

		for (int k = 0; k < nb_w; k++) {
			for (int l = 0; l < nb_w; l++) {
				wL[l] = 1;
				if (ieR >= 0) {
					model->get_num_flux(wL, wR, vn, flux);
				} else {
					model->get_num_flux_boundary(wL, wR, vn, flux);
				}
				wL[l] = 0;
				for (int ii = 0; ii < NQV; ii++) {
					int ipg = mesh->fnode2enode[ifa * 2 * NQV + ii] + ieL * NQN;

					for (int jj = 0; jj < NQV; jj++) {
						int jpg =
							mesh->fnode2enode[ifa * 2 * NQV + jj] + ieL * NQN;
						int i = simulation_get_varindex(nb_w, ipg, k);
						int j = simulation_get_varindex(nb_w, jpg, l);

						gdn_real val =
							flux[k] * t10_get_face_scale_entry(ii, jj);
						if (val != 0.)
							sparse_coo_insert(dg_flux, i, j, val);
					}
				}

				for (int ii = 0; ii < NQV; ii++) {
					int ipg =
						mesh->fnode2enode[ifa * 2 * NQV + NQV + ii] + ieR * NQN;

					for (int jj = 0; jj < NQV; jj++) {
						int jpg =
							mesh->fnode2enode[ifa * 2 * NQV + jj] + ieL * NQN;
						int i = simulation_get_varindex(nb_w, ipg, k);
						int j = simulation_get_varindex(nb_w, jpg, l);
						gdn_real val =
							flux[k] * t10_get_face_scale_entry(ii, jj);
						if (ieR >= 0 && val != 0.)
							sparse_coo_insert(dg_flux, i, j, -val);
					}
				}
			}
		}

		for (int k = 0; k < nb_w; k++) {
			for (int l = 0; l < nb_w; l++) {
				wR[l] = 1;
				if (ieR >= 0) {
					model->get_num_flux(wL, wR, vn, flux);
				} else {
					model->get_num_flux_boundary(wL, wR, vn, flux);
				}
				wR[l] = 0;
				for (int ii = 0; ii < NQV; ii++) {
					int ipg = mesh->fnode2enode[ifa * 2 * NQV + ii] + ieL * NQN;

					for (int jj = 0; jj < NQV; jj++) {
						int jpg;
						if (ieR >= 0) {
							jpg = mesh->fnode2enode[ifa * 2 * NQV + NQV + jj] +
								  ieR * NQN;
						} else {
							jpg = mesh->fnode2enode[ifa * 2 * NQV + NQV + jj] +
								  start + facount * NQV;
						}

						int i = simulation_get_varindex(nb_w, ipg, k);
						int j = simulation_get_varindex(nb_w, jpg, l);
						gdn_real val =
							flux[k] * t10_get_face_scale_entry(ii, jj);
						if (val != 0.)
							sparse_coo_insert(dg_flux, i, j, val);
					}
				}

				for (int ii = 0; ii < NQV; ii++) {
					int ipg =
						mesh->fnode2enode[ifa * 2 * NQV + NQV + ii] + ieR * NQN;

					for (int jj = 0; jj < NQV; jj++) {
						int jpg = mesh->fnode2enode[ifa * 2 * NQV + NQV + jj] +
								  ieR * NQN;
						int i = simulation_get_varindex(nb_w, ipg, k);
						int j = simulation_get_varindex(nb_w, jpg, l);
						gdn_real val =
							flux[k] * t10_get_face_scale_entry(ii, jj);
						if (ieR >= 0 && val != 0.)
							sparse_coo_insert(dg_flux, i, j, -val);
					}
				}
			}
		}
	}
	free(wL);
	free(wR);
}

void simulation_assemble_flux_kinetic(gdn_simulation *simu,
									  gdn_sparse *dg_flux[])
{
	gdn_mesh *mesh = simu->tt;
	gdn_model *model = simu->mdl;

	int facount = -1;
	int start = mesh->nbelems * NQN;
	gdn_real flux;

	for (int ifa = 0; ifa < mesh->nbfaces; ifa++) {
		int ieL = mesh->face2elem[E2F * ifa + 0];
		int ifL = mesh->face2elem[E2F * ifa + 1];
		int ieR = mesh->face2elem[E2F * ifa + 2];

		if (ieR < 0)
			facount++;

		gdn_real cr = 1;
		if (ifL == 3)
			cr = sqrt(3);

		gdn_real physnode[NNO][3];

		for (int iloc = 0; iloc < NNO; iloc++) {
			int ino = mesh->elem2qnode[ieL * NQN + iloc];
			physnode[iloc][0] = mesh->node[ino * 3 + 0];
			physnode[iloc][1] = mesh->node[ino * 3 + 1];
			physnode[iloc][2] = mesh->node[ino * 3 + 2];
		}

		gdn_real dtau[3][3], codtau[3][3], xref[3], vn[3];
		gdn_tetra_ref2phy(physnode, xref, NULL, ifL, NULL, dtau, codtau, NULL,
						  vn);

		for (int ii = 0; ii < 3; ii++)
			vn[ii] *= cr;

		for (int iv = 0; iv < model->nb_v; iv++) { /*Loop over velocities */
			gdn_real *v = &(model->vi[3 * iv]);

			if (ieR >= 0)
				kinetic_num_flux_upwind(1., 0., vn, v, &flux);
			else
				kinetic_num_flux_boundary(1., 0., vn, v, &flux);

			for (int ii = 0; ii < NQV; ii++) {
				int ipg = mesh->fnode2enode[ifa * 2 * NQV + ii] + ieL * NQN;
				for (int jj = 0; jj < NQV; jj++) {
					int jpg = mesh->fnode2enode[ifa * 2 * NQV + jj] + ieL * NQN;
					int i = simulation_get_varindex(1, ipg, 0);
					int j = simulation_get_varindex(1, jpg, 0);
					gdn_real val = flux * t10_get_face_scale_entry(ii, jj);
					if (val != 0.)
						sparse_coo_insert(dg_flux[iv], i, j, val);
				}
			}

			for (int ii = 0; ii < NQV; ii++) {
				int ipg =
					mesh->fnode2enode[ifa * 2 * NQV + NQV + ii] + ieR * NQN;
				for (int jj = 0; jj < NQV; jj++) {
					int jpg = mesh->fnode2enode[ifa * 2 * NQV + jj] + ieL * NQN;
					int i = simulation_get_varindex(1, ipg, 0);
					int j = simulation_get_varindex(1, jpg, 0);
					gdn_real val = flux * t10_get_face_scale_entry(ii, jj);
					if (ieR >= 0 && val != 0.)
						sparse_coo_insert(dg_flux[iv], i, j, -val);
				}
			}

			if (ieR >= 0) {
				kinetic_num_flux_upwind(0., 1., vn, v, &flux);
			} else {
				kinetic_num_flux_boundary(0., 1., vn, v, &flux);
			}
			for (int ii = 0; ii < NQV; ii++) {
				int ipg = mesh->fnode2enode[ifa * 2 * NQV + ii] + ieL * NQN;
				for (int jj = 0; jj < NQV; jj++) {
					int jpg;
					if (ieR >= 0) {
						jpg = mesh->fnode2enode[ifa * 2 * NQV + NQV + jj] +
							  ieR * NQN;
					} else {
						jpg = mesh->fnode2enode[ifa * 2 * NQV + NQV + jj] +
							  start + facount * NQV;
					}
					int i = simulation_get_varindex(1, ipg, 0);
					int j = simulation_get_varindex(1, jpg, 0);
					gdn_real val = flux * t10_get_face_scale_entry(ii, jj);
					if (val != 0.) {
						sparse_coo_insert(dg_flux[iv], i, j, val);
					}
				}
			}

			for (int ii = 0; ii < NQV; ii++) {
				int ipg =
					mesh->fnode2enode[ifa * 2 * NQV + NQV + ii] + ieR * NQN;
				for (int jj = 0; jj < NQV; jj++) {
					int jpg =
						mesh->fnode2enode[ifa * 2 * NQV + NQV + jj] + ieR * NQN;
					int i = simulation_get_varindex(1, ipg, 0);
					int j = simulation_get_varindex(1, jpg, 0);
					gdn_real val = flux * t10_get_face_scale_entry(ii, jj);
					if (ieR >= 0 && val != 0.)
						sparse_coo_insert(dg_flux[iv], i, j, -val);
				}
			}
		}
	}
}

gdn_real simulation_error_l2(gdn_simulation *simu)
{
	const gdn_mesh *mesh = simu->tt;
	const gdn_model *model = simu->mdl;
	const int m = simu->m;
	const int nb_w = model->nb_w;

	gdn_real f_mic[m];
	gdn_real *w_loc = NULL;
	gdn_real w_ex[nb_w];
	gdn_real *error_w = (gdn_real *)calloc(nb_w, sizeof(gdn_real));

	gdn_real xphy[3];

	if (simu->use_kinetic_scheme) {
		w_loc = (gdn_real *)calloc(nb_w, sizeof(gdn_real));
	}

	for (int ie = 0; ie < mesh->nbelems; ie++) {
		gdn_real det;
		mesh_get_elem_det(mesh, ie, &det);

		for (int iloc = 0; iloc < NQN; iloc++) {
			int ipg = ie * NQN + iloc;
			int ino = mesh->elem2qnode[ipg];
			xphy[0] = mesh->node[ino * 3 + 0];
			xphy[1] = mesh->node[ino * 3 + 1];
			xphy[2] = mesh->node[ino * 3 + 2];

			model->get_imposed_data(xphy, simu->t, w_ex);

			if (simu->use_kinetic_scheme) {
				simulation_get_f(simu, ipg, f_mic);
				model->get_w(model, f_mic, w_loc);
			} else {
				int imem = simulation_get_varindex(m, ipg, 0);
				w_loc = &(simu->wn[imem]);
			}
			/* Compute and sum the difference */
			for (int iw = 0; iw < nb_w; iw++) {
				gdn_real diff = w_loc[iw] - w_ex[iw];
				error_w[iw] += diff * diff * t10_get_wpg_ref_entry(iloc) * det;
			}
		}
	}

	/* Compute the l2 error */
	gdn_real error_l2 = 0;
	for (int iw = 0; iw < nb_w; iw++) {
		error_w[iw] = sqrt(fabs(error_w[iw]));
		error_l2 += error_w[iw] / nb_w;
	}
	free(error_w);

	if (simu->use_kinetic_scheme) {
		if (w_loc) {
			free(w_loc);
		}
	}
	return error_l2;
}

void simulation_dump_info(const int num_run, gdn_simulation *simu,
						  const gdn_real error, const gdn_real order)
{
	if (num_run > 0) {
		fflush(stdout);
		printf(
			"\r%-30s %-12d %-12d %-12d %-6.12f  %-6.12f  %-6.12f  %-6.12f  %-6.12f  %-6.12f  %-6.12f  %-6.12f\n",
			"fooo", simu->tt->nbelems, simu->m * simu->tt->nbelems * NQN,
			simu->itermax, simu->tmax, simu->cfl, simu->dt, simu->hmin,
			simu->cpu_running_t, simu->cpu_assembly_t, error, order);

	} else {
		fflush(stdout);
		printf(
			"\r%-30s %-12s %-12s %-12s %-14s  %-14s  %-14s  %-14s  %-14s  %-14s  %-14s  %-14s\n",
			"filename", "elements", "DOF", "iter", "tmax", "cfl", "dt", "hmin",
			"cpu_run_t", "cpu_asm_t", "error", "order");

		printf(
			"%-30s %-12d %-12d %-12d %-6.12f  %-6.12f  %-6.12f  %-6.12f  %-6.12f  %-6.12f  %-6.12f  %-14s\n",
			"fooo", simu->tt->nbelems, simu->m * simu->tt->nbelems * NQN,
			simu->itermax, simu->tmax, simu->cfl, simu->dt, simu->hmin,
			simu->cpu_running_t, simu->cpu_assembly_t, error, "---");
	}
}