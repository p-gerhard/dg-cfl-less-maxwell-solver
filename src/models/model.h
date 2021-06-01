/*
 * Copyright (c) 
 * 2019-2021 Philippe Helluy <helluy@math.unistra.fr>
 * 2020-2021 Pierre Gerhard <pierre.gerhard@gmail.com>
 *
 * gdon3d is free software; you can redistribute it and/or modify
 * it under the terms of the GPLv3 license. See LICENSE for details.
 */
#ifndef MODEL_H
#define MODEL_H

#include <stdbool.h>
#include <gdon3d.h>

/* Typdef for pointers of function in gdn_model struct */
typedef void (*gdn_get_num_flux)(const gdn_real *wL, const gdn_real *wR,
								 const gdn_real *vnorm, gdn_real *flux);

typedef void (*gdn_get_num_flux_boundary)(const gdn_real *wL,
										  const gdn_real *wR,
										  const gdn_real *vnorm,
										  gdn_real *flux);

typedef void (*gdn_get_imposed_data)(const gdn_real *x, const gdn_real t,
									 gdn_real *w);

typedef void (*gdn_get_proj_pos)(const gdn_real *w, const gdn_real *n,
								 gdn_real *res);

typedef void (*gdn_get_proj_neg)(const gdn_real *w, const gdn_real *n,
								 gdn_real *res);

typedef void (*gdn_get_w)(const void *self, const gdn_real *f, gdn_real *w);

typedef void (*gdn_get_feq)(const void *self, const gdn_real *w, gdn_real *feq);

typedef void (*gdn_get_fR)(const void *self, const gdn_real *fL,
						   const gdn_real *vn, const gdn_real *wI,
						   gdn_real *fR);

typedef struct gdn_model {
	char *name;
	int nb_v;
	int nb_w;
	gdn_real *vi;
	gdn_real lambda;

	/* Macro scale functions */
	gdn_get_num_flux get_num_flux;
	gdn_get_num_flux_boundary get_num_flux_boundary;
	gdn_get_imposed_data get_imposed_data;
	gdn_get_proj_pos get_proj_pos;
	gdn_get_proj_neg get_proj_neg;

	/* Micro scale functions */
	bool use_relax_scheme;
	gdn_real omega;
	gdn_get_w get_w;
	gdn_get_feq get_feq;
	gdn_get_fR get_fR;

} gdn_model;

extern void (*gdn_flux_phy_macro)(const gdn_real *wL, const gdn_real *wR,
								  const gdn_real *vnorm, gdn_real *flux);

extern void (*gdn_boundary_flux)(const gdn_real *wL, const gdn_real *wR,
								 const gdn_real *vnorm, gdn_real *flux);

extern void (*gdn_imposed_data)(const gdn_real *x, const gdn_real t,
								gdn_real *w);

extern void (*gdn_imposed_data_micro)(const gdn_real *x, const gdn_real *v,
									  const gdn_real t, gdn_real *w);

#endif