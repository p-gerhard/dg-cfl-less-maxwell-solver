/*
 * Copyright (c) 
 * 2019-2021 Philippe Helluy <helluy@math.unistra.fr>
 * 2020-2021 Pierre Gerhard <pierre.gerhard@gmail.com>
 *
 * gdon3d is free software; you can redistribute it and/or modify
 * it under the terms of the GPLv3 license. See LICENSE for details.
 */
#ifndef MAXWELL_H
#define MAXWELL_H

#include <gdon3d.h>

void maxwell_imposed_stationary(const gdn_real *x, const gdn_real t,
								gdn_real *w);

void maxwell_imposed_macro_cos(const gdn_real *x, const gdn_real t,
							   gdn_real *w);

void maxwell_imposed_macro_cavity(const gdn_real *x, const gdn_real t,
								  gdn_real *w);

void maxwell_num_flux_upwind(const gdn_real *wL, const gdn_real *wR,
							 const gdn_real *vnorm, gdn_real *flux);

void maxwell_num_flux_boundary_metal(const gdn_real *wL, const gdn_real *wR,
									 const gdn_real *vnorm, gdn_real *flux);

void maxwell_get_proj_neg_prod_w(const gdn_real *w, const gdn_real *n,
								 gdn_real *res);

void maxwell_get_proj_neg_null_prod_w(const gdn_real *w, const gdn_real *n,
									  gdn_real *res);

void maxwell_get_proj_pos_prod_w(const gdn_real *w, const gdn_real *n,
								 gdn_real *res);

void maxwell_get_proj_pos_null_prod_w(const gdn_real *w, const gdn_real *n,
									  gdn_real *res);
#endif