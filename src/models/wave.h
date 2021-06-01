/*
 * Copyright (c) 
 * 2019-2021 Philippe Helluy <helluy@math.unistra.fr>
 * 2020-2021 Pierre Gerhard <pierre.gerhard@gmail.com>
 *
 * gdon3d is free software; you can redistribute it and/or modify
 * it under the terms of the GPLv3 license. See LICENSE for details.
 */
#ifndef WAVE_H
#define WAVE_H

#include <gdon3d.h>

void wave_imposed_macro_compact(const gdn_real *x, const gdn_real t,
								gdn_real *w);

void wave_imposed_macro_plane_wave(const gdn_real *x, const gdn_real t,
								   gdn_real *w);

void wave_flux_phy(const gdn_real *w, const gdn_real *vnorm, gdn_real *flux);

void wave_num_flux_upwind(const gdn_real *wL, const gdn_real *wR,
						  const gdn_real *vnorm, gdn_real *flux);

void wave_num_flux_centered(const gdn_real *wL, const gdn_real *wR,
							const gdn_real *vnorm, gdn_real *flux);

void wave_num_flux_rusanov(const gdn_real *wL, const gdn_real *wR,
						   const gdn_real *vnorm, gdn_real *flux);

void wave_get_proj_neg_prod_w(const gdn_real *w, const gdn_real *n,
							  gdn_real *res);

void wave_get_proj_neg_null_prod_w(const gdn_real *w, const gdn_real *n,
								   gdn_real *res);

void wave_get_proj_pos_prod_w(const gdn_real *w, const gdn_real *n,
							  gdn_real *res);

void wave_get_proj_pos_null_prod_w(const gdn_real *w, const gdn_real *n,
								   gdn_real *res);

/* Two-dimensional second order wave equation routines */

void wave_2d_imposed_macro_plane_wave(const gdn_real *x, const gdn_real t,
									  gdn_real *w);

void wave_2d_flux_phy(const gdn_real *w, const gdn_real *vnorm, gdn_real *flux);

void wave_2d_num_flux_upwind(const gdn_real *wL, const gdn_real *wR,
							 const gdn_real *vnorm, gdn_real *flux);

void wave_2d_num_flux_centered(const gdn_real *wL, const gdn_real *wR,
							   const gdn_real *vnorm, gdn_real *flux);

void wave_2d_num_flux_rusanov(const gdn_real *wL, const gdn_real *wR,
							  const gdn_real *vnorm, gdn_real *flux);

void wave_2d_get_proj_neg_prod_w(const gdn_real *w, const gdn_real *n,
								 gdn_real *res);

void wave_2d_get_proj_neg_null_prod_w(const gdn_real *w, const gdn_real *n,
									  gdn_real *res);

void wave_2d_get_proj_pos_prod_w(const gdn_real *w, const gdn_real *n,
								 gdn_real *res);

void wave_2d_get_proj_pos_null_prod_w(const gdn_real *w, const gdn_real *n,
									  gdn_real *res);
#endif