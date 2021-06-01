/*
 * Copyright (c) 
 * 2019-2021 Philippe Helluy <helluy@math.unistra.fr>
 * 2020-2021 Pierre Gerhard <pierre.gerhard@gmail.com>
 *
 * gdon3d is free software; you can redistribute it and/or modify
 * it under the terms of the GPLv3 license. See LICENSE for details.
 */
#ifndef TRANSPORT_H
#define TRANSPORT_H

#include <gdon3d.h>

void transport_imposed_compact_macro(const gdn_real *x, const gdn_real t,
									 gdn_real *w);

void transport_imposed_gaussian_macro(const gdn_real *x, const gdn_real t,
									  gdn_real *w);

void transport_num_flux_upwind(const gdn_real *wL, const gdn_real *wR,
							   const gdn_real *vnorm, gdn_real *flux);
#endif
