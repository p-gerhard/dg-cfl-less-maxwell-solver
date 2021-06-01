/*
 * Copyright (c) 
 * 2019-2021 Philippe Helluy <helluy@math.unistra.fr>
 * 2020-2021 Pierre Gerhard <pierre.gerhard@gmail.com>
 *
 * gdon3d is free software; you can redistribute it and/or modify
 * it under the terms of the GPLv3 license. See LICENSE for details.
 */
#ifndef TETRAGEOMETRY_H
#define TETRAGEOMETRY_H

#include <gdon3d.h>

gdn_real geometry_dot_product(const gdn_real u[3], const gdn_real v[3]);

void geometry_get_tetra_area_and_vol(gdn_real (*node)[3], gdn_real *surf, 
        gdn_real *vol);

void gdn_tetra_ref2phy(gdn_real physnode[4][3],
             gdn_real xref[3],
             gdn_real dphiref[3],
             int ifa,
             gdn_real xphy[3],
             gdn_real dtau[3][3],
             gdn_real codtau[3][3],
             gdn_real dphi[3],
             gdn_real vn[3]);


void gdn_tetra_phy2ref(gdn_real physnode[4][3],gdn_real xphy[3],gdn_real xref[3]);

void gdn_tetra_robust_phy2ref(gdn_real physnode[4][3],gdn_real xphy[3],gdn_real xref[3]);

#endif
