/*
 * Copyright (c) 
 * 2019-2021 Philippe Helluy <helluy@math.unistra.fr>
 * 2020-2021 Pierre Gerhard <pierre.gerhard@gmail.com>
 *
 * gdon3d is free software; you can redistribute it and/or modify
 * it under the terms of the GPLv3 license. See LICENSE for details.
 */
#ifndef DG_H
#define DG_H

#include <simulation/simulation.h>

void dg_solve(gdn_simulation *simu);

void dg_init_sol_macro(gdn_simulation *simu);

void dg_update_boundary(gdn_simulation *simu, const gdn_real t);

void dg_apply_relaxation(gdn_simulation *simu);

#endif