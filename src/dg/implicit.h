/*
 * Copyright (c) 
 * 2019-2021 Philippe Helluy <helluy@math.unistra.fr>
 * 2020-2021 Pierre Gerhard <pierre.gerhard@gmail.com>
 *
 * gdon3d is free software; you can redistribute it and/or modify
 * it under the terms of the GPLv3 license. See LICENSE for details.
 */
#ifndef IMPLICIT_H
#define IMPLICIT_H

#include <simulation/simulation.h>

void dg_implicit_solve_theta_scheme_relax_v1(gdn_simulation *simu);

void dg_implicit_solve_theta_scheme_relax_v2(gdn_simulation *simu);

#endif
