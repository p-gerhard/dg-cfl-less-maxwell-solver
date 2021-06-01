/*
 * Copyright (c) 
 * 2019-2021 Philippe Helluy <helluy@math.unistra.fr>
 * 2020-2021 Pierre Gerhard <pierre.gerhard@gmail.com>
 *
 * gdon3d is free software; you can redistribute it and/or modify
 * it under the terms of the GPLv3 license. See LICENSE for details.
 */
#ifndef UTILS_CONVERGENCE_H
#define UTILS_CONVERGENCE_H

void utils_l2_time_convergence(gdn_model *model, const int nb_run,
							   const char *msh_filename, const gdn_real init_dt,
							   const gdn_real tmax, const bool export_xdmf);


#endif