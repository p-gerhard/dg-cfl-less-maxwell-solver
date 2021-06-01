/*
 * Copyright (c) 
 * 2019-2021 Philippe Helluy <helluy@math.unistra.fr>
 * 2020-2021 Pierre Gerhard <pierre.gerhard@gmail.com>
 *
 * gdon3d is free software; you can redistribute it and/or modify
 * it under the terms of the GPLv3 license. See LICENSE for details.
 */
#ifndef T4_H
#define T4_H

#include <gdon3d.h>
#include "element.h"

#define T4_ELEM_CODE 4
#define T4_PHY_DIM 3
#define T4_FACE_PER_ELEM 4
#define T4_NODE_PER_ELEM 4
#define T4_NODE_PER_FACE 3

/*
 * Nodes (T4)  : - 0, 1, 2, 3
 * 
 * Edges       : - 0, 1
 *               - 0, 2
 *               - 0, 3
 *               - 1, 3
 *               - 1, 2
 *               - 2, 3
 *               
 * Faces       : - F0 : 0, 3, 2 (Back)
 *               - F1 : 0, 1, 3 (Bottom)
 *               - F2 : 0, 2, 1 (Left)
 *               - F3 : 1, 2, 3 (Right)
 *
 *          2                                 
 *          |\__
 *          |\  \___
 *          | \     \__
 *          |  \       \__
 *          |   \         \__
 *          |     \          \__
 *          |      \            \_
 *          |       \         ___/ \ 3
 *          |        \   ____/     |
 *          |       ___\/          |
 *          |  ____/    \          |
 *          |_/          \         |
 *         0 \__          \        |
 *              \__         \      |
 *                 \__       \     |
 *                    \__     \    |
 *                       \__   \   |
 *                          \__ \  |
 *                             \_\ |
 *                                \|
 *                                  1
 *          
 */

const Element *t4_get_element(void);

void t4_get_ref_phi(gdn_real x, gdn_real y, gdn_real z, gdn_real *phi);
void t4_get_ref_grad_phi(gdn_real x, gdn_real y, gdn_real z,
						 gdn_real *grad_phi);

void t4_get_ref_node_coord(const int i, gdn_real coord[3]);
gdn_real t4_get_ref_node_entry(const int i, const int j);
gdn_real t4_get_ref_normal_entry(const int i, const int j);
int t4_get_ref_face2node_entry(const int i, const int j);
gdn_real t4_get_wpg_ref_entry(const int i);
#endif