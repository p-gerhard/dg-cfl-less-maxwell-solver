/*
 * Copyright (c) 
 * 2019-2021 Philippe Helluy <helluy@math.unistra.fr>
 * 2020-2021 Pierre Gerhard <pierre.gerhard@gmail.com>
 *
 * gdon3d is free software; you can redistribute it and/or modify
 * it under the terms of the GPLv3 license. See LICENSE for details.
 */
#ifndef T10_H
#define T10_H

#include <gdon3d.h>
#include "element.h"

#define T10_ELEM_CODE 10
#define T10_PHY_DIM 3
#define T10_FACE_PER_ELEM 4
#define T10_NODE_PER_ELEM 10
#define T10_NODE_PER_FACE 6

/*
 * Nodes (T4)  : - 0, 1, 2, 3
 * Nodes (T10) : - 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
 * 
 * Edges       : - 0, 4, 1
 *               - 0, 6, 2
 *               - 0, 7, 3
 *               - 1, 8, 3
 *               - 1, 5, 2
 *               - 2, 9, 3
 *               
 * Faces       : - F0 : 0, 3, 2, 7, 9, 6 (Back)
 *               - F1 : 0, 1, 3, 4, 8, 7 (Bottom)
 *               - F2 : 0, 2, 1, 6, 5, 4 (Left)
 *               - F3 : 1, 2, 3, 5, 9, 8 (Right)
 *
 *          2                                 
 *          |\__
 *          |\  \___
 *          | \     \__
 *          |  \       \__9
 *          |   \         \__
 *          |     \          \__
 *         6|      \            \_
 *          |       \         ___/ \ 3
 *          |        \ 7 ____/     |
 *          |       ___\/          |
 *          |  ____/    \5         |
 *          |_/          \         |
 *         0 \__          \        |
 *              \__         \      |8
 *                 \__       \     |
 *                  4 \__     \    |
 *                       \__   \   |
 *                          \__ \  |
 *                             \_\ |
 *                                \|
 *                                  1
 *          
 */

#define ECODE 4 /* Elements type code : 4 - tetrahedra (T4). */
#define NFA 4 /* Number of element's faces.                */
#define NNO 4 /* Number of element's (T4) nodes.           */
#define NVE 3 /* Number of vertices per face (T4).         */
#define E2F 4 /* umbering of elements to face.             */
#define NQN 10 /* Number of element's (T10) nodes.          */
#define NQV 6 /* Number of nodes per face (T10).           */

const Element *t10_get_element(void);

void t10_get_ref_phi(gdn_real x, gdn_real y, gdn_real z, gdn_real *phi);
void t10_get_ref_grad_phi(gdn_real x, gdn_real y, gdn_real z,
						  gdn_real *grad_phi);
gdn_real t10_get_ref_normal_entry(const int i, const int j);

gdn_real t10_get_ref_node_entry(const int i, const int j);
void t10_get_ref_node_coord(const int i, gdn_real coord[3]);

int t10_get_ref_face2node_entry(const int i, const int j);
gdn_real t10_get_wpg_ref_entry(const int i);
gdn_real t10_get_mass_ref_entry(const int i, const int j);
gdn_real t10_get_inv_mass_ref_entry(const int i, const int j);
gdn_real t10_get_internal_ref_entry(const int i, const int j, const int id_dir);
gdn_real t10_get_flux_ref_entry(const int i, const int j, const int id_face);
gdn_real t10_get_face_scale_entry(const int i, const int j);
#endif