/*
 * Copyright (c) 
 * 2019-2021 Philippe Helluy <helluy@math.unistra.fr>
 * 2020-2021 Pierre Gerhard <pierre.gerhard@gmail.com>
 *
 * gdon3d is free software; you can redistribute it and/or modify
 * it under the terms of the GPLv3 license. See LICENSE for details.
 */
#include "t4.h"
#include "element.h"
#include <gdon3d.h>

static const Element t4_element = { .elem_code = T4_ELEM_CODE,
									.phy_dim = T4_PHY_DIM,
									.face_per_elem = T4_FACE_PER_ELEM,
									.node_per_elem = T4_NODE_PER_ELEM,
									.node_per_face = T4_NODE_PER_FACE,
									.get_ref_phi = t4_get_ref_phi,
									.get_ref_grad_phi = t4_get_ref_grad_phi,
									.get_ref_normal_entry =
										t4_get_ref_normal_entry };

static const gdn_real t4_ref_node[T4_NODE_PER_ELEM][T4_PHY_DIM] = { { 0, 0, 0 },
																	{ 1, 0, 0 },
																	{ 0, 1, 0 },
																	{ 0, 0,
																	  1 } };

static const gdn_real t4_ref_normal[T4_FACE_PER_ELEM][T4_PHY_DIM] = {
	{ -1, 0, 0 },
	{ 0, -1, 0 },
	{ 0, 0, -1 },
	{ 0.5773502691896258420811705, 0.5773502691896258420811705,
	  0.5773502691896258420811705 }
};

static const int t4_ref_face2node[T4_NODE_PER_ELEM][T4_PHY_DIM] = { { 0, 3, 2 },
																	{ 0, 1, 3 },
																	{ 0, 2, 1 },
																	{ 1, 2,
																	  3 } };

static const gdn_real t4_wpg_ref[T4_NODE_PER_ELEM] = {
	0.41666666666666666667e-1, 0.41666666666666666667e-1,
	0.41666666666666666667e-1, 0.41666666666666666667e-1
};

const Element *t4_get_element(void)
{
	return &t4_element;
}

void t4_get_ref_phi(gdn_real x, gdn_real y, gdn_real z, gdn_real *phi)
{
	phi[0] = -x - y - z + 1;
	phi[1] = x;
	phi[2] = y;
	phi[3] = z;
}

void t4_get_ref_grad_phi(gdn_real x, gdn_real y, gdn_real z, gdn_real *grad_phi)
{
	grad_phi[0] = -1;
	grad_phi[1] = -1;
	grad_phi[2] = -1;
	grad_phi[3] = 1;
	grad_phi[4] = 0;
	grad_phi[5] = 0;
	grad_phi[6] = 0;
	grad_phi[7] = 1;
	grad_phi[8] = 0;
	grad_phi[9] = 0;
	grad_phi[10] = 0;
	grad_phi[11] = 1;
}

void t4_get_ref_node_coord(const int i, gdn_real coord[3])
{
	coord[0] = t4_ref_node[i][0];
	coord[1] = t4_ref_node[i][1];
	coord[2] = t4_ref_node[i][2];
}

gdn_real t4_get_ref_node_entry(const int i, const int j)
{
	return t4_ref_node[i][j];
}

gdn_real t4_get_ref_normal_entry(const int i, const int j)
{
	return t4_ref_normal[i][j];
}

int t4_get_ref_face2node_entry(const int i, const int j)
{
	return t4_ref_face2node[i][j];
}

gdn_real t4_get_wpg_ref_entry(const int i)
{
	return t4_wpg_ref[i];
}