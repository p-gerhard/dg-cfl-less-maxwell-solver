/*
 * Copyright (c) 
 * 2019-2021 Philippe Helluy <helluy@math.unistra.fr>
 * 2020-2021 Pierre Gerhard <pierre.gerhard@gmail.com>
 *
 * gdon3d is free software; you can redistribute it and/or modify
 * it under the terms of the GPLv3 license. See LICENSE for details.
 */
#ifndef ELEMENT_H
#define ELEMENT_H

#include <gdon3d.h>

typedef struct Element {
	int elem_code;
	int phy_dim;
	int face_per_elem;
	int node_per_elem;
	int node_per_face;
	void (*get_ref_phi)(gdn_real x, gdn_real y, gdn_real z, gdn_real *phi);
	void (*get_ref_grad_phi)(gdn_real x, gdn_real y, gdn_real z,
							 gdn_real *grad_phi);
	gdn_real (*get_ref_normal_entry)(int i, int j);
} Element;

#endif