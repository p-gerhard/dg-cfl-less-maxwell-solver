/*
 * Copyright (c) 
 * 2019-2021 Philippe Helluy <helluy@math.unistra.fr>
 * 2020-2021 Pierre Gerhard <pierre.gerhard@gmail.com>
 *
 * gdon3d is free software; you can redistribute it and/or modify
 * it under the terms of the GPLv3 license. See LICENSE for details.
 */
#undef NDEBUG
#include <assert.h>
#include <stdlib.h>

#include "element.h"
#include <gdon3d.h>
#include "t4.h"

void element_ref2phy(Element *elem, gdn_real *(physnode)[3], gdn_real xref[3],
					 gdn_real dphiref[3], int ifa, gdn_real xphy[3],
					 gdn_real dtau[3][3], gdn_real codtau[3][3],
					 gdn_real dphi[3], gdn_real vnds[3])
{
	int nb_nodes = elem->node_per_elem;

	/* Compute the mapping and its Jacobian */
	gdn_real x = xref[0];
	gdn_real y = xref[1];
	gdn_real z = xref[2];

	/* Gradient of the shape functions and value (4th component) of the
   * shape functions. Mapping function for order 1 tetrahedron
   * @TODO : To be improved for T10 tetrahedra
   */

	if (xphy != NULL) {
		gdn_real *phi = (gdn_real *)malloc(nb_nodes * sizeof(gdn_real));
		elem->get_ref_phi(x, y, z, phi);

		for (int k = 0; k < 3; k++) {
			xphy[k] = 0;
			for (int i = 0; i < nb_nodes; i++) {
				xphy[k] += physnode[i][k] * phi[i];
			}
		}
		free(phi);
	}

	if (dtau != NULL) {
		gdn_real *grad_phi =
			(gdn_real *)malloc(3 * nb_nodes * sizeof(gdn_real));
		elem->get_ref_grad_phi(x, y, z, grad_phi);
		for (unsigned int ii = 0; ii < 3; ii++) {
			for (unsigned int jj = 0; jj < 3; jj++) {
				dtau[ii][jj] = 0;
			}
			for (unsigned int i = 0; i < nb_nodes; i++) {
				for (unsigned int jj = 0; jj < 3; jj++) {
					dtau[ii][jj] += physnode[i][ii] * grad_phi[3 * i + jj];
				}
			}
		}
		free(grad_phi);
	}

	if (codtau != NULL) {
		assert(dtau != NULL);
		codtau[0][0] = dtau[1][1] * dtau[2][2] - dtau[1][2] * dtau[2][1];
		codtau[0][1] = -dtau[1][0] * dtau[2][2] + dtau[1][2] * dtau[2][0];
		codtau[0][2] = dtau[1][0] * dtau[2][1] - dtau[1][1] * dtau[2][0];
		codtau[1][0] = -dtau[0][1] * dtau[2][2] + dtau[0][2] * dtau[2][1];
		codtau[1][1] = dtau[0][0] * dtau[2][2] - dtau[0][2] * dtau[2][0];
		codtau[1][2] = -dtau[0][0] * dtau[2][1] + dtau[0][1] * dtau[2][0];
		codtau[2][0] = dtau[0][1] * dtau[1][2] - dtau[0][2] * dtau[1][1];
		codtau[2][1] = -dtau[0][0] * dtau[1][2] + dtau[0][2] * dtau[1][0];
		codtau[2][2] = dtau[0][0] * dtau[1][1] - dtau[0][1] * dtau[1][0];
	}

	if (dphi != NULL) {
		assert(codtau != NULL);
		for (unsigned int ii = 0; ii < 3; ii++) {
			dphi[ii] = 0;
			for (unsigned int jj = 0; jj < 3; jj++) {
				dphi[ii] += codtau[ii][jj] * dphiref[jj];
			}
		}
	}

	if (vnds != NULL) {
		assert(codtau != NULL);
		assert(ifa >= 0);
		if (ifa >= 0) {
			for (unsigned int ii = 0; ii < 3; ii++) {
				vnds[ii] = 0.0;
				for (unsigned int jj = 0; jj < 3; jj++) {
					vnds[ii] +=
						codtau[ii][jj] * elem->get_ref_normal_entry(ifa, jj);
				}
			}
		}
	}
}
