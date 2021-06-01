/*
 * Copyright (c) 
 * 2019-2021 Philippe Helluy <helluy@math.unistra.fr>
 * 2020-2021 Pierre Gerhard <pierre.gerhard@gmail.com>
 *
 * gdon3d is free software; you can redistribute it and/or modify
 * it under the terms of the GPLv3 license. See LICENSE for details.
 */
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

#include <gdon3d.h>
#include "t4.h"
#include "tetrageometry.h"

#define _ITERNEWTON 20
#define _NTHETA 100

void vect_from_pt(const gdn_real p1[3], const gdn_real p2[3], gdn_real v[3])
{
	v[0] = p2[0] - p1[0];
	v[1] = p2[1] - p1[1];
	v[2] = p2[2] - p1[2];
}

gdn_real vect_get_norm(gdn_real v[3])
{
	return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

void vect_normalize(gdn_real v[3])
{
	gdn_real norm = vect_get_norm(v);
	if (norm != 0) {
		gdn_real r = 1.0 / vect_get_norm(v);
		v[0] *= r;
		v[1] *= r;
		v[2] *= r;
	}
}

gdn_real point_get_dist(gdn_real p1[3], gdn_real p2[3])
{
	gdn_real d[3] = { p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2] };

	return vect_get_norm(d);
}

gdn_real geometry_dot_product(const gdn_real v[3], const gdn_real u[3])
{
	return v[0] * u[0] + v[1] * u[1] + v[2] * u[2];
}

gdn_real geometry_get_dist_point_to_face(gdn_real pt[3], gdn_real face[3][3])
{
	/* Normal vector to the face */
	gdn_real nx = (face[1][1] - face[0][1]) * (face[2][2] - face[0][2]) -
				  (face[1][2] - face[0][2]) * (face[2][1] - face[0][1]);

	gdn_real ny = (face[1][2] - face[0][2]) * (face[2][0] - face[0][0]) -
				  (face[1][0] - face[0][0]) * (face[2][2] - face[0][2]);

	gdn_real nz = (face[1][0] - face[0][0]) * (face[2][1] - face[0][1]) -
				  (face[1][1] - face[0][1]) * (face[2][0] - face[0][0]);

	gdn_real vn[3] = { nx, ny, nz };

	/* Vector relating the point to the face */
	gdn_real vp[3] = { pt[0] - face[0][0], pt[1] - face[0][1],
					   pt[2] - face[0][2] };

	gdn_real normvn = 0;

	/* Distance to the face */
	gdn_real fdist = 0;

	for (unsigned int id = 0; id < 3; id++) {
		normvn += vn[id] * vn[id];
		fdist += vp[id] * vn[id];
	}

	normvn = sqrt(normvn);
	fdist /= normvn;
	return fabs(fdist);
}

/* Computes the angle between two faces of a cell
 * a and b are on intersection of the faces, 
 * c is on first face and d on second face
 */
gdn_real geometry_get_angle_between_faces(gdn_real *a, gdn_real *b, gdn_real *c,
										  gdn_real *d)
{
	gdn_real ab[3], ac[3], ad[3];
	vect_from_pt(a, b, ab);
	vect_from_pt(a, c, ac);
	vect_from_pt(a, d, ad);

	gdn_real alpha1, alpha2, beta;
	alpha1 = geometry_dot_product(ab, ac);
	alpha2 = geometry_dot_product(ab, ad);
	beta = geometry_dot_product(ab, ab);

	gdn_real u[3], v[3];
	for (unsigned int ii = 0; ii < 3; ii++) {
		u[ii] = alpha1 * ab[ii] - beta * ac[ii];
		v[ii] = alpha2 * ab[ii] - beta * ad[ii];
	}

	gdn_real u_dot_v, u_dot_u, v_dot_v;

	u_dot_v = geometry_dot_product(u, v);
	u_dot_u = geometry_dot_product(u, u);
	v_dot_v = geometry_dot_product(v, v);

	gdn_real theta = acos(u_dot_v / (sqrt(u_dot_u) * sqrt(v_dot_v)));
	return theta;
}

void geometry_get_tetra_area_and_vol(gdn_real (*node)[3], gdn_real *surf,
									 gdn_real *vol)
{
	/* Source https://keisan.casio.com/exec/system/1329962711 */
	gdn_real a[6] = {
		point_get_dist(node[3], node[0]), point_get_dist(node[3], node[1]),
		point_get_dist(node[3], node[2]), point_get_dist(node[0], node[1]),
		point_get_dist(node[1], node[2]), point_get_dist(node[0], node[2])
	};

	/* Compute tetrahedron surface */
	gdn_real s1 = (a[0] + a[1] + a[3]) / 2;
	gdn_real s2 = (a[1] + a[2] + a[4]) / 2;
	gdn_real s3 = (a[2] + a[5] + a[0]) / 2;
	gdn_real s4 = (a[3] + a[4] + a[5]) / 2;

	*surf = sqrt(s1 * (s1 - a[0]) * (s1 - a[1]) * (s1 - a[3])) +
			sqrt(s2 * (s2 - a[1]) * (s2 - a[2]) * (s2 - a[4])) +
			sqrt(s3 * (s3 - a[2]) * (s3 - a[5]) * (s3 - a[0])) +
			sqrt(s4 * (s4 - a[3]) * (s4 - a[4]) * (s4 - a[5]));

	/* Compute tetrahedron volume */
	gdn_real a1 = a[0] * a[0];
	gdn_real a2 = a[1] * a[1];
	gdn_real a3 = a[2] * a[2];
	gdn_real a4 = a[3] * a[3];
	gdn_real a5 = a[4] * a[4];
	gdn_real a6 = a[5] * a[5];

	gdn_real v2 = a1 * a5 * (a2 + a3 + a4 + a6 - a1 - a5) +
				  a2 * a6 * (a1 + a3 + a4 + a5 - a2 - a6) +
				  a3 * a4 * (a1 + a2 + a5 + a6 - a3 - a4) - a1 * a2 * a4 -
				  a2 * a3 * a5 - a1 * a3 * a6 - a4 * a5 * a6;
	*vol = sqrt(v2 / 144.);

	/* Compute minimale distance between two nodes of a given tetrahedron */
	// gdn_real hmin = a[0];
	// for (int k = 1; k < 6; k++) {
	//   gdn_real val = a[k];
	//   if (val < hmin)
	//     hmin = val;
	// }
}

void gdn_tetra_ref2phy(gdn_real physnode[4][3], gdn_real xref[3],
					   gdn_real dphiref[3], int ifa, gdn_real xphy[3],
					   gdn_real dtau[3][3], gdn_real codtau[3][3],
					   gdn_real dphi[3], gdn_real vnds[3])
{
	/* @TODO : make it generic for element */
	int nno = T4_NODE_PER_ELEM;

	/* Compute the mapping and its Jacobian */
	gdn_real x = xref[0];
	gdn_real y = xref[1];
	gdn_real z = xref[2];

	/* Gradient of the shape functions and value (4th component) of the
   * shape functions. Mapping function for order 1 tetrahedron
   * @TODO : To be improved for T10 tetrahedra 
   */

	gdn_real gradphi[4][4];
	gradphi[0][0] = -1; /* Contribution of node 0 in x direction */
	gradphi[0][1] = -1; /* Contribution of node 0 in y direction */
	gradphi[0][2] = -1; /* Contribution of node 0 in z direction */
	gradphi[0][3] = 1 - x - y - z;
	gradphi[1][0] = 1;
	gradphi[1][1] = 0;
	gradphi[1][2] = 0;
	gradphi[1][3] = x;
	gradphi[2][0] = 0;
	gradphi[2][1] = 1;
	gradphi[2][2] = 0;
	gradphi[2][3] = y;
	gradphi[3][0] = 0;
	gradphi[3][1] = 0;
	gradphi[3][2] = 1;
	gradphi[3][3] = z;

	if (xphy != NULL) {
		for (unsigned int ii = 0; ii < 3; ii++) {
			xphy[ii] = 0;
			for (unsigned int i = 0; i < nno; i++) {
				xphy[ii] += gradphi[i][3] * physnode[i][ii];
			}
		}
	}

	if (dtau != NULL) {
		for (unsigned int ii = 0; ii < 3; ii++) {
			for (unsigned int jj = 0; jj < 3; jj++) {
				dtau[ii][jj] = 0;
			}
			for (unsigned int i = 0; i < nno; i++) {
				for (unsigned int jj = 0; jj < 3; jj++) {
					dtau[ii][jj] += physnode[i][ii] * gradphi[i][jj];
				}
			}
		}
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
						codtau[ii][jj] * t4_get_ref_normal_entry(ifa, jj);
				}
			}
		}
	}
}

void gdn_tetra_phy2ref(gdn_real physnode[4][3], gdn_real xphy[3],
					   gdn_real xref[3])
{
	gdn_real dtau[3][3], codtau[3][3], dxref[3], dxphy[3];
	int ifa = -1;

	xref[0] = 0.5;
	xref[1] = 0.5;
	xref[2] = 0.5;

	gdn_real *codtau0 = codtau[0];
	gdn_real *codtau1 = codtau[1];
	gdn_real *codtau2 = codtau[2];

	for (unsigned int iter = 0; iter < _ITERNEWTON; ++iter) {
		gdn_tetra_ref2phy(physnode, xref, 0, ifa, dxphy, dtau, codtau, 0, 0);
		dxphy[0] -= xphy[0];
		dxphy[1] -= xphy[1];
		dxphy[2] -= xphy[2];
		gdn_real overdet = 1.0 / geometry_dot_product(dtau[0], codtau[0]);
		//assert(overdet > 0);

		for (unsigned int ii = 0; ii < 3; ii++) {
			dxref[ii] = 0;
			dxref[ii] = codtau0[ii] * dxphy[0] + codtau1[ii] * dxphy[1] +
						codtau2[ii] * dxphy[2];

			xref[ii] -= dxref[ii] * overdet;
		}
	}
}

void gdn_tetra_robust_phy2ref(gdn_real physnode[4][3], gdn_real xphy[3],
							  gdn_real xref[3])
{
	gdn_real dtau[3][3], codtau[3][3], dxref[3], dxphy[3], xphy0[3];
	int ifa = -1;

	/* Construct a point xphy0 for which we know the inverse map */
	xref[0] = 0.5;
	xref[1] = 0.5;
	xref[2] = 0.5;

	gdn_tetra_ref2phy(physnode, xref, NULL, ifa, xphy0, NULL, NULL, NULL, NULL);

	/* Homotopy path :
   * theta = 0 -> xphy0
   * theta = 1 -> xphy
   */

	gdn_real dtheta = 1. / _NTHETA;

	for (unsigned int itheta = 0; itheta <= _NTHETA; itheta++) {
		gdn_real theta = itheta * dtheta;
		/* Intermediate point to find */
		gdn_real xphy1[3];
		for (unsigned int ii = 0; ii < 3; ii++) {
			xphy1[ii] = theta * xphy[ii] + (1 - theta) * xphy0[ii];
		}

		for (unsigned int iter = 0; iter < _ITERNEWTON; ++iter) {
			gdn_tetra_ref2phy(physnode, xref, 0, ifa, dxphy, dtau, codtau, 0,
							  0);
			dxphy[0] -= xphy1[0];
			dxphy[1] -= xphy1[1];
			dxphy[2] -= xphy1[2];
			gdn_real det = geometry_dot_product(dtau[0], codtau[0]);
			//assert(det > 0);

			for (unsigned int ii = 0; ii < 3; ii++) {
				dxref[ii] = 0;
				for (unsigned int jj = 0; jj < 3; jj++) {
					dxref[ii] += codtau[jj][ii] * dxphy[jj];
				}
				xref[ii] -= dxref[ii] / det;
			}
		}

		bool is_in_elem =
			(xref[0] >= _GDN_XREF_LOWBOUND) && (xref[0] <= _GDN_XREF_UPBOUND) &&
			(xref[1] >= _GDN_XREF_LOWBOUND) && (xref[1] <= _GDN_XREF_UPBOUND) &&
			(xref[2] >= _GDN_XREF_LOWBOUND) && (xref[2] <= _GDN_XREF_UPBOUND);

		if (!is_in_elem)
			return;
	}
}