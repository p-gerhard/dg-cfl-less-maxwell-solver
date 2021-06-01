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
#include <math.h>

#include <gdon3d.h>
#include "transport.h"

/* Advection vector */
static const gdn_real V_TRANSPORT[3] = { 1., 0, 0 };

/* Constants for compact and gaussian functions */
static const gdn_real G_ORDER = 5.;
static const gdn_real G_SCALE = 0.20;
static const gdn_real G_CENTER[3] = { 0.3, 0.5, 0.5 };

/**
 * \fn transport_imposed_compact_macro
 * \brief Imposed compact function of diameter G_SCALE and center G_CENTER.
 *        This solution is constructed using the method of characteristics, 
 *        therby it is an exact solution for all (x, t) of the transport 
 *        equation.
 * \param[in] x    Position vector
 * \param[in] t    Time
 * \param[inout] w Macro state
 */
void transport_imposed_compact_macro(const gdn_real *x, const gdn_real t,
									 gdn_real *w)
{
	const gdn_real X = x[0] - G_CENTER[0] - t * V_TRANSPORT[0];
	const gdn_real Y = x[1] - G_CENTER[1] - t * V_TRANSPORT[1];
	const gdn_real Z = x[2] - G_CENTER[2] - t * V_TRANSPORT[2];
	const gdn_real r2 = G_SCALE * G_SCALE;
	const gdn_real d = (X * X + Y * Y + Z * Z) / r2;

	if (d >= 1) {
		w[0] = 0;
	} else {
		w[0] = pow(1. - d, G_ORDER);
	}
}

/**
 * \fn transport_imposed_gaussian_macro
 * \brief Imposed 3D Gaussian function of diameter G_SCALE and center G_CENTER.
 *        This solution is constructed using the method of characteristics, 
 *        therby it is an exact solution for all (x, t) of the transport 
 *        equation.
 * \param[in] x    Position vector
 * \param[in] t    Time
 * \param[inout] w Macro state
 */
void transport_imposed_gaussian_macro(const gdn_real *x, const gdn_real t,
									  gdn_real *w)
{
	const gdn_real X = (x[0] - G_CENTER[0] - t * V_TRANSPORT[0]) / G_SCALE;
	const gdn_real Y = (x[1] - G_CENTER[1] - t * V_TRANSPORT[1]) / G_SCALE;
	const gdn_real Z = (x[2] - G_CENTER[2] - t * V_TRANSPORT[2]) / G_SCALE;

	const gdn_real A = 1.0 / (15.7496099458 * G_SCALE * G_SCALE * G_SCALE);

	w[0] = A * exp(-0.5 * (X * X + Y * Y + Z * Z));
}

/**
 * \fn transport_num_flux_upwind
 * \brief Upwind numerical flux for transport equation.
 * \param[in] wL      Left state
 * \param[in] wR      Right state
 * \param[in] vnorm   Outgoing normal vector
 * \param[inout] flux Numerical flux
 */
void transport_num_flux_upwind(const gdn_real *wL, const gdn_real *wR,
							   const gdn_real *vnorm, gdn_real *flux)
{
	const gdn_real vn = V_TRANSPORT[0] * vnorm[0] + V_TRANSPORT[1] * vnorm[1] +
						V_TRANSPORT[2] * vnorm[2];

	const gdn_real vnp = vn > 0 ? vn : 0;
	const gdn_real vnm = vn - vnp;

	flux[0] = vnp * wL[0] + vnm * wR[0];
}