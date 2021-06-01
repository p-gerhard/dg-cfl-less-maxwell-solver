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

#include <maths/lebedev.h>
#include <gdon3d.h>
#include "wave.h"

/* Constant for 2d and 3d homogene wave equation */
static const gdn_real WAVE_SPEED = 1.;

/* Constant for plane 3d wave solution */
static const gdn_real PLANE_WAVE_AMPLITUDE = 1.;
static const gdn_real PLANE_WAVE_K[3] = { 0.5, 0.5, 0.5 };

/* Constant for plane 2d wave solution */
static const gdn_real TWO_DIM_PLANE_WAVE_AMPLITUDE = 1.;
static const gdn_real TWO_DIM_PLANE_WAVE_K[2] = { 1., 1. };

/* Constants for 3d compact support function g(x,y,z) */
static const gdn_real G_ORDER = 5.;
static const gdn_real G_SCALE = 0.20;
static const gdn_real G_CENTER[3] = { 0.5, 0.5, 0.5 };

/**
 * \fn get_g_data
 * \brief Helper function that compute the value of g(x, y z) and it's 
 *        spatial gradient.
 * \param[in] x         Position vector
 * \param[inout] g      Value of function g
 * \param[inout] grad_g Gradient (spatial) vector of g
 */
static void get_g_data(const gdn_real *x, gdn_real *g, gdn_real grad_g[3])
{
	const gdn_real X = x[0] - G_CENTER[0];
	const gdn_real Y = x[1] - G_CENTER[1];
	const gdn_real Z = x[2] - G_CENTER[2];
	const gdn_real r2 = G_SCALE * G_SCALE;
	const gdn_real d = (X * X + Y * Y + Z * Z) / r2;

	if (d >= 1) {
		*g = 0;
		grad_g[0] = 0;
		grad_g[1] = 0;
		grad_g[2] = 0;
	} else {
		/* Compute g */
		*g = pow(1. - d, G_ORDER);

		/* Compute spatial gradient of g */
		gdn_real t1 = -G_ORDER * pow(1. - d, G_ORDER - 1);
		gdn_real dddx = (2 * x[0] - 2 * G_CENTER[0]) / r2;
		gdn_real dddy = (2 * x[1] - 2 * G_CENTER[1]) / r2;
		gdn_real dddz = (2 * x[2] - 2 * G_CENTER[2]) / r2;

		grad_g[0] = t1 * dddx;
		grad_g[1] = t1 * dddy;
		grad_g[2] = t1 * dddz;
	}
}

/**
 * \fn wave_imposed_macro_compact
 * \brief Initial condition for the conservative form of the wave equation.
 *        We want to solve in 3D (X, t) the 2nd order wave equation :
 * 
 *               u_{tt} - c^2 u_{XX}   = 0,
 *                           u(X, t=0) = 0,            (1)
 *                         u_t(X, t=0) = g(x),
 * 
 *        with g smooth (G_ORDER) regular compact (G_SCALE) solution. The exact 
 *        solution of this problem can be computed using "Spherical Means 
 *        Method". The previous system into a system of conservation law using :
 *        the following data layout :
 * 
 *               w = [u_t, u_x, u_y, u_z]^T.
 * 
 *        The current function compute w(X, t) for any X and t as an exact 
 *        solution of the system (1).
 *         
 * \param[in] x    Position vector
 * \param[in] t    Time
 * \param[inout] w Exact solution
 */
void wave_imposed_macro_compact(const gdn_real *x, const gdn_real t,
								gdn_real *w)
{
	gdn_real dudt = 0;
	gdn_real dudx = 0;
	gdn_real dudy = 0;
	gdn_real dudz = 0;
	const gdn_real r = WAVE_SPEED * t;

	/* Spherical (S2) Means Method is done using a using Lebedev quadrature 
   * WARNING : In our case we have :
   *     \sum_{i=0}^{N_LEB - 1} w_i = 1
   * which implies a (4 * Pi) simplification at the end.
   */

	const lebedev *leb_set = &quad_leb_029;

#pragma omp for nowait
	for (int i = 0; i < leb_set->N; i++) {
		gdn_real g;
		gdn_real grad_g[3] = { 0, 0, 0 };
		gdn_real y[3] = { 0, 0, 0 };

		gdn_real vx = leb_set->vi[3 * i + 0];
		gdn_real vy = leb_set->vi[3 * i + 1];
		gdn_real vz = leb_set->vi[3 * i + 2];
		gdn_real wi = leb_set->wi[i];

		y[0] = x[0] + r * vx;
		y[1] = x[1] + r * vy;
		y[2] = x[2] + r * vz;

		get_g_data(y, &g, grad_g);

		dudt +=
			wi * (g + r * (grad_g[0] * vx + grad_g[1] * vy + grad_g[2] * vz));
		dudx += wi * grad_g[0];
		dudy += wi * grad_g[1];
		dudz += wi * grad_g[2];
	}

	w[0] = dudt;
	w[1] = t * dudx;
	w[2] = t * dudy;
	w[3] = t * dudz;
}

/**
 * \fn wave_imposed_macro_plane_wave
 * \brief (3D) Plane wave solution of the conservative 2nd order wave equation :
 *               u_{tt} - c^2 u_{XX} = 0,
 *        Plane wave form : 
 *               u(x,y,z,t) = cos(k . x - omega * t + phi)
 *        Data layout :
 *               w = [u_t, u_x, u_y, u_z]^T.
 * \param[in] x    Position vector
 * \param[in] t    Time
 * \param[inout] w Exact solution
 */
void wave_imposed_macro_plane_wave(const gdn_real *x, const gdn_real t,
								   gdn_real *w)
{
	const gdn_real k_dot_x = x[0] * PLANE_WAVE_K[0] + x[1] * PLANE_WAVE_K[1] +
							 x[2] * PLANE_WAVE_K[2];

	gdn_real omega = sqrt(PLANE_WAVE_K[0] * PLANE_WAVE_K[0] +
						  PLANE_WAVE_K[1] * PLANE_WAVE_K[1] +
						  PLANE_WAVE_K[2] * PLANE_WAVE_K[2]) *
					 WAVE_SPEED;

	//  const gdn_real t0 = PLANE_WAVE_AMPLITUDE * sin(k_dot_x - omega * t);
	const gdn_real t0 =
		PLANE_WAVE_AMPLITUDE * (k_dot_x - omega * t) * (k_dot_x - omega * t);

	w[0] = omega * t0;
	w[1] = -PLANE_WAVE_K[0] * t0;
	w[2] = -PLANE_WAVE_K[1] * t0;
	w[3] = -PLANE_WAVE_K[2] * t0;
}

/**
 * \fn wave_flux_phy
 * \brief Compute the physical flux for the conservative form of the 2nd order 
 *        wave equation :
 *               u_{tt} - c^2 u_{xx} = 0.
 * \param[in] w        State vector
 * \param[inout] vnorm Outgoing normal vector
 * \param[inout] flux  Physical flux
 */
void wave_flux_phy(const gdn_real *w, const gdn_real *vnorm, gdn_real *flux)
{
	gdn_real c2 = WAVE_SPEED * WAVE_SPEED;

	flux[0] = -c2 * (vnorm[0] * w[1] + vnorm[1] * w[2] + vnorm[2] * w[3]);
	flux[1] = -w[0] * vnorm[0];
	flux[2] = -w[0] * vnorm[1];
	flux[3] = -w[0] * vnorm[2];
}

/**
 * \fn wave_num_flux_upwind
 * \brief Upwind numerical flux for the conservative form of the 2nd order 
 *        wave equation :
 *               u_{tt} - c^2 u_{xx} = 0.
 * \param[in] wL       Left state vector
 * \param[in] wR       Right state vector
 * \param[inout] vnorm Outgoing normal vector
 * \param[inout] flux  Numerical flux
 */
void wave_num_flux_upwind(const gdn_real *wL, const gdn_real *wR,
						  const gdn_real *vnorm, gdn_real *flux)
{
	gdn_real t1;
	gdn_real t11;
	gdn_real t13;
	gdn_real t14;
	gdn_real t17;
	gdn_real t18;
	gdn_real t2;
	gdn_real t21;
	gdn_real t22;
	gdn_real t30;
	gdn_real t4;
	gdn_real t42;
	gdn_real t45;
	gdn_real t5;
	gdn_real t6;
	gdn_real t7;
	gdn_real t8;
	gdn_real t9;
	t1 = wL[0];
	t2 = wR[0];
	t4 = vnorm[0];
	t5 = t4 * t4;
	t6 = vnorm[1];
	t7 = t6 * t6;
	t8 = vnorm[2];
	t9 = t8 * t8;
	t11 = sqrt(t5 + t7 + t9);
	t13 = wR[1];
	t14 = wL[1];
	t17 = wR[2];
	t18 = wL[2];
	t21 = wR[3];
	t22 = wL[3];
	flux[0] =
		-WAVE_SPEED *
		(t11 * (-t1 + t2) + WAVE_SPEED * (t4 * (t13 + t14) + t6 * (t18 + t17) +
										  (t21 + t22) * t8)) /
		0.2e1;
	t30 = 0.1e1 / t11;
	t42 = t11 * (t1 + t2) +
		  (t4 * (t13 - t14) + t6 * (t17 - t18) + (t21 - t22) * t8) * WAVE_SPEED;
	flux[1] = -t42 * t4 * t30 / 0.2e1;
	t45 = t42 * t30;
	flux[2] = -t6 * t45 / 0.2e1;
	flux[3] = -t8 * t45 / 0.2e1;
}

/**
 * \fn wave_num_flux_centered
 * \brief Centered numerical flux for the conservative form of the 2nd order 
 *        wave equation :
 *               u_{tt} - c^2 u_{xx} = 0.
 * 
 *        Centered numerical flux F* : 
 * 
 *               F*(wL, wR, n) = 0.5 * (F(wL) + F(wR))
 * 
 *        F(.) the physical flux. 
 * 
 * \param[in] wL       Left state vector
 * \param[in] wR       Right state vector
 * \param[inout] vnorm Outgoing normal vector
 * \param[inout] flux  Numerical flux
 */
void wave_num_flux_centered(const gdn_real *wL, const gdn_real *wR,
							const gdn_real *vnorm, gdn_real *flux)
{
	gdn_real t1;
	gdn_real t16;
	gdn_real t2;
	gdn_real t27;
	gdn_real t9;
	t1 = WAVE_SPEED * WAVE_SPEED;
	t2 = vnorm[0];
	t9 = vnorm[1];
	t16 = vnorm[2];
	flux[0] =
		-(0.5e0 * t2 * wL[1] + 0.5e0 * wR[1] * t2 + 0.5e0 * wL[2] * t9 +
		  0.5e0 * wR[2] * t9 + 0.5e0 * wL[3] * t16 + 0.5e0 * wR[3] * t16) *
		t1;

	t27 = wL[0] + wR[0];
	flux[1] = -0.5e0 * t27 * t2;
	flux[2] = -0.5e0 * t27 * t9;
	flux[3] = -0.5e0 * t27 * t16;
}

/**
 * \fn wave_num_flux_rusanov
 * \brief Rusanov numerical flux for the conservative form of the 2nd order 
 *        wave equation :
 *               u_{tt} - c^2 u_{xx} = 0.
 *        Rusanov numerical flux F* : 
 * 
 *               F*(wL, wR, n) = 0.5 *(F(wL) + F(wR)).n - 0.5 * lambda (wR - wL)
 * 
 *        with lambda = max(|Jac(wR)|, |Jac(wL)|) (here = WAVE_SPEED) and
 *        F(.) the physical flux.  
 * \param[in] wL       Left state vector
 * \param[in] wR       Right state vector
 * \param[inout] vnorm Outgoing normal vector
 * \param[inout] flux  Numerical flux
 */
void wave_num_flux_rusanov(const gdn_real *wL, const gdn_real *wR,
						   const gdn_real *vnorm, gdn_real *flux)
{
	flux[0] = 0.5 * ((wL[1] + wR[1]) * vnorm[0] + (wL[2] + wR[2]) * vnorm[1] +
					 (wL[3] + wR[3]) * vnorm[2]);

	gdn_real t0 = 0.5 * (wL[0] + wR[0]);

	flux[1] = t0 * vnorm[0];
	flux[2] = t0 * vnorm[1];
	flux[3] = t0 * vnorm[2];

	flux[0] = WAVE_SPEED * (flux[0] - 0.5 * (wR[0] - wL[0]));
	flux[1] = WAVE_SPEED * (flux[1] - 0.5 * (wR[1] - wL[1]));
	flux[2] = WAVE_SPEED * (flux[2] - 0.5 * (wR[2] - wL[2]));
	flux[3] = WAVE_SPEED * (flux[3] - 0.5 * (wR[3] - wL[3]));
}

/**
 * \fn wave_get_proj_neg_prod_w
 * \brief Return the projection of the state w over the negative space
 *        {-} of A_i.
 *        For more details see Maple worksheet doc/model_wave.mw.
 * \param[in] w      State vector
 * \param[in] n      Outgoing normal vector
 * \param[inout] res Projection of w
 */
void wave_get_proj_neg_prod_w(const gdn_real *w, const gdn_real *n,
							  gdn_real *res)
{
	gdn_real t1;
	gdn_real t10;
	gdn_real t11;
	gdn_real t12;
	gdn_real t15;
	gdn_real t18;
	gdn_real t2;
	gdn_real t23;
	gdn_real t28;
	gdn_real t3;
	gdn_real t31;
	gdn_real t32;
	gdn_real t34;
	gdn_real t35;
	gdn_real t4;
	gdn_real t40;
	gdn_real t44;
	gdn_real t5;
	gdn_real t6;
	gdn_real t7;
	gdn_real t8;
	gdn_real t9;
	t1 = w[0];
	t2 = n[0];
	t3 = t2 * t2;
	t4 = n[1];
	t5 = t4 * t4;
	t6 = n[2];
	t7 = t6 * t6;
	t8 = t3 + t5 + t7;
	t9 = sqrt(t8);
	t10 = 0.1e1 / t9;
	t11 = t10 * WAVE_SPEED;
	t12 = w[1];
	t15 = w[2];
	t18 = w[3];
	res[0] = t12 * t2 * t11 / 0.2e1 + t15 * t4 * t11 / 0.2e1 +
			 t18 * t6 * t11 / 0.2e1 + t1 / 0.2e1;
	t23 = t10 / WAVE_SPEED;
	t28 = 0.1e1 / t8 / 0.2e1;
	t31 = t4 * t2;
	t32 = t15 * t28;
	t34 = t6 * t2;
	t35 = t18 * t28;
	res[1] = t1 * t2 * t23 / 0.2e1 + t12 * t28 * t3 + t32 * t31 + t35 * t34;
	t40 = t12 * t28;
	t44 = t6 * t4;
	res[2] = t1 * t4 * t23 / 0.2e1 + t40 * t31 + t15 * t28 * t5 + t35 * t44;
	res[3] = t1 * t6 * t23 / 0.2e1 + t40 * t34 + t32 * t44 + t18 * t28 * t7;
}

/**
 * \fn wave_get_proj_neg_null_prod_w
 * \brief Return the projection of the state w over the negative and null space
 *        {-, 0} of A_i.
 *        For more details see Maple worksheet doc/model_wave.mw.
 * \param[in] w      State vector
 * \param[in] n      Outgoing normal vector
 * \param[inout] res Projection of w
 */
void wave_get_proj_neg_null_prod_w(const gdn_real *w, const gdn_real *n,
								   gdn_real *res)
{
	gdn_real t1;
	gdn_real t10;
	gdn_real t11;
	gdn_real t12;
	gdn_real t15;
	gdn_real t18;
	gdn_real t2;
	gdn_real t23;
	gdn_real t27;
	gdn_real t28;
	gdn_real t3;
	gdn_real t31;
	gdn_real t34;
	gdn_real t35;
	gdn_real t37;
	gdn_real t38;
	gdn_real t4;
	gdn_real t43;
	gdn_real t45;
	gdn_real t49;
	gdn_real t5;
	gdn_real t6;
	gdn_real t7;
	gdn_real t8;
	gdn_real t9;
	t1 = w[0];
	t2 = n[0];
	t3 = t2 * t2;
	t4 = n[1];
	t5 = t4 * t4;
	t6 = n[2];
	t7 = t6 * t6;
	t8 = t3 + t5 + t7;
	t9 = sqrt(t8);
	t10 = 0.1e1 / t9;
	t11 = t10 * WAVE_SPEED;
	t12 = w[1];
	t15 = w[2];
	t18 = w[3];
	res[0] = t12 * t2 * t11 / 0.2e1 + t15 * t4 * t11 / 0.2e1 +
			 t18 * t6 * t11 / 0.2e1 + t1 / 0.2e1;
	t23 = t10 / WAVE_SPEED;
	t27 = 0.2e1 * t5;
	t28 = 0.2e1 * t7;
	t31 = 0.1e1 / t8 / 0.2e1;
	t34 = t4 * t2;
	t35 = t15 * t31;
	t37 = t6 * t2;
	t38 = t18 * t31;
	res[1] = t1 * t2 * t23 / 0.2e1 + t12 * t31 * (t3 + t27 + t28) - t35 * t34 -
			 t38 * t37;
	t43 = t12 * t31;
	t45 = 0.2e1 * t3;
	t49 = t6 * t4;
	res[2] = t1 * t4 * t23 / 0.2e1 - t43 * t34 + t15 * t31 * (t45 + t5 + t28) -
			 t38 * t49;
	res[3] = t1 * t6 * t23 / 0.2e1 - t43 * t37 - t35 * t49 +
			 t18 * t31 * (t45 + t27 + t7);
}

/**
 * \fn wave_get_proj_prod_prod_w
 * \brief Return the projection of the state w over the positive space
 *        {-} of A_i.
 *        For more details see Maple worksheet doc/model_wave.mw.
 * \param[in] w      State vector
 * \param[in] n      Outgoing normal vector
 * \param[inout] res Projection of w
 */
void wave_get_proj_pos_prod_w(const gdn_real *w, const gdn_real *n,
							  gdn_real *res)
{
	gdn_real t1;
	gdn_real t10;
	gdn_real t11;
	gdn_real t12;
	gdn_real t15;
	gdn_real t18;
	gdn_real t2;
	gdn_real t23;
	gdn_real t28;
	gdn_real t3;
	gdn_real t31;
	gdn_real t32;
	gdn_real t34;
	gdn_real t35;
	gdn_real t4;
	gdn_real t40;
	gdn_real t44;
	gdn_real t5;
	gdn_real t6;
	gdn_real t7;
	gdn_real t8;
	gdn_real t9;
	t1 = w[0];
	t2 = n[0];
	t3 = t2 * t2;
	t4 = n[1];
	t5 = t4 * t4;
	t6 = n[2];
	t7 = t6 * t6;
	t8 = t3 + t5 + t7;
	t9 = sqrt(t8);
	t10 = 0.1e1 / t9;
	t11 = t10 * WAVE_SPEED;
	t12 = w[1];
	t15 = w[2];
	t18 = w[3];
	res[0] = -t12 * t2 * t11 / 0.2e1 - t15 * t4 * t11 / 0.2e1 -
			 t18 * t6 * t11 / 0.2e1 + t1 / 0.2e1;
	t23 = t10 / WAVE_SPEED;
	t28 = 0.1e1 / t8 / 0.2e1;
	t31 = t4 * t2;
	t32 = t15 * t28;
	t34 = t6 * t2;
	t35 = t18 * t28;
	res[1] = -t1 * t2 * t23 / 0.2e1 + t12 * t28 * t3 + t32 * t31 + t35 * t34;
	t40 = t12 * t28;
	t44 = t6 * t4;
	res[2] = -t1 * t4 * t23 / 0.2e1 + t40 * t31 + t15 * t28 * t5 + t35 * t44;
	res[3] = -t1 * t6 * t23 / 0.2e1 + t40 * t34 + t32 * t44 + t18 * t28 * t7;
}

/**
 * \fn wave_get_proj_pos_null_prod_w
 * \brief Return the projection of the state w over the positive and null space
 *        {+, 0} of A_i.
 *        For more details see Maple worksheet doc/model_wave.mw.
 * \param[in] w      State vector
 * \param[in] n      Outgoing normal vector
 * \param[inout] res Projection of w
 */
void wave_get_proj_pos_null_prod_w(const gdn_real *w, const gdn_real *n,
								   gdn_real *res)
{
	gdn_real t1;
	gdn_real t10;
	gdn_real t11;
	gdn_real t12;
	gdn_real t15;
	gdn_real t18;
	gdn_real t2;
	gdn_real t23;
	gdn_real t27;
	gdn_real t28;
	gdn_real t3;
	gdn_real t31;
	gdn_real t34;
	gdn_real t35;
	gdn_real t37;
	gdn_real t38;
	gdn_real t4;
	gdn_real t43;
	gdn_real t45;
	gdn_real t49;
	gdn_real t5;
	gdn_real t6;
	gdn_real t7;
	gdn_real t8;
	gdn_real t9;
	t1 = w[0];
	t2 = n[0];
	t3 = t2 * t2;
	t4 = n[1];
	t5 = t4 * t4;
	t6 = n[2];
	t7 = t6 * t6;
	t8 = t3 + t5 + t7;
	t9 = sqrt(t8);
	t10 = 0.1e1 / t9;
	t11 = t10 * WAVE_SPEED;
	t12 = w[1];
	t15 = w[2];
	t18 = w[3];
	res[0] = -t12 * t2 * t11 / 0.2e1 - t15 * t4 * t11 / 0.2e1 -
			 t18 * t6 * t11 / 0.2e1 + t1 / 0.2e1;
	t23 = t10 / WAVE_SPEED;
	t27 = 0.2e1 * t5;
	t28 = 0.2e1 * t7;
	t31 = 0.1e1 / t8 / 0.2e1;
	t34 = t4 * t2;
	t35 = t15 * t31;
	t37 = t6 * t2;
	t38 = t18 * t31;
	res[1] = -t1 * t2 * t23 / 0.2e1 + t12 * t31 * (t3 + t27 + t28) - t35 * t34 -
			 t38 * t37;
	t43 = t12 * t31;
	t45 = 0.2e1 * t3;
	t49 = t6 * t4;
	res[2] = -t1 * t4 * t23 / 0.2e1 - t43 * t34 + t15 * t31 * (t45 + t5 + t28) -
			 t38 * t49;
	res[3] = -t1 * t6 * t23 / 0.2e1 - t43 * t37 - t35 * t49 +
			 t18 * t31 * (t45 + t27 + t7);
}

/* Two-dimensional second order wave equation routines */

/**
 * \fn wave_2d_imposed_macro_plane_wave
 * \brief (2D) Plane wave solution of the conservative 2nd order wave equation :
 *               u_{tt} - c^2 u_{XX} = 0,
 *        Plane wave form : 
 *               u(x,y,t) = cos(k . x - omega * t + phi)
 *        Data layout :
 *               w = [u_t, u_x, u_y]^T.
 * \param[in] x    Position vector
 * \param[in] t    Time
 * \param[inout] w Exact solution
 */
void wave_2d_imposed_macro_plane_wave(const gdn_real *x, const gdn_real t,
									  gdn_real *w)
{
	const gdn_real k_dot_x =
		x[0] * TWO_DIM_PLANE_WAVE_K[0] + x[1] * TWO_DIM_PLANE_WAVE_K[1];

	const gdn_real omega =
		sqrt(TWO_DIM_PLANE_WAVE_K[0] * TWO_DIM_PLANE_WAVE_K[0] +
			 TWO_DIM_PLANE_WAVE_K[1] * TWO_DIM_PLANE_WAVE_K[1]) *
		WAVE_SPEED;

	const gdn_real t0 = TWO_DIM_PLANE_WAVE_AMPLITUDE * sin(k_dot_x - omega * t);

	w[0] = omega * t0;
	w[1] = -TWO_DIM_PLANE_WAVE_K[0] * t0;
	w[2] = -TWO_DIM_PLANE_WAVE_K[1] * t0;
}

/**
 * \fn wave_2d_flux_phy
 * \brief Compute the physical flux for the conservative form of the 
 *        two-dimensional 2nd order wave equation :
 *               u_{tt} - c^2 u_{xx} = 0.
 * \param[in] w        State vector
 * \param[inout] vnorm Outgoing normal vector
 * \param[inout] flux  Physical flux
 */
void wave_2d_phy_flux(const gdn_real *w, const gdn_real *vnorm, gdn_real *flux)
{
	gdn_real t11;
	gdn_real t2;
	gdn_real t5;
	gdn_real t9;
	t2 = WAVE_SPEED * WAVE_SPEED;
	t5 = vnorm[0];
	t9 = vnorm[1];
	flux[0] = -t2 * t5 * w[1] - t2 * t9 * w[2];
	t11 = w[0];
	flux[1] = -t11 * t5;
	flux[2] = -t11 * t9;
}

/**
 * \fn wave_2d_num_flux_upwind
 * \brief Upwind numerical flux for the conservative form of the two-dimensional 
 *        2nd order wave equation :
 *               u_{tt} - c^2 u_{xx} = 0.
 * \param[in] wL       Left state vector
 * \param[in] wR       Right state vector
 * \param[inout] vnorm Outgoing normal vector
 * \param[inout] flux  Numerical flux
 */
void wave_2d_num_flux_upwind(const gdn_real *wL, const gdn_real *wR,
							 const gdn_real *vnorm, gdn_real *flux)
{
	gdn_real t10;
	gdn_real t12;
	gdn_real t13;
	gdn_real t16;
	gdn_real t17;
	gdn_real t2;
	gdn_real t25;
	gdn_real t3;
	gdn_real t34;
	gdn_real t5;
	gdn_real t6;
	gdn_real t7;
	gdn_real t8;
	t2 = wL[0];
	t3 = wR[0];
	t5 = vnorm[0];
	t6 = t5 * t5;
	t7 = vnorm[1];
	t8 = t7 * t7;
	t10 = sqrt(t6 + t8);
	t12 = wL[1];
	t13 = wR[1];
	t16 = wL[2];
	t17 = wR[2];
	flux[0] = -(t10 * (-t2 + t3) +
				WAVE_SPEED * (t5 * (t12 + t13) + (t16 + t17) * t7)) *
			  WAVE_SPEED / 0.2e1;
	t25 = 0.1e1 / t10;
	t34 = t10 * (-t2 - t3) + WAVE_SPEED * (t5 * (t12 - t13) + (t16 - t17) * t7);
	flux[1] = t5 * t34 * t25 / 0.2e1;
	flux[2] = t34 * t7 * t25 / 0.2e1;
}

/**
 * \fn wave_2d_num_flux_centered
 * \brief Centered numerical flux for the conservative form of the 
 *        two-dimensional 2nd order wave equation :
 *               u_{tt} - c^2 u_{xx} = 0.
 * 
 *        Centered numerical flux F* : 
 * 
  *               F*(wL, wR, n) = 0.5 * (F(wL) + F(wR))
 * 
 *        F(.) the physical flux. 
 * 
 * \param[in] wL       Left state vector
 * \param[in] wR       Right state vector
 * \param[inout] vnorm Outgoing normal vector
 * \param[inout] flux  Numerical flux
 */
void wave_2d_num_flux_centered(const gdn_real *wL, const gdn_real *wR,
							   const gdn_real *vnorm, gdn_real *flux)
{
	gdn_real t10;
	gdn_real t2;
	gdn_real t21;
	gdn_real t3;
	t2 = WAVE_SPEED * WAVE_SPEED;
	t3 = vnorm[0];
	t10 = vnorm[1];
	flux[0] = -(0.5e0 * wL[1] * t3 + 0.5e0 * wR[1] * t3 + 0.5e0 * wL[2] * t10 +
				0.5e0 * wR[2] * t10) *
			  t2;

	t21 = wL[0] + wR[0];
	flux[1] = -0.5e0 * t21 * t3;
	flux[2] = -0.5e0 * t10 * t21;
}

/**
 * \fn wave_num_flux_rusanov
 * \brief Rusanov numerical flux for the conservative form of the 
 *        two-dimensional  2nd order wave equation : 
 *               u_{tt} - c^2 u_{xx} = 0.
 *        Rusanov numerical flux F* : 
 * 
 *               F*(wL, wR, n) = 0.5 *(F(wL) + F(wR)).n - 0.5 * lambda (wR - wL)
 * 
 *        with lambda = max(|Jac(wR)|, |Jac(wL)|) (here = WAVE_SPEED) and
 *        F(.) the physical flux.  
 * \param[in] wL       Left state vector
 * \param[in] wR       Right state vector
 * \param[inout] vnorm Outgoing normal vector
 * \param[inout] flux  Numerical flux
 */
void wave_2d_num_flux_rusanov(const gdn_real *wL, const gdn_real *wR,
							  const gdn_real *vnorm, gdn_real *flux)
{
	gdn_real t10;
	gdn_real t12;
	gdn_real t14;
	gdn_real t17;
	gdn_real t20;
	gdn_real t22;
	gdn_real t25;
	gdn_real t3;
	gdn_real t5;
	gdn_real t7;
	t3 = 0.5e0 * wL[1];
	t5 = 0.5e0 * wR[1];
	t7 = vnorm[0];
	t10 = 0.5e0 * wL[2];
	t12 = 0.5e0 * wR[2];
	t14 = vnorm[1];
	t17 = WAVE_SPEED * WAVE_SPEED;
	t20 = 0.5e0 * wL[0];
	t22 = 0.5e0 * wR[0];
	flux[0] =
		t17 * (t7 * (-t3 - t5) + t14 * (-t10 - t12)) + WAVE_SPEED * (t20 - t22);
	t25 = -t20 - t22;
	flux[1] = t7 * t25 + WAVE_SPEED * (t3 - t5);
	flux[2] = t14 * t25 + WAVE_SPEED * (t10 - t12);
}

/**
 * \fn wave_2d_get_proj_neg_prod_w
 * \brief Return the projection of the state w over the negative space
 *        {-} of A_i.
 *        For more details see Maple worksheet doc/model_wave.mw.
 * \param[in] w      State vector
 * \param[in] n      Outgoing normal vector
 * \param[inout] res Projection of w
 */
void wave_2d_get_proj_neg_prod_w(const gdn_real *w, const gdn_real *n,
								 gdn_real *res)
{
	gdn_real t10;
	gdn_real t11;
	gdn_real t14;
	gdn_real t18;
	gdn_real t2;
	gdn_real t24;
	gdn_real t27;
	gdn_real t3;
	gdn_real t4;
	gdn_real t5;
	gdn_real t6;
	gdn_real t7;
	gdn_real t8;
	gdn_real t9;
	t2 = w[0];
	t3 = n[0];
	t4 = t3 * t3;
	t5 = n[1];
	t6 = t5 * t5;
	t7 = t4 + t6;
	t8 = sqrt(t7);
	t9 = 0.1e1 / t8;
	t10 = WAVE_SPEED * t9;
	t11 = w[1];
	t14 = w[2];
	res[0] = t11 * t3 * t10 / 0.2e1 + t14 * t5 * t10 / 0.2e1 + t2 / 0.2e1;
	t18 = 0.1e1 / WAVE_SPEED;
	t24 = 0.1e1 / t7 / 0.2e1;
	t27 = t3 * t5;
	res[1] = t2 * t3 * t18 * t9 / 0.2e1 + t11 * t24 * t4 + t14 * t24 * t27;
	res[2] = t2 * t18 * t5 * t9 / 0.2e1 + t11 * t24 * t27 + t14 * t24 * t6;
}

/**
 * \fn wave_2d_get_proj_neg_null_prod_w
 * \brief Return the projection of the state w over the negative and null space
 *        {-, 0} of A_i.
 *        For more details see Maple worksheet doc/model_wave.mw.
 * \param[in] w      State vector
 * \param[in] n      Outgoing normal vector
 * \param[inout] res Projection of w
 */
void wave_2d_get_proj_neg_null_prod_w(const gdn_real *w, const gdn_real *n,
									  gdn_real *res)
{
	gdn_real t10;
	gdn_real t11;
	gdn_real t14;
	gdn_real t18;
	gdn_real t2;
	gdn_real t26;
	gdn_real t29;
	gdn_real t3;
	gdn_real t4;
	gdn_real t5;
	gdn_real t6;
	gdn_real t7;
	gdn_real t8;
	gdn_real t9;
	t2 = w[0];
	t3 = n[0];
	t4 = t3 * t3;
	t5 = n[1];
	t6 = t5 * t5;
	t7 = t4 + t6;
	t8 = sqrt(t7);
	t9 = 0.1e1 / t8;
	t10 = WAVE_SPEED * t9;
	t11 = w[1];
	t14 = w[2];
	res[0] = t11 * t3 * t10 / 0.2e1 + t14 * t5 * t10 / 0.2e1 + t2 / 0.2e1;
	t18 = 0.1e1 / WAVE_SPEED;
	t26 = 0.1e1 / t7 / 0.2e1;
	t29 = t3 * t5;
	res[1] = t2 * t3 * t18 * t9 / 0.2e1 + t11 * t26 * (t4 + 0.2e1 * t6) -
			 t14 * t26 * t29;
	res[2] = t2 * t18 * t5 * t9 / 0.2e1 - t11 * t26 * t29 +
			 t14 * t26 * (0.2e1 * t4 + t6);
}

/**
 * \fn wave_2d_get_proj_prod_prod_w
 * \brief Return the projection of the state w over the positive space
 *        {-} of A_i.
 *        For more details see Maple worksheet doc/model_wave.mw.
 * \param[in] w      State vector
 * \param[in] n      Outgoing normal vector
 * \param[inout] res Projection of w
 */
void wave_2d_get_proj_pos_prod_w(const gdn_real *w, const gdn_real *n,
								 gdn_real *res)
{
	gdn_real t10;
	gdn_real t11;
	gdn_real t14;
	gdn_real t18;
	gdn_real t2;
	gdn_real t24;
	gdn_real t27;
	gdn_real t3;
	gdn_real t4;
	gdn_real t5;
	gdn_real t6;
	gdn_real t7;
	gdn_real t8;
	gdn_real t9;
	t2 = w[0];
	t3 = n[0];
	t4 = t3 * t3;
	t5 = n[1];
	t6 = t5 * t5;
	t7 = t4 + t6;
	t8 = sqrt(t7);
	t9 = 0.1e1 / t8;
	t10 = WAVE_SPEED * t9;
	t11 = w[1];
	t14 = w[2];
	res[0] = -t11 * t3 * t10 / 0.2e1 - t14 * t5 * t10 / 0.2e1 + t2 / 0.2e1;
	t18 = 0.1e1 / WAVE_SPEED;
	t24 = 0.1e1 / t7 / 0.2e1;
	t27 = t3 * t5;
	res[1] = -t2 * t3 * t18 * t9 / 0.2e1 + t11 * t24 * t4 + t14 * t24 * t27;
	res[2] = -t2 * t18 * t5 * t9 / 0.2e1 + t11 * t24 * t27 + t14 * t24 * t6;
}

/**
 * \fn wave_2d_get_proj_pos_null_prod_w
 * \brief Return the projection of the state w over the positive and null space
 *        {+, 0} of A_i.
 *        For more details see Maple worksheet doc/model_wave.mw.
 * \param[in] w      State vector
 * \param[in] n      Outgoing normal vector
 * \param[inout] res Projection of w
 */
void wave_2d_get_proj_pos_null_prod_w(const gdn_real *w, const gdn_real *n,
									  gdn_real *res)
{
	gdn_real t10;
	gdn_real t11;
	gdn_real t14;
	gdn_real t18;
	gdn_real t2;
	gdn_real t26;
	gdn_real t29;
	gdn_real t3;
	gdn_real t4;
	gdn_real t5;
	gdn_real t6;
	gdn_real t7;
	gdn_real t8;
	gdn_real t9;
	t2 = w[0];
	t3 = n[0];
	t4 = t3 * t3;
	t5 = n[1];
	t6 = t5 * t5;
	t7 = t4 + t6;
	t8 = sqrt(t7);
	t9 = 0.1e1 / t8;
	t10 = WAVE_SPEED * t9;
	t11 = w[1];
	t14 = w[2];
	res[0] = -t11 * t3 * t10 / 0.2e1 - t14 * t5 * t10 / 0.2e1 + t2 / 0.2e1;
	t18 = 0.1e1 / WAVE_SPEED;
	t26 = 0.1e1 / t7 / 0.2e1;
	t29 = t3 * t5;
	res[1] = -t2 * t3 * t18 * t9 / 0.2e1 + t11 * t26 * (t4 + 0.2e1 * t6) -
			 t14 * t26 * t29;
	res[2] = -t2 * t18 * t5 * t9 / 0.2e1 - t11 * t26 * t29 +
			 t14 * t26 * (0.2e1 * t4 + t6);
}