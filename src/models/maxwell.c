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
#include "maxwell.h"

/**
 * \fn maxwell_imposed_macro_cos
 * \brief Imposed oscillatory solution of frequence nu for Maxwell's equations.
 *        This solution is constructed using the method of characteristics, 
 *        therby it is an exact solution for all (x, t) of Maxwell's equations.
 * \param[in] x    Position vector
 * \param[in] t    Time
 * \param[inout] w Macro state
 */
void maxwell_imposed_macro_cos(const gdn_real *x, const gdn_real t, gdn_real *w)
{
	// const gdn_real nu = 0.5;
	const gdn_real nu = 5;
	const gdn_real s = cos(nu * M_PI * (x[0] - t));
	// const gdn_real s = (x[0] - t) * (x[0] - t);
	w[0] = 0;
	w[1] = 0;
	w[2] = s;
	w[3] = 0;
	w[4] = -s;
	w[5] = 0;
}

/**
 * \fn maxwell_imposed_macro_cavity
 * \brief Imposed 3D cavity solution for Maxwell's equations.
 *        The cavity solution is constant over the time, therby it is an exact 
 *        solution for all (x, t) of Maxwell's equations.
 * \param[in] x    Position vector
 * \param[in] t    Time
 * \param[inout] w Macro state
 */
void maxwell_imposed_macro_cavity(const gdn_real *x, const gdn_real t,
								  gdn_real *w)
{
	w[0] = 0;
	w[1] = 0;
	w[2] = sin(0.3141592654e1 * x[0]) * sin(0.3141592654e1 * x[1]) *
		   cos(sqrt(0.2e1) * 0.3141592654e1 * t);

	w[3] = -sin(sqrt(0.2e1) * 0.3141592654e1 * t) * sqrt(0.2e1) *
		   sin(0.3141592654e1 * x[0]) * cos(0.3141592654e1 * x[1]) / 0.2e1;

	w[4] = sin(sqrt(0.2e1) * 0.3141592654e1 * t) * sqrt(0.2e1) *
		   cos(0.3141592654e1 * x[0]) * sin(0.3141592654e1 * x[1]) / 0.2e1;

	w[5] = 0;
}

/**
 * \fn maxwell_num_flux_upwind
 * \brief Upwind numerical flux for Maxwell's equations.
 *        Data layout: w = {Ex, Ey, Ez, Hx, Hy, Hz}. Let
 *        {{E}} = ( ER + EL ) / 2  
 *        and  
 *        [[E]] = ( ER - EL ) / 2.
 *        
 *        The first three components of the flux are :
 *        -n x {{H}} + n x n x [[E]] / r
 *        and the three last are :
 *         n x {{E}} + n x n x [[H]] / r
 * \param[in] wL      Left state
 * \param[in] wR      Right state
 * \param[in] vnorm   Outgoing normal vector
 * \param[inout] flux Numerical flux
 */
void maxwell_num_flux_upwind(const gdn_real *wL, const gdn_real *wR,
							 const gdn_real *vnorm, gdn_real *flux)
{
	const gdn_real nx = vnorm[0];
	const gdn_real ny = vnorm[1];
	const gdn_real nz = vnorm[2];

	const gdn_real overr = 1.0 / (sqrt(nx * nx + ny * ny + nz * nz) + 1e-16);
	const gdn_real nxy = overr * nx * ny;
	const gdn_real nxz = overr * nx * nz;
	const gdn_real nyz = overr * ny * nz;
	const gdn_real nxx = overr * nx * nx;
	const gdn_real nyy = overr * ny * ny;
	const gdn_real nzz = overr * nz * nz;

	const gdn_real Esx = 0.5 * (wR[0] + wL[0]);
	const gdn_real Esy = 0.5 * (wR[1] + wL[1]);
	const gdn_real Esz = 0.5 * (wR[2] + wL[2]);

	const gdn_real Hsx = 0.5 * (wR[3] + wL[3]);
	const gdn_real Hsy = 0.5 * (wR[4] + wL[4]);
	const gdn_real Hsz = 0.5 * (wR[5] + wL[5]);

	const gdn_real Edx = 0.5 * (wR[0] - wL[0]);
	const gdn_real Edy = 0.5 * (wR[1] - wL[1]);
	const gdn_real Edz = 0.5 * (wR[2] - wL[2]);

	const gdn_real Hdx = 0.5 * (wR[3] - wL[3]);
	const gdn_real Hdy = 0.5 * (wR[4] - wL[4]);
	const gdn_real Hdz = 0.5 * (wR[5] - wL[5]);

	/* E flux */
	flux[0] = nz * Hsy - ny * Hsz - (nyy + nzz) * Edx + nxy * Edy + nxz * Edz;
	flux[1] = -nz * Hsx + nx * Hsz + nxy * Edx - (nxx + nzz) * Edy + nyz * Edz;
	flux[2] = -nx * Hsy + ny * Hsx + nxz * Edx + nyz * Edy - (nxx + nyy) * Edz;

	/* H flux */
	flux[3] = -nz * Esy + ny * Esz - (nyy + nzz) * Hdx + nxy * Hdy + nxz * Hdz;
	flux[4] = nz * Esx - nx * Esz + nxy * Hdx - (nxx + nzz) * Hdy + nyz * Hdz;
	flux[5] = -ny * Esx + nx * Esy + nxz * Hdx + nyz * Hdy - (nxx + nyy) * Hdz;
}

/**
 * \fn maxwell_num_flux_boundary_metal
 * \brief Metal numerical flux (on boundaries)  for Maxwell's equations
 *        Data layout: w = {Ex, Ey, Ez, Hx, Hy, Hz}. The flux is given by :       
 *        (-n x HL, 0)^T (see Thomas Strub phd manuscrit)
 *        the right state wR and the electrical simulation E are not used.
 * \param[in] wL      Left state
 * \param[in] wR      Right state
 * \param[in] vnorm   Outgoing normal vector
 * \param[inout] flux Numerical flux
 */
void maxwell_num_flux_boundary_metal(const gdn_real *wL, const gdn_real *wR,
									 const gdn_real *vnorm, gdn_real *flux)
{
	const gdn_real nx = vnorm[0];
	const gdn_real ny = vnorm[1];
	const gdn_real nz = vnorm[2];

	const gdn_real Hx = wL[3];
	const gdn_real Hy = wL[4];
	const gdn_real Hz = wL[5];

	flux[0] = nz * Hy - ny * Hz;
	flux[1] = -nz * Hx + nx * Hz;
	flux[2] = -nx * Hy + ny * Hx;

	flux[3] = 0;
	flux[4] = 0;
	flux[5] = 0;
}

/**
 * \fn maxwell_get_proj_neg_prod_w
 * \brief Return the projection of the state w over the negative space
 *        {-} of A_i.
 *        For more details see Maple worksheet doc/model_maxwell.mw.
 * \param[in] w      State vector
 * \param[in] n      Outgoing normal vector
 * \param[inout] res Projection of w
 */
void maxwell_get_proj_neg_prod_w(const gdn_real *w, const gdn_real *n,
								 gdn_real *res)
{
	gdn_real t1;
	gdn_real t10;
	gdn_real t11;
	gdn_real t12;
	gdn_real t14;
	gdn_real t15;
	gdn_real t16;
	gdn_real t18;
	gdn_real t19;
	gdn_real t2;
	gdn_real t20;
	gdn_real t22;
	gdn_real t23;
	gdn_real t24;
	gdn_real t25;
	gdn_real t28;
	gdn_real t29;
	gdn_real t3;
	gdn_real t32;
	gdn_real t35;
	gdn_real t37;
	gdn_real t39;
	gdn_real t4;
	gdn_real t42;
	gdn_real t48;
	gdn_real t59;
	gdn_real t6;
	gdn_real t61;
	gdn_real t67;
	gdn_real t7;
	gdn_real t8;
	t1 = n[1];
	t2 = t1 * t1;
	t3 = n[2];
	t4 = t3 * t3;
	t6 = n[0];
	t7 = t6 * t6;
	t8 = t7 + t2 + t4;
	t10 = 0.1e1 / t8 / 0.2e1;
	t11 = t10 * (t2 + t4);
	t12 = w[0];
	t14 = t6 * t1;
	t15 = w[1];
	t16 = t15 * t10;
	t18 = t6 * t3;
	t19 = w[2];
	t20 = t19 * t10;
	t22 = sqrt(t8);
	t23 = 0.1e1 / t22;
	t24 = t3 * t23;
	t25 = w[4];
	t28 = t1 * t23;
	t29 = w[5];
	res[0] = t12 * t11 - t16 * t14 - t20 * t18 - t25 * t24 / 0.2e1 +
			 t28 * t29 / 0.2e1;
	t32 = t12 * t10;
	t35 = t10 * (t7 + t4);
	t37 = t1 * t3;
	t39 = w[3];
	t42 = t6 * t23;
	res[1] = -t32 * t14 + t15 * t35 - t20 * t37 + t39 * t24 / 0.2e1 -
			 t29 * t42 / 0.2e1;
	t48 = t10 * (t7 + t2);
	res[2] = -t32 * t18 - t16 * t37 + t19 * t48 - t39 * t28 / 0.2e1 +
			 t25 * t42 / 0.2e1;
	t59 = t25 * t10;
	t61 = t29 * t10;
	res[3] = t15 * t24 / 0.2e1 - t19 * t28 / 0.2e1 + t39 * t11 - t59 * t14 -
			 t61 * t18;
	t67 = t39 * t10;
	res[4] = -t12 * t24 / 0.2e1 + t19 * t42 / 0.2e1 - t67 * t14 + t25 * t35 -
			 t61 * t37;
	res[5] = t12 * t28 / 0.2e1 - t15 * t42 / 0.2e1 - t67 * t18 - t59 * t37 +
			 t29 * t48;
}

/**
 * \fn maxwell_get_proj_neg_null_prod_w
 * \brief Return the projection of the state w over the negative and null space
 *        {-, 0} of A_i.
 *        For more details see Maple worksheet doc/model_maxwell.mw.
 * \param[in] w      State vector
 * \param[in] n      Outgoing normal vector
 * \param[inout] res Projection of w
 */
void maxwell_get_proj_neg_null_prod_w(const gdn_real *w, const gdn_real *n,
									  gdn_real *res)
{
	gdn_real t1;
	gdn_real t11;
	gdn_real t12;
	gdn_real t13;
	gdn_real t15;
	gdn_real t16;
	gdn_real t17;
	gdn_real t19;
	gdn_real t2;
	gdn_real t20;
	gdn_real t21;
	gdn_real t23;
	gdn_real t24;
	gdn_real t25;
	gdn_real t26;
	gdn_real t29;
	gdn_real t30;
	gdn_real t33;
	gdn_real t37;
	gdn_real t39;
	gdn_real t4;
	gdn_real t41;
	gdn_real t44;
	gdn_real t5;
	gdn_real t51;
	gdn_real t6;
	gdn_real t62;
	gdn_real t64;
	gdn_real t7;
	gdn_real t70;
	gdn_real t9;
	t1 = n[0];
	t2 = t1 * t1;
	t4 = n[1];
	t5 = t4 * t4;
	t6 = n[2];
	t7 = t6 * t6;
	t9 = t2 + t5 + t7;
	t11 = 0.1e1 / t9 / 0.2e1;
	t12 = t11 * (0.2e1 * t2 + t5 + t7);
	t13 = w[0];
	t15 = t1 * t4;
	t16 = w[1];
	t17 = t16 * t11;
	t19 = t1 * t6;
	t20 = w[2];
	t21 = t20 * t11;
	t23 = sqrt(t9);
	t24 = 0.1e1 / t23;
	t25 = t6 * t24;
	t26 = w[4];
	t29 = t4 * t24;
	t30 = w[5];
	res[0] = t13 * t12 + t17 * t15 + t21 * t19 - t26 * t25 / 0.2e1 +
			 t29 * t30 / 0.2e1;
	t33 = t13 * t11;
	t37 = t11 * (t2 + 0.2e1 * t5 + t7);
	t39 = t6 * t4;
	t41 = w[3];
	t44 = t1 * t24;
	res[1] = t33 * t15 + t16 * t37 + t21 * t39 + t41 * t25 / 0.2e1 -
			 t30 * t44 / 0.2e1;
	t51 = t11 * (t2 + t5 + 0.2e1 * t7);
	res[2] = t33 * t19 + t17 * t39 + t20 * t51 - t41 * t29 / 0.2e1 +
			 t26 * t44 / 0.2e1;
	t62 = t26 * t11;
	t64 = t30 * t11;
	res[3] = t16 * t25 / 0.2e1 - t20 * t29 / 0.2e1 + t41 * t12 + t62 * t15 +
			 t64 * t19;
	t70 = t41 * t11;
	res[4] = -t13 * t25 / 0.2e1 + t20 * t44 / 0.2e1 + t70 * t15 + t26 * t37 +
			 t64 * t39;
	res[5] = t13 * t29 / 0.2e1 - t16 * t44 / 0.2e1 + t70 * t19 + t62 * t39 +
			 t30 * t51;
}

/**
 * \fn maxwell_get_proj_prod_prod_w
 * \brief Return the projection of the state w over the positive space
 *        {-} of A_i.
 *        For more details see Maple worksheet doc/model_maxwell.mw.
 * \param[in] w      State vector
 * \param[in] n      Outgoing normal vector
 * \param[inout] res Projection of w
 */
void maxwell_get_proj_pos_prod_w(const gdn_real *w, const gdn_real *n,
								 gdn_real *res)
{
	gdn_real t1;
	gdn_real t10;
	gdn_real t11;
	gdn_real t12;
	gdn_real t14;
	gdn_real t15;
	gdn_real t16;
	gdn_real t18;
	gdn_real t19;
	gdn_real t2;
	gdn_real t20;
	gdn_real t22;
	gdn_real t23;
	gdn_real t24;
	gdn_real t25;
	gdn_real t28;
	gdn_real t29;
	gdn_real t3;
	gdn_real t32;
	gdn_real t35;
	gdn_real t37;
	gdn_real t39;
	gdn_real t4;
	gdn_real t42;
	gdn_real t48;
	gdn_real t59;
	gdn_real t6;
	gdn_real t61;
	gdn_real t67;
	gdn_real t7;
	gdn_real t8;
	t1 = n[1];
	t2 = t1 * t1;
	t3 = n[2];
	t4 = t3 * t3;
	t6 = n[0];
	t7 = t6 * t6;
	t8 = t7 + t2 + t4;
	t10 = 0.1e1 / t8 / 0.2e1;
	t11 = t10 * (t2 + t4);
	t12 = w[0];
	t14 = t6 * t1;
	t15 = w[1];
	t16 = t15 * t10;
	t18 = t6 * t3;
	t19 = w[2];
	t20 = t19 * t10;
	t22 = sqrt(t8);
	t23 = 0.1e1 / t22;
	t24 = t3 * t23;
	t25 = w[4];
	t28 = t1 * t23;
	t29 = w[5];
	res[0] = t12 * t11 - t16 * t14 - t20 * t18 + t25 * t24 / 0.2e1 -
			 t28 * t29 / 0.2e1;
	t32 = t12 * t10;
	t35 = t10 * (t7 + t4);
	t37 = t1 * t3;
	t39 = w[3];
	t42 = t6 * t23;
	res[1] = -t32 * t14 + t15 * t35 - t20 * t37 - t39 * t24 / 0.2e1 +
			 t29 * t42 / 0.2e1;
	t48 = t10 * (t7 + t2);
	res[2] = -t32 * t18 - t16 * t37 + t19 * t48 + t39 * t28 / 0.2e1 -
			 t25 * t42 / 0.2e1;
	t59 = t25 * t10;
	t61 = t29 * t10;
	res[3] = -t15 * t24 / 0.2e1 + t19 * t28 / 0.2e1 + t39 * t11 - t59 * t14 -
			 t61 * t18;
	t67 = t39 * t10;
	res[4] = t12 * t24 / 0.2e1 - t19 * t42 / 0.2e1 - t67 * t14 + t25 * t35 -
			 t61 * t37;
	res[5] = -t12 * t28 / 0.2e1 + t15 * t42 / 0.2e1 - t67 * t18 - t59 * t37 +
			 t29 * t48;
}

/**
 * \fn maxwell_get_proj_pos_null_prod_w
 * \brief Return the projection of the state w over the positive and null space
 *        {+, 0} of A_i.
 *        For more details see Maple worksheet doc/model_maxwell.mw.
 * \param[in] w      State vector
 * \param[in] n      Outgoing normal vector
 * \param[inout] res Projection of w
 */
void maxwell_get_proj_pos_null_prod_w(const gdn_real *w, const gdn_real *n,
									  gdn_real *res)
{
	gdn_real t1;
	gdn_real t11;
	gdn_real t12;
	gdn_real t13;
	gdn_real t15;
	gdn_real t16;
	gdn_real t17;
	gdn_real t19;
	gdn_real t2;
	gdn_real t20;
	gdn_real t21;
	gdn_real t23;
	gdn_real t24;
	gdn_real t25;
	gdn_real t26;
	gdn_real t29;
	gdn_real t30;
	gdn_real t33;
	gdn_real t37;
	gdn_real t39;
	gdn_real t4;
	gdn_real t41;
	gdn_real t44;
	gdn_real t5;
	gdn_real t51;
	gdn_real t6;
	gdn_real t62;
	gdn_real t64;
	gdn_real t7;
	gdn_real t70;
	gdn_real t9;
	t1 = n[0];
	t2 = t1 * t1;
	t4 = n[1];
	t5 = t4 * t4;
	t6 = n[2];
	t7 = t6 * t6;
	t9 = t2 + t5 + t7;
	t11 = 0.1e1 / t9 / 0.2e1;
	t12 = t11 * (0.2e1 * t2 + t5 + t7);
	t13 = w[0];
	t15 = t1 * t4;
	t16 = w[1];
	t17 = t16 * t11;
	t19 = t1 * t6;
	t20 = w[2];
	t21 = t20 * t11;
	t23 = sqrt(t9);
	t24 = 0.1e1 / t23;
	t25 = t6 * t24;
	t26 = w[4];
	t29 = t4 * t24;
	t30 = w[5];
	res[0] = t13 * t12 + t17 * t15 + t21 * t19 + t26 * t25 / 0.2e1 -
			 t29 * t30 / 0.2e1;
	t33 = t13 * t11;
	t37 = t11 * (t2 + 0.2e1 * t5 + t7);
	t39 = t6 * t4;
	t41 = w[3];
	t44 = t1 * t24;
	res[1] = t33 * t15 + t16 * t37 + t21 * t39 - t41 * t25 / 0.2e1 +
			 t30 * t44 / 0.2e1;
	t51 = t11 * (t2 + t5 + 0.2e1 * t7);
	res[2] = t33 * t19 + t17 * t39 + t20 * t51 + t41 * t29 / 0.2e1 -
			 t26 * t44 / 0.2e1;
	t62 = t26 * t11;
	t64 = t30 * t11;
	res[3] = -t16 * t25 / 0.2e1 + t20 * t29 / 0.2e1 + t41 * t12 + t62 * t15 +
			 t64 * t19;
	t70 = t41 * t11;
	res[4] = t13 * t25 / 0.2e1 - t20 * t44 / 0.2e1 + t70 * t15 + t26 * t37 +
			 t64 * t39;
	res[5] = -t13 * t29 / 0.2e1 + t16 * t44 / 0.2e1 + t70 * t19 + t62 * t39 +
			 t30 * t51;
}