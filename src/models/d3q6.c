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
#include "model.h"
#include "d3q6.h"

/* D3Q6 kinetic model */
#define D3Q6_LAMBDA 3.
#define D3Q6_NB_V 6

static gdn_real d3q6_vi[3 * D3Q6_NB_V] = {
	D3Q6_LAMBDA,  0, 0, 0, D3Q6_LAMBDA,	 0, 0, 0, D3Q6_LAMBDA,
	-D3Q6_LAMBDA, 0, 0, 0, -D3Q6_LAMBDA, 0, 0, 0, -D3Q6_LAMBDA
};

static inline gdn_real dot3(const gdn_real u[3], const gdn_real v[3])
{
	return ((u[0] * v[0]) + (u[1] * v[1]) + (u[2] * v[2]));
}

/**
 * \fn d3q6_set_model
 * \brief For the kinetic model : "d3q6". Helper function to set gdn_model struct 
 *        with d3q6 variables and functions.
 * \param[inout] m struct.
 */
void d3q6_set_model(void *self)
{
	gdn_model *m = self;
	m->nb_v = D3Q6_NB_V;
	m->lambda = D3Q6_LAMBDA;
	m->vi = d3q6_vi;
	m->get_w = d3q6_get_w;
	m->get_feq = d3q6_get_feq;
	m->get_fR = d3q6_get_fR;
	m->omega = 1.999999999999;
}

/**
 * \fn d3q6_get_nb_v
 * \brief For the kinetic model : "d3q6". Get the number of velocities (nb_v).
 * \return Number of velocities (nb_v)
 */
int d3q6_get_nb_v(void)
{
	return D3Q6_NB_V;
}

/**
 * \fn d3q6_get_lambda
 * \brief For the kinetic model : "d3q6". Get the scaling factor (lambda).
 * \return Scaling factor (lambda)
 */
gdn_real d3q6_get_lambda(void)
{
	return D3Q6_LAMBDA;
}

/**
 * \fn d3q6__get_lambda
 * \brief For the kinetic model : "d3q6". Get a pointer to acces the velocities
 * (vi).
 * \return Pointer on velocity array (vi)
 */
gdn_real *d3q6_get_vi(void)
{
	return d3q6_vi;
}

/**
 * \fn d3q6_get_w
 * \brief Compute the macro state (w) from the distribution function (f).
 * \param[in] m     Model struct
 * \param[in] f     Distribution function
 * \param[inout] w  Computed macro state
 */
void d3q6_get_w(const void *self, const gdn_real *f, gdn_real *w)
{
	const gdn_model *m = self;
	/* Reconstruct macro quantities from micro */
	const int nb_v = m->nb_v;
	const int nb_w = m->nb_w;

	for (int iw = 0; iw < nb_w; iw++) {
		w[iw] = 0;
		for (int iv = 0; iv < nb_v; iv++)
			w[iw] += f[iw * nb_v + iv];
	}
}

/**
 * \fn d3q6_apply_M
 * \brief For the kinetic model : "d3q6". Apply the transformation 
 *        M(f) -> (w, q).
 * \param[in] m    Model struct
 * \param[in] f    Distribution function
 * \param[inout] w Moment of order 0 of the distribution function f
 * \param[inout] q Moment of order 1 & 2 (cross) of the distribution function f
 */
void d3q6_apply_M(const void *self, const gdn_real *f, gdn_real *w, gdn_real *q)
{
	/*Moments : M : {1, X, Y, Z, X^2 - Y^2, X^2 - Z^2} */
	const gdn_model *m = self;

	const int nb_v = m->nb_v;
	const int nb_w = m->nb_w;
	const int nb_m = m->nb_v - 1;
	const gdn_real *vi = m->vi;

	for (int iw = 0; iw < nb_w; iw++) {
		w[iw] = 0;
		q[nb_m * iw + 0] = 0;
		q[nb_m * iw + 1] = 0;
		q[nb_m * iw + 2] = 0;
		q[nb_m * iw + 3] = 0;
		q[nb_m * iw + 4] = 0;

		for (int iv = 0; iv < nb_v; iv++) {
			w[iw] += f[iw * nb_v + iv];
			q[nb_m * iw + 0] += f[iw * nb_v + iv] * vi[3 * iv + 0];
			q[nb_m * iw + 1] += f[iw * nb_v + iv] * vi[3 * iv + 1];
			q[nb_m * iw + 2] += f[iw * nb_v + iv] * vi[3 * iv + 2];

			q[nb_m * iw + 3] +=
				f[iw * nb_v + iv] * (+vi[3 * iv + 0] * vi[3 * iv + 0] -
									 vi[3 * iv + 1] * vi[3 * iv + 1]);

			q[nb_m * iw + 4] +=
				f[iw * nb_v + iv] * (+vi[3 * iv + 0] * vi[3 * iv + 0] -
									 vi[3 * iv + 2] * vi[3 * iv + 2]);
		}
	}
}

/**
 * \fn d3q6_apply_M_inv
 * \brief For the kinetic model : "d3q6". Apply the transformation 
 *        M^-1(w, q) -> (f).
 * \param[in] m    Model struct.
 * \param[in] w    Moment of order 0 of the distribution function f
 * \param[in] q    Moment of order 1 & 2 (cross) the distribution function f
 * \param[inout] f Distribution function
 */
void d3q6_apply_M_inv(const void *self, const gdn_real *w, const gdn_real *q,
					  gdn_real *f)
{
	const gdn_model *m = self;

	const int nb_v = m->nb_v;
	const int nb_w = m->nb_w;
	const int nb_m = m->nb_v - 1;

	const gdn_real l = m->lambda;
	const gdn_real l2 = m->lambda * m->lambda;

	/* DEGUG */
	// for (int iw = 0; iw < nb_w; iw++) {
	//   gdn_real q_1[3] = {q[nb_m * iw + 0],
	//                      q[nb_m * iw + 1],
	//                      q[nb_m * iw + 2]};

	//   for (int iv = 0; iv < nb_v; iv ++) {
	//     gdn_real v[3]   = {vi[3 * iv + 0],
	//                        vi[3 * iv + 1],
	//                        vi[3 * iv + 2]};

	//     f[iw * nb_v + iv] = w[iw] / 6. + dot3(q_1, v) / (2. * l2);
	//   }
	// }

	// return;

	/* Since we do not multiply per vi on the numerator we divide per l and l2 */
	for (int iw = 0; iw < nb_w; iw++) {
		f[iw * nb_v + 0] = +w[iw] / 6. + q[nb_m * iw + 0] / (2. * l) +
						   q[nb_m * iw + 3] / (6. * l2) +
						   q[nb_m * iw + 4] / (6. * l2);

		f[iw * nb_v + 1] = +w[iw] / 6. + q[nb_m * iw + 1] / (2. * l) -
						   q[nb_m * iw + 3] / (3. * l2) +
						   q[nb_m * iw + 4] / (6. * l2);

		f[iw * nb_v + 2] = +w[iw] / 6. + q[nb_m * iw + 2] / (2. * l) +
						   q[nb_m * iw + 3] / (6. * l2) -
						   q[nb_m * iw + 4] / (3. * l2);

		f[iw * nb_v + 3] = +w[iw] / 6. - q[nb_m * iw + 0] / (2. * l) +
						   q[nb_m * iw + 3] / (6. * l2) +
						   q[nb_m * iw + 4] / (6. * l2);

		f[iw * nb_v + 4] = +w[iw] / 6. - q[nb_m * iw + 1] / (2. * l) -
						   q[nb_m * iw + 3] / (3. * l2) +
						   q[nb_m * iw + 4] / (6. * l2);

		f[iw * nb_v + 5] = +w[iw] / 6. - q[nb_m * iw + 2] / (2. * l) +
						   q[nb_m * iw + 3] / (6. * l2) -
						   q[nb_m * iw + 4] / (3. * l2);
	}
}

/**
 * \fn d3q6_get_qM
 * \brief For the kinetic model : "d3q6". Compute the physical macro flux 
 *        qM = (A^i).w of a given macro state w.
 * \param[in] m     Model struct
 * \param[in] w     Macro state
 * \param[inout] qM Physical macro flux
 */
void d3q6_get_qM(const void *self, const gdn_real *w, gdn_real *qM)
{
	const gdn_model *m = self;
	const int nb_w = m->nb_w;
	const int nb_m = m->nb_v - 1;

	/* WARNING : Allocation on stack using VLA */
	gdn_real flux_1[nb_w];
	gdn_real flux_2[nb_w];
	gdn_real flux_3[nb_w];

	const gdn_real n1[3] = { 1, 0, 0 };
	const gdn_real n2[3] = { 0, 1, 0 };
	const gdn_real n3[3] = { 0, 0, 1 };

	m->get_num_flux(w, w, n1, flux_1);
	m->get_num_flux(w, w, n2, flux_2);
	m->get_num_flux(w, w, n3, flux_3);

	for (int iw = 0; iw < nb_w; iw++) {
		qM[nb_m * iw + 0] = flux_1[iw];
		qM[nb_m * iw + 1] = flux_2[iw];
		qM[nb_m * iw + 2] = flux_3[iw];
		qM[nb_m * iw + 3] = 0; /* Imposing null value at equilibrium */
		qM[nb_m * iw + 4] = 0; /* Imposing null value at equilibrium */
	}
}

/**
 * \fn d3q6_apply_Q
 * \brief For the kinetic model : "d3q6". Apply the transformation 
 *        Q(f, qM(w(f))) -> (w, y).
 * \param[in] m    Model struct
 * \param[in] f    Distribution function
 * \param[inout] w Moment of order 0 of the distribution function f
 * \param[inout] y Distance between the moment of order 1 & 2 (cross) (q) and 
 *                 the physical flux (qM) associated to w
 */
void d3q6_apply_Q(const void *self, const gdn_real *f, gdn_real *w, gdn_real *y)
{
	const gdn_model *m = self;
	const int nb_w = m->nb_w;
	const int nb_m = m->nb_v - 1;

	/* Warning allocation on stack using VLA */
	gdn_real q[nb_m * nb_w];
	gdn_real qM[nb_m * nb_w];

	d3q6_apply_M(m, f, w, q);
	d3q6_get_qM(m, w, qM);

	for (int iw = 0; iw < nb_w; iw++) {
		y[nb_m * iw + 0] = q[nb_m * iw + 0] - qM[nb_m * iw + 0];
		y[nb_m * iw + 1] = q[nb_m * iw + 1] - qM[nb_m * iw + 1];
		y[nb_m * iw + 2] = q[nb_m * iw + 2] - qM[nb_m * iw + 2];
		y[nb_m * iw + 3] = q[nb_m * iw + 3] - qM[nb_m * iw + 3];
		y[nb_m * iw + 4] = q[nb_m * iw + 4] - qM[nb_m * iw + 4];
	}
}

/**
 * \fn d3q6_apply_Q_inv
 * \brief For the kinetic model : "d3q6". Apply the transformation 
 *        Q^-1(w, y) -> (f).
 * \param[in] m    Model struct
 * \param[in] w    Macro state
 * \param[in] y    Distance between the moment of order 1 & 2 (cross) (q) and 
 *                 the physical flux qM associated to w
 * \param[inout] f Distribution function
 */
void d3q6_apply_Q_inv(const void *self, const gdn_real *w, const gdn_real *y,
					  gdn_real *f)
{
	const gdn_model *m = self;
	const int nb_w = m->nb_w;
	const int nb_m = m->nb_v - 1;

	/* Warning allocation on stack using VLA */
	gdn_real q[nb_m * nb_w];
	gdn_real qM[nb_m * nb_w];

	d3q6_get_qM(m, w, qM);

	for (int iw = 0; iw < nb_w; iw++) {
		q[nb_m * iw + 0] = y[nb_m * iw + 0] + qM[nb_m * iw + 0];
		q[nb_m * iw + 1] = y[nb_m * iw + 1] + qM[nb_m * iw + 1];
		q[nb_m * iw + 2] = y[nb_m * iw + 2] + qM[nb_m * iw + 2];
		q[nb_m * iw + 3] = y[nb_m * iw + 3] + qM[nb_m * iw + 3];
		q[nb_m * iw + 4] = y[nb_m * iw + 4] + qM[nb_m * iw + 4];
	}
	d3q6_apply_M_inv(m, w, q, f);
}

/**
 * \fn d3q6_get_feq
 * \brief For the kinetic model : "d3q6". Compute the equilibrium distribution
 *        function feq from a macro state w.
 * \param[in] w       Macro state
 * \param[inout] feq  Equilibrium distribution function
 */
void d3q6_get_feq(const void *self, const gdn_real *w, gdn_real *feq)
{
	const gdn_model *m = self;
	const int nb_w = m->nb_w;
	const int nb_m = m->nb_v - 1;
	/* Warning allocation on stack using VLA */
	gdn_real qM[nb_m * nb_w];

	d3q6_get_qM(m, w, qM);

	d3q6_apply_M_inv(m, w, qM, feq);
}

/**
 * \fn d3q6_get_fR
 * \brief   For the kinetic model : "d3q6". Compute the value of the right 
 *          distribution function fR on the boundary. 
 *          We Apply the transformation Q(fL) -> (wL, yL), we get wR from wL 
 *          and w_inc using projectors and finally we apply the transformation 
 *          Q^-1(wR, yL) -> (fR).
 * \param[in] m      Model struct
 * \param[in] vn     Outgoing normal vector
 * \param[in] w_inc  Moment of order 0
 * \param[in] qM     Physical macro flux
 * \param[inout] y   Quantity of interest
 * \param[inout] f   Distribution function
 */
void d3q6_get_fR(const void *self, const gdn_real *fL, const gdn_real *vn,
				 const gdn_real *wI, gdn_real *fR)
{
	const gdn_model *m = self;
	const int nb_w = m->nb_w;
	const int nb_m = m->nb_v - 1;

	const bool use_feq = true;
	if (use_feq) {
		/* Compute f_eq from w and the physical flux */
		d3q6_get_feq(m, wI, fR);
		return;
	}

	/* Warning allocation on stack using VLA */
	gdn_real wL[nb_w];
	gdn_real wR[nb_w];
	gdn_real proj_pos_prod_wL[nb_w];
	gdn_real proj_neg_prod_wI[nb_w];
	gdn_real yL[nb_m * nb_w];

	/* Compute Q(fL) -> (wL, yL) */
	d3q6_apply_Q(m, fL, wL, yL);

	/* Compute wR using projectors */
	m->get_proj_pos(wL, vn, proj_pos_prod_wL); /* P^{+,0}.wL */
	m->get_proj_neg(wI, vn, proj_neg_prod_wI); /* P^{-}.wI   */

	for (int iw = 0; iw < nb_w; iw++) {
		wR[iw] = proj_pos_prod_wL[iw] + proj_neg_prod_wI[iw];
	}

	/* Compute (with yR = yL)  Q^-1(wR, yR) -> fR */
	d3q6_apply_Q_inv(m, wR, yL, fR);
}