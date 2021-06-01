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
#include "d2q4.h"

/* D2Q4 kinetic model */
#define D2Q4_LAMBDA 4.
#define D2Q4_NB_V 4

/* Real D2Q4 */
// static gdn_real d2q4_vi[2 * D2Q4_NB_V] =
//   {  D2Q4_LAMBDA,            0,
//                0,  D2Q4_LAMBDA,
//     -D2Q4_LAMBDA,            0,
//                0, -D2Q4_LAMBDA};

/* D2Q4 using 3D representation (z coordinate = 0) */
static gdn_real d2q4_vi[3 * D2Q4_NB_V] = {
	D2Q4_LAMBDA, 0, 0, 0, D2Q4_LAMBDA, 0, -D2Q4_LAMBDA, 0, 0, 0, -D2Q4_LAMBDA, 0
};

/**
 * \fn d2q4_set_model
 * \brief For the kinetic model : "d2q4". Helper function to set gdn_model struct 
 *        with d2q4 variables and functions.
 * \param[inout] m struct.
 */
void d2q4_set_model(void *self)
{
	gdn_model *m = self;
	m->nb_v = D2Q4_NB_V;
	m->lambda = D2Q4_LAMBDA;
	m->vi = d2q4_vi;
	m->get_w = d2q4_get_w;
	m->get_feq = d2q4_get_feq;
	m->get_fR = d2q4_get_fR;
	m->omega = 1.999999999999;
}

/**
 * \fn d2q4_get_nb_v
 * \brief For the kinetic model : "d2q4". Get the number of velocities (nb_v).
 * \return Number of velocities (nb_v)
 */
int d2q4_get_nb_v(void)
{
	return D2Q4_NB_V;
}

/**
 * \fn d2q4_get_lambda
 * \brief For the kinetic model : "d2q4". Get the scaling factor (lambda).
 * \return Scaling factor (lambda)
 */
gdn_real d2q4_get_lambda(void)
{
	return D2Q4_LAMBDA;
}

/**
 * \fn d2q4__get_lambda
 * \brief For the kinetic model : "d2q4". Get a pointer to acces the velocities
 * (vi).
 * \return Pointer on velocity array (vi)
 */
gdn_real *d2q4_get_vi(void)
{
	return d2q4_vi;
}

/**
 * \fn d2q4_get_w
 * \brief Compute the macro state (w) from the distribution function (f).
 * \param[in] m     Model struct
 * \param[in] f     Distribution function
 * \param[inout] w  Computed macro state
 */
void d2q4_get_w(const void *self, const gdn_real *f, gdn_real *w)
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
 * \fn d2q4_apply_M
 * \brief For the kinetic model : "d2q4". Apply the transformation 
 *        M(f) -> (w, q).
 * \param[in] m    Model struct
 * \param[in] f    Distribution function
 * \param[inout] w Moment of order 0 of the distribution function f
 * \param[inout] q Moment of order 1 & 2 (cross) of the distribution function f
 */
void d2q4_apply_M(const void *self, const gdn_real *f, gdn_real *w, gdn_real *q)
{
	/*Moments : M : {1, X, Y, X^2 - Y^2} */
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

		for (int iv = 0; iv < nb_v; iv++) {
			w[iw] += f[iw * nb_v + iv];
			q[nb_m * iw + 0] += f[iw * nb_v + iv] * vi[3 * iv + 0];
			q[nb_m * iw + 1] += f[iw * nb_v + iv] * vi[3 * iv + 1];

			q[nb_m * iw + 2] +=
				f[iw * nb_v + iv] * (+vi[3 * iv + 0] * vi[3 * iv + 0] -
									 vi[3 * iv + 1] * vi[3 * iv + 1]);
		}
	}
}

/**
 * \fn d2q4_apply_M_inv
 * \brief For the kinetic model : "d2q4". Apply the transformation 
 *        M^-1(w, q) -> (f).
 * \param[in] m    Model struct.
 * \param[in] w    Moment of order 0 of the distribution function f
 * \param[in] q    Moment of order 1 & 2 (cross) the distribution function f
 * \param[inout] f Distribution function
 */
void d2q4_apply_M_inv(const void *self, const gdn_real *w, const gdn_real *q,
					  gdn_real *f)
{
	const gdn_model *m = self;

	const int nb_v = m->nb_v;
	const int nb_w = m->nb_w;
	const int nb_m = m->nb_v - 1;

	const gdn_real l = m->lambda;
	const gdn_real l2 = m->lambda * m->lambda;

	/* Since we do not multiply per vi on the numerator we divide per l and l2 */
	for (int iw = 0; iw < nb_w; iw++) {
		f[iw * nb_v + 0] = +w[iw] / 4. + q[nb_m * iw + 0] / (2. * l) +
						   q[nb_m * iw + 2] / (4. * l2);

		f[iw * nb_v + 1] = +w[iw] / 4. + q[nb_m * iw + 1] / (2. * l) -
						   q[nb_m * iw + 2] / (4. * l2);

		f[iw * nb_v + 2] = +w[iw] / 4. - q[nb_m * iw + 0] / (2. * l) +
						   q[nb_m * iw + 2] / (4. * l2);

		f[iw * nb_v + 3] = +w[iw] / 4. - q[nb_m * iw + 1] / (2. * l) -
						   q[nb_m * iw + 2] / (4. * l2);
	}
}

/**
 * \fn d2q4_get_qM
 * \brief For the kinetic model : "d2q4". Compute the physical macro flux 
 *        qM = (A^i).w of a given macro state w.
 * \param[in] m     Model struct
 * \param[in] w     Macro state
 * \param[inout] qM Physical macro flux
 */
void d2q4_get_qM(const void *self, const gdn_real *w, gdn_real *qM)
{
	const gdn_model *m = self;
	const int nb_w = m->nb_w;
	const int nb_m = m->nb_v - 1;

	/* WARNING : Allocation on stack using VLA */
	gdn_real flux_1[nb_w];
	gdn_real flux_2[nb_w];

	const gdn_real n1[3] = { 1, 0, 0 };
	const gdn_real n2[3] = { 0, 1, 0 };

	m->get_num_flux(w, w, n1, flux_1);
	m->get_num_flux(w, w, n2, flux_2);

	for (int iw = 0; iw < nb_w; iw++) {
		qM[nb_m * iw + 0] = flux_1[iw];
		qM[nb_m * iw + 1] = flux_2[iw];
		qM[nb_m * iw + 2] = 0; /* Imposing null value at equilibrium */
	}
}

/**
 * \fn d2q4_apply_Q
 * \brief For the kinetic model : "d2q4". Apply the transformation 
 *        Q(f, qM(w(f))) -> (w, y).
 * \param[in] m    Model struct
 * \param[in] f    Distribution function
 * \param[inout] w Moment of order 0 of the distribution function f
 * \param[inout] y Distance between the moment of order 1 & 2 (cross) (q) and 
 *                 the physical flux (qM) associated to w
 */
void d2q4_apply_Q(const void *self, const gdn_real *f, gdn_real *w, gdn_real *y)
{
	const gdn_model *m = self;
	const int nb_w = m->nb_w;
	const int nb_m = m->nb_v - 1;

	/* Warning allocation on stack using VLA */
	gdn_real q[nb_m * nb_w];
	gdn_real qM[nb_m * nb_w];

	d2q4_apply_M(m, f, w, q);
	d2q4_get_qM(m, w, qM);

	for (int iw = 0; iw < nb_w; iw++) {
		y[nb_m * iw + 0] = q[nb_m * iw + 0] - qM[nb_m * iw + 0];
		y[nb_m * iw + 1] = q[nb_m * iw + 1] - qM[nb_m * iw + 1];
		y[nb_m * iw + 2] = q[nb_m * iw + 2] - qM[nb_m * iw + 2];
	}
}

/**
 * \fn d2q4_apply_Q_inv
 * \brief For the kinetic model : "d2q4". Apply the transformation 
 *        Q^-1(w, y) -> (f).
 * \param[in] m    Model struct
 * \param[in] w    Macro state
 * \param[in] y    Distance between the moment of order 1 & 2 (cross) (q) and 
 *                 the physical flux qM associated to w
 * \param[inout] f Distribution function
 */
void d2q4_apply_Q_inv(const void *self, const gdn_real *w, const gdn_real *y,
					  gdn_real *f)
{
	const gdn_model *m = self;
	const int nb_w = m->nb_w;
	const int nb_m = m->nb_v - 1;

	/* Warning allocation on stack using VLA */
	gdn_real q[nb_m * nb_w];
	gdn_real qM[nb_m * nb_w];

	d2q4_get_qM(m, w, qM);

	for (int iw = 0; iw < nb_w; iw++) {
		q[nb_m * iw + 0] = y[nb_m * iw + 0] + qM[nb_m * iw + 0];
		q[nb_m * iw + 1] = y[nb_m * iw + 1] + qM[nb_m * iw + 1];
		q[nb_m * iw + 2] = y[nb_m * iw + 2] + qM[nb_m * iw + 2];
	}
	d2q4_apply_M_inv(m, w, q, f);
}

/**
 * \fn d2q4_get_feq
 * \brief For the kinetic model : "d2q4". Compute the equilibrium distribution
 *        function feq from a macro state w.
 * \param[in] w       Macro state
 * \param[inout] feq  Equilibrium distribution function
 */
void d2q4_get_feq(const void *self, const gdn_real *w, gdn_real *feq)
{
	const gdn_model *m = self;
	const int nb_w = m->nb_w;
	const int nb_m = m->nb_v - 1;
	/* Warning allocation on stack using VLA */
	gdn_real qM[nb_m * nb_w];

	d2q4_get_qM(m, w, qM);

	d2q4_apply_M_inv(m, w, qM, feq);
}

/**
 * \fn d2q4_get_fR
 * \brief   For the kinetic model : "d2q4". Compute the value of the right 
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
void d2q4_get_fR(const void *self, const gdn_real *fL, const gdn_real *vn,
				 const gdn_real *wI, gdn_real *fR)
{
	const gdn_model *m = self;
	const int nb_w = m->nb_w;
	const int nb_m = m->nb_v - 1;

	const bool use_feq = true;
	if (use_feq) {
		/* Compute f_eq from w and the physical flux */
		d2q4_get_feq(m, wI, fR);
		return;
	}

	/* Warning allocation on stack using VLA */
	gdn_real wL[nb_w];
	gdn_real wR[nb_w];
	gdn_real proj_pos_prod_wL[nb_w];
	gdn_real proj_neg_prod_wI[nb_w];
	gdn_real yL[nb_m * nb_w];

	/* Compute Q(fL) -> (wL, yL) */
	d2q4_apply_Q(m, fL, wL, yL);

	/* Compute wR using projectors */
	m->get_proj_pos(wL, vn, proj_pos_prod_wL); /* P^{+,0}.wL */
	m->get_proj_neg(wI, vn, proj_neg_prod_wI); /* P^{-}.wI   */

	for (int iw = 0; iw < nb_w; iw++) {
		wR[iw] = proj_pos_prod_wL[iw] + proj_neg_prod_wI[iw];
	}

	/* Compute (with yR = yL)  Q^-1(wR, yR) -> fR */
	d2q4_apply_Q_inv(m, wR, yL, fR);
}