/*
 * Copyright (c) 
 * 2019-2021 Philippe Helluy <helluy@math.unistra.fr>
 * 2020-2021 Pierre Gerhard <pierre.gerhard@gmail.com>
 *
 * gdon3d is free software; you can redistribute it and/or modify
 * it under the terms of the GPLv3 license. See LICENSE for details.
 */
#ifndef D2Q4_H
#define D2Q4_H

#include <gdon3d.h>

void d2q4_set_model(void *self);

int d2q4_get_nb_v(void);

gdn_real d2q4_get_lambda(void);

gdn_real *d2q4_get_vi(void);

void d2q4_get_w(const void *self, const gdn_real *f, gdn_real *w);

void d2q4_apply_M(const void *self, const gdn_real *f, gdn_real *w, 
                  gdn_real *q);

void d2q4_apply_M_inv(const void *self, const gdn_real *w, const gdn_real *q,
                      gdn_real *f);                  

void d2q4_get_qM(const void *self, const gdn_real *w, gdn_real *qM);

void d2q4_apply_Q(const void *self, const gdn_real *f, gdn_real *w, 
                  gdn_real *y);

void d2q4_apply_Q_inv(const void *self, const gdn_real *w, const gdn_real *y, 
                      gdn_real *f);

void d2q4_get_feq(const void *self, const gdn_real *w, gdn_real *feq);

void d2q4_get_fR(const void *self, const gdn_real *fL, const gdn_real *vn, 
                 const gdn_real *wI, gdn_real *fR);            
#endif
 
