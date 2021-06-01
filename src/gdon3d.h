/*
 * Copyright (c) 
 * 2019-2021 Philippe Helluy <helluy@math.unistra.fr>
 * 2020-2021 Pierre Gerhard <pierre.gerhard@gmail.com>
 *
 * gdon3d is free software; you can redistribute it and/or modify
 * it under the terms of the GPLv3 license. See LICENSE for details.
 */
#ifndef GDON3D_H
#define GDON3D_H

#ifdef __cplusplus
extern "C" {
#endif

#include <complex.h>
#include <float.h>
#include <limits.h>

/* Setting up double or simple precision */
#ifndef _DOUBLE_PRECISION
  #define gdn_real float
  #define gdn_complex complex
  #ifndef _COMPLEX
    #define gdn_value float
  #else
    #define gdn_value complex
  #endif
#else
  #define gdn_real double
  #define gdn_complex double _Complex
  #ifndef _COMPLEX
    #define gdn_value double
  #else
    #define gdn_value double complex
  #endif
#endif

/* Limits */
#ifdef _DOUBLE_PRECISION
  #define _VERY_SMALL (1e-12)
  #define _SMALL (1e-6)
#else
  #define  _VERY_SMALL (1e-6)
  #define _SMALL (1e-4)
#endif

/* Limits for robust_phy2ref */
#ifdef _DOUBLE_PRECISION
  #define _GDN_XREF_BOUND_THRESHOLD DBL_EPSILON
#else
  #define _GDN_XREF_BOUND_THRESHOLD FLT_EPSILON
#endif

#define _GDN_XREF_UPBOUND (1.0+_GDN_XREF_BOUND_THRESHOLD)
#define _GDN_XREF_LOWBOUND (-_GDN_XREF_BOUND_THRESHOLD)

#ifdef __cplusplus
}
#endif
#endif


