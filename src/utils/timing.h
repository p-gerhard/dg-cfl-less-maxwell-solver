/*
 * Copyright (c) 
 * 2019-2021 Philippe Helluy <helluy@math.unistra.fr>
 * 2020-2021 Pierre Gerhard <pierre.gerhard@gmail.com>
 *
 * gdon3d is free software; you can redistribute it and/or modify
 * it under the terms of the GPLv3 license. See LICENSE for details.
 */
#ifndef UTILS_TIMING_H
#define UTILS_TIMING_H

#include <sys/time.h>

void tic(struct timeval *start);

gdn_real toc(struct timeval *start);

#endif