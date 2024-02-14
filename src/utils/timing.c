/*
 * Copyright (c) 
 * 2019-2021 Philippe Helluy <helluy@math.unistra.fr>
 * 2020-2021 Pierre Gerhard <pierre.gerhard@gmail.com>
 *
 * gdon3d is free software; you can redistribute it and/or modify
 * it under the terms of the GPLv3 license. See LICENSE for details.
 */
#include <stdio.h>
#include <sys/time.h>

#include <gdon3d.h>

#include "timing.h"

void tic(struct timeval *start)
{
	gettimeofday(start, NULL);
}

gdn_real toc(struct timeval *start)
{
	struct timeval end;
	gettimeofday(&end, NULL);

	const gdn_real delta = ((end.tv_sec - start->tv_sec) * 1000000u +
							end.tv_usec - start->tv_usec) /
						   1.e6;

	return delta;
}