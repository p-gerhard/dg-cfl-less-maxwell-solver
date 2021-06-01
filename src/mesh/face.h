/*
 * Copyright (c) 
 * 2019-2021 Philippe Helluy <helluy@math.unistra.fr>
 * 2020-2021 Pierre Gerhard <pierre.gerhard@gmail.com>
 *
 * gdon3d is free software; you can redistribute it and/or modify
 * it under the terms of the GPLv3 license. See LICENSE for details.
 */
#ifndef FACE_H
#define FACE_H

#include <gdon3d.h>

typedef struct Face3 {
	int node[3];
	int left;
	int locfaceleft;
} Face3;

void mesh_face_sort_node_face3(Face3 *face);

int mesh_face_compare_face3(const void *a, const void *b);

#endif