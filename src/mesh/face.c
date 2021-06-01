/*
 * Copyright (c) 
 * 2019-2021 Philippe Helluy <helluy@math.unistra.fr>
 * 2020-2021 Pierre Gerhard <pierre.gerhard@gmail.com>
 *
 * gdon3d is free software; you can redistribute it and/or modify
 * it under the terms of the GPLv3 license. See LICENSE for details.
 */
#include <stdlib.h>

#include "face.h"
#include <gdon3d.h>

static int compare_int(const void *a, const void *b)
{
	return (*(int *)a - *(int *)b);
}

// Sort the nodes list of the face
void mesh_face_sort_node_face3(Face3 *face)
{
	qsort(face->node, 3, sizeof(int), compare_int);
}

int mesh_face_compare_face3(const void *a, const void *b)
{
	Face3 *f1 = (Face3 *)a;
	Face3 *f2 = (Face3 *)b;

	int r = f1->node[0] - f2->node[0];

	if (r == 0)
		r = f1->node[1] - f2->node[1];

	if (r == 0)
		r = f1->node[2] - f2->node[2];

	return r;
};
