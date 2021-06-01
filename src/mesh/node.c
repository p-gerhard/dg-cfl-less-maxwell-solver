/*
 * Copyright (c) 
 * 2019-2021 Philippe Helluy <helluy@math.unistra.fr>
 * 2020-2021 Pierre Gerhard <pierre.gerhard@gmail.com>
 *
 * gdon3d is free software; you can redistribute it and/or modify
 * it under the terms of the GPLv3 license. See LICENSE for details.
 */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <gdon3d.h>
#include "node.h"

#define NODE_COMPARE_TOL 1e-8

static int compare_real(const gdn_real a, const gdn_real b)
{
	if (a < b)
		return -1;

	if (a > b)
		return 1;

	if (fabs(a - b) < NODE_COMPARE_TOL)
		return 0;

	/* WARNING : CHECK THIS RETURN */
	return 0;
}

static int compare_int(const int a, const int b)
{
	return a - b;
}

int mesh_node_compare_coordinate(const void *a, const void *b)
{
	Node *n1 = (Node *)a;
	Node *n2 = (Node *)b;

	int r = compare_real(n1->x[0], n2->x[0]);

	if (r == 0)
		r = compare_real(n1->x[1], n2->x[1]);

	if (r == 0)
		r = compare_real(n1->x[2], n2->x[2]);

	return r;
}

int mesh_node_compare_elem_id(const void *a, const void *b)
{
	Node *n1 = (Node *)a;
	Node *n2 = (Node *)b;

	int r = compare_int(n1->id_elem, n2->id_elem);

	if (r == 0)
		r = compare_int(n1->id_node, n2->id_node);

	return r;
}

void mesh_node_print_node(Node *node, int node_index)
{
	printf("%-6d"
		   "%-6d"
		   "%-8d"
		   "%-12f"
		   "%-12f"
		   "%-12f"
		   "%-7d"
		   "%-7d\n",
		   node_index, node->rank_e, node->rank_c, node->x[0], node->x[1],
		   node->x[2], node->id_elem, node->n_copy);
}

void mesh_node_print_node_list(Node *list_node, int list_size, char *list_name)
{
	printf("\n\nPrinting list      : %s\n", list_name);
	printf("Number of elements : %d\n", list_size);

	printf("%-6s"
		   "%-6s"
		   "%-8s"
		   "%-12s"
		   "%-12s"
		   "%-12s"
		   "%-7s"
		   "%-7s\n",
		   "index", "rank_e", "rank_c", "x", "y", "z", "id_elem", "n_copy");

	for (int node_index = 0; node_index < list_size; node_index++) {
		mesh_node_print_node(&list_node[node_index], node_index);
	}
}