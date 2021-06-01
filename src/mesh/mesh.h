/*
 * Copyright (c) 
 * 2019-2021 Philippe Helluy <helluy@math.unistra.fr>
 * 2020-2021 Pierre Gerhard <pierre.gerhard@gmail.com>
 *
 * gdon3d is free software; you can redistribute it and/or modify
 * it under the terms of the GPLv3 license. See LICENSE for details.
 */
#ifndef MESH_H
#define MESH_H

#include <stdbool.h>
#include <gdon3d.h>

typedef struct gdn_mesh {
	char file_name[1024];
	int nbelems;
	int nbnodes;
	int nbqnodes_in;
	int nbqnodes_out;
	int nbfaces;
	int nboundaryfaces;
	int *boundaryface;
	int nmacrointerfaces;
	int *macrointerface;
	bool is_read;
	bool is_built;
	bool has_physical_groups;
	bool physical_group_is_alloc;
	bool phys_group_id_is_alloc;
	int nb_phys_groups;
	int *phys_group_id;
	int *physical_group;
	int *elem2node;
	int *elem2qnode;
	int *elem2elem;
	int *face2elem;
	int *fnode2enode;
	int max_node2elem;
	int *node2elem;
	gdn_real *node;
	gdn_real *qnode;
	gdn_real xmin[3];
	gdn_real xmax[3];
} gdn_mesh;

int mesh_get_nb_elem(const gdn_mesh *msh);

int mesh_get_nb_face(const gdn_mesh *msh);

int mesh_get_nb_loc_qnode(const gdn_mesh *msh);

int mesh_get_gauss_point_id(const gdn_mesh *msh, const int id_elem,
							const int id_loc_node);

gdn_real mesh_get_mass_ref_ij(const gdn_mesh *msh, const int id_loc_node_i,
							  const int id_loc_node_j);

gdn_real mesh_get_inv_mass_ref_ij(const gdn_mesh *msh, const int id_loc_node_i,
								  const int id_loc_node_j);

gdn_real mesh_get_internal_ref_ijk(const gdn_mesh *msh, const int id_loc_node_i,
								   const int id_loc_node_j, const int dir_k);

int mesh_get_ipg_elem_left_from_face(const gdn_mesh *msh, const int id_face,
									 const int id_loc_node_on_face);

int mesh_get_ipg_elem_right_from_face(const gdn_mesh *msh, const int id_face,
									  const int id_loc_node_on_face);

bool mesh_get_bool_is_boundary_face(const gdn_mesh *msh, const int id_face);

void mesh_get_elem_node(const gdn_mesh *msh, const int id_elem,
						gdn_real (*t4_node)[3]);

void mesh_get_elem_det(const gdn_mesh *msh, const int id_elem, gdn_real *det);

void mesh_get_elem_codtau(const gdn_mesh *msh, const int id_elem,
						  gdn_real codtau[3][3]);

void mesh_get_face_vn(const gdn_mesh *msh, const int id_face, gdn_real *vn);

gdn_real mesh_get_t4_hmin(const gdn_mesh *msh);

void gdn_gdn_mesh_free(gdn_mesh *m);

int gdn_gdn_mesh_is_empty(gdn_mesh *m);

void mesh_read_from_msh22_file(gdn_mesh *m, const char *filename);

void mesh_build_connectivity(gdn_mesh *m);

void gdn_gdn_mesh_print(gdn_mesh *m);
#endif
