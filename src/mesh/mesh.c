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
#include <errno.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "face.h"
#include <gdon3d.h>
#include "mesh.h"
#include "node.h"
#include "t4.h"
#include "t10.h"
#include "tetrageometry.h"

gdn_mesh null_gdn_mesh = { 0 };

gdn_real mesh_get_mass_ref_entry(const gdn_mesh *msh, const int i, const int j)
{
	return t10_get_mass_ref_entry(i, j);
}

gdn_real mesh_get_inv_mass_ref_entry(const gdn_mesh *msh, const int i,
									 const int j)
{
	return t10_get_inv_mass_ref_entry(i, j);
}

gdn_real mesh_get_internal_ref_entry(const gdn_mesh *msh, const int i,
									 const int j, const int id_dir)
{
	return t10_get_internal_ref_entry(i, j, id_dir);
}

int mesh_get_nb_elem(const gdn_mesh *msh)
{
	return msh->nbelems;
}

int mesh_get_nb_face(const gdn_mesh *msh)
{
	return msh->nbfaces;
}

int mesh_get_nb_loc_qnode(const gdn_mesh *msh)
{
	return NQN;
}

int mesh_get_gauss_point_id(const gdn_mesh *msh, const int id_elem,
							const int id_loc_node)
{
	int nb_loc_qnode = mesh_get_nb_loc_qnode(msh);
	int ipg = id_elem * nb_loc_qnode + id_loc_node;
	return ipg;
}

void mesh_get_elem_node(const gdn_mesh *msh, const int ie,
						gdn_real (*t4_node)[3])
{
	for (int iloc = 0; iloc < NNO; iloc++) {
		int ino = msh->elem2qnode[ie * NQN + iloc];
		t4_node[iloc][0] = msh->node[ino * 3 + 0];
		t4_node[iloc][1] = msh->node[ino * 3 + 1];
		t4_node[iloc][2] = msh->node[ino * 3 + 2];
	}
}

void mesh_get_elem_det(const gdn_mesh *msh, const int id_elem, gdn_real *det)
{
	gdn_real dtau[3][3];
	gdn_real codtau[3][3];
	gdn_real xref[3];
	gdn_real t4_node[NNO][3];

	mesh_get_elem_node(msh, id_elem, t4_node);
	gdn_tetra_ref2phy(t4_node, xref, NULL, -1, NULL, dtau, codtau, NULL, NULL);
	*det = geometry_dot_product(dtau[0], codtau[0]);
}

void mesh_get_elem_codtau(const gdn_mesh *msh, const int id_elem,
						  gdn_real codtau[3][3])
{
	gdn_real dtau[3][3];
	gdn_real xref[3];
	gdn_real t4_node[NNO][3];

	mesh_get_elem_node(msh, id_elem, t4_node);
	gdn_tetra_ref2phy(t4_node, xref, NULL, -1, NULL, dtau, codtau, NULL, NULL);
}

void mesh_get_face_vn(const gdn_mesh *msh, const int id_face, gdn_real *vn)
{
	gdn_real dtau[3][3];
	gdn_real codtau[3][3];
	gdn_real xref[3];
	gdn_real t4_node[NNO][3];

	int id_element_left = msh->face2elem[E2F * id_face + 0];
	int id_loc_face_left = msh->face2elem[E2F * id_face + 1];

	mesh_get_elem_node(msh, id_element_left, t4_node);
	gdn_tetra_ref2phy(t4_node, xref, NULL, id_loc_face_left, NULL, dtau, codtau,
					  NULL, vn);

	gdn_real cr = 1;

	if (id_face == 3) {
		cr = sqrt(3);
	}

	for (int ii = 0; ii < 3; ii++) {
		vn[ii] *= cr;
	}
}

int mesh_get_ipg_elem_right_from_face(const gdn_mesh *msh, const int id_face,
									  const int id_loc_node_on_face)
{
	int id_elem_right = msh->face2elem[E2F * id_face + 2];

	int iloc = msh->fnode2enode[id_face * 2 * NQV + NQV + id_loc_node_on_face];
	int ipg = id_elem_right * NQN + iloc;
	return ipg;
}

int mesh_get_ipg_elem_left_from_face(const gdn_mesh *msh, const int id_face,
									 const int id_loc_node_on_face)
{
	int id_elem_left = msh->face2elem[E2F * id_face + 0];

	int iloc = msh->fnode2enode[id_face * 2 * NQV + id_loc_node_on_face];
	int ipg = id_elem_left * NQN + iloc;
	return ipg;
}

gdn_real mesh_get_t4_hmin(const gdn_mesh *msh)
{
	gdn_real hmin = FLT_MAX;

	gdn_real t4_node[NNO][3];

	int nb_elem = mesh_get_nb_elem(msh);

	gdn_real xref[3], dtau[3][3], codtau[3][3], vnds[3];
	gdn_real surf_1 = 0;
	gdn_real surf_2 = 0;

	gdn_real vol_1 = 0;
	gdn_real vol_2 = 0;

	for (int ie = 0; ie < nb_elem; ie++) {
		mesh_get_elem_node(msh, ie, t4_node);

		surf_1 = 0;
		vol_1 = 0;
		surf_2 = 0;
		vol_2 = 0;

		/* Compute volume of the tetrahedron */
		gdn_tetra_ref2phy(t4_node, xref, NULL, -1, NULL, dtau, codtau, NULL,
						  NULL);
		gdn_real det = geometry_dot_product(dtau[0], codtau[0]);

		for (int ii = 0; ii < NNO; ii++) {
			vol_1 += t4_get_wpg_ref_entry(ii) * det;
		}

		/* Compute Area of the tetrahedron */
		/* Surface scaling */
		gdn_real test[4] = { 0.5, 0.5, 0.5, 0.86602540378443865 };
		for (int ifa = 0; ifa < 4; ifa++) {
			gdn_tetra_ref2phy(t4_node, xref, NULL, ifa, NULL, dtau, codtau,
							  NULL, vnds);
			gdn_real norm =
				sqrt(vnds[0] * vnds[0] + vnds[1] * vnds[1] + vnds[2] * vnds[2]);
			surf_1 += norm * test[ifa];
		}

		/* Other methode */
		geometry_get_tetra_area_and_vol(t4_node, &surf_2, &vol_2);
		assert(fabs(surf_1 - surf_2) < 1e-12);
		assert(fabs(vol_1 - vol_2) < 1e-12);

		gdn_real loc_h = vol_1 / surf_1;

		/* Check if structured */
		// if(ie > 0) {
		//   assert(fabs(hmin - loc_h) < 1e-12);
		// }
		/* Compute local hmin */
		if (loc_h < hmin) {
			hmin = loc_h;
		}
		// printf("h_min=%1.12f, s=%1.12f, v=%1.12f\n", loc_h, surf_1, vol_1);
	}
	return hmin;
}

bool mesh_get_bool_is_boundary_face(const gdn_mesh *msh, const int id_face)
{
	bool is_boundary_face = false;

	int id_elem_right = msh->face2elem[E2F * id_face + 2];

	if (id_elem_right < 0) {
		is_boundary_face = true;
	}

	return is_boundary_face;
}

void gdn_node_get_t10_elem_order(gdn_mesh *m, Node *t10_node_lst)
{
	for (int ie = 0; ie < m->nbelems; ie++) {
		/* Get nodes id of the current T4 */
		int id[NNO] = { m->elem2node[ie * NNO + 0], m->elem2node[ie * NNO + 1],
						m->elem2node[ie * NNO + 2],
						m->elem2node[ie * NNO + 3] };

		/* Get the coordinates xyz of the four T4 nodes */
		gdn_real t4_node[NNO][3] = {
			{ m->node[id[0] * 3], m->node[id[0] * 3 + 1],
			  m->node[id[0] * 3 + 2] },
			{ m->node[id[1] * 3], m->node[id[1] * 3 + 1],
			  m->node[id[1] * 3 + 2] },
			{ m->node[id[2] * 3], m->node[id[2] * 3 + 1],
			  m->node[id[2] * 3 + 2] },
			{ m->node[id[3] * 3], m->node[id[3] * 3 + 1],
			  m->node[id[3] * 3 + 2] }
		};

		/* Fill t10_node_lst using element order  */
		for (int iqn = 0; iqn < NQN; iqn++) {
			int idx = NQN * ie + iqn;

			t10_node_lst[idx].id_elem = ie;
			t10_node_lst[idx].id_node = iqn;
			t10_node_lst[idx].rank_e = idx;
			t10_node_lst[idx].rank_c = 0;
			t10_node_lst[idx].n_copy = 1;

			/* Node coordinates */
			gdn_real t10_ref_node_coord[3];
			t10_get_ref_node_coord(iqn, t10_ref_node_coord);
			gdn_tetra_ref2phy(t4_node, t10_ref_node_coord, NULL, 0,
							  t10_node_lst[idx].x, NULL, NULL, NULL, NULL);
		}
	}
}

void gdn_gdn_mesh_free(gdn_mesh *m)
{
	m->nbelems = 0;
	m->nbnodes = 0;
	m->nbqnodes_in = 0;
	m->nbfaces = 0;
	m->nboundaryfaces = 0;

	if (m->boundaryface) {
		free(m->boundaryface);
		m->boundaryface = NULL;
	}

	m->nmacrointerfaces = 0;

	if (m->macrointerface) {
		free(m->macrointerface);
		m->macrointerface = NULL;
	}

	m->is_read = false;
	m->is_built = false;
	m->has_physical_groups = false;
	m->physical_group_is_alloc = false;
	m->phys_group_id_is_alloc = false;
	m->nb_phys_groups = 0;

	if (m->phys_group_id) {
		free(m->phys_group_id);
		m->phys_group_id = NULL;
	}

	if (m->physical_group) {
		free(m->physical_group);
		m->physical_group = NULL;
	}

	if (m->elem2node) {
		free(m->elem2node);
		m->elem2node = NULL;
	}

	if (m->elem2qnode) {
		free(m->elem2qnode);
		m->elem2qnode = NULL;
	}

	if (m->elem2elem) {
		free(m->elem2elem);
		m->elem2elem = NULL;
	}

	if (m->face2elem) {
		free(m->face2elem);
		m->face2elem = NULL;
	}

	if (m->fnode2enode) {
		free(m->fnode2enode);
		m->fnode2enode = NULL;
	}

	m->max_node2elem = 0;

	if (m->node2elem) {
		free(m->node2elem);
		m->node2elem = NULL;
	}

	if (m->node) {
		free(m->node);
		m->node = NULL;
	}
}

int fp_move_to_pattern_line(const char *pattern, FILE *f)
{
	char *line_buff = NULL;
	size_t line_buff_size = 0;
	ssize_t ret = 1;

	/* Read lines in f until pattern is found or EOF is reached */
	do {
		if ((ret = getline(&line_buff, &line_buff_size, f)) == -1) {
			break;
		}
	} while (strcmp(line_buff, pattern) != 0);

	int ret_code = 0;
	/* Reached EOF */
	if (ret == -1 && errno == 0) {
		printf(
			"warning: fp_move_to_pattern_line: didn't find the pattern : %s \n",
			pattern);
		ret_code = 0;
		/* Other error */
	} else if (ret == -1 && errno != 0) {
		fprintf(stderr, "getline: %s\n", strerror(errno));
		ret_code = -1;
	} else {
		assert(strcmp(line_buff, pattern) == 0);
		ret_code = 1;
	}

	free(line_buff);
	return ret_code;
}

int gmsh_v2_extract_meshformat(float *version, bool *is_ascii, int *data_size,
							   FILE *f)
{
	char *line_buff = NULL;
	size_t line_buff_size = 0;
	ssize_t ret = 1;
	int delim = (int)' ';

	/* Move file pointer to the $MeshFormat line */
	assert(fp_move_to_pattern_line("$MeshFormat\n", f));

	/* Field : version-number (float) */
	if ((ret = getdelim(&line_buff, &line_buff_size, delim, f)) != -1) {
		*version = atof(line_buff);
	}

	/* Field : file-type (int) */
	if ((ret = getdelim(&line_buff, &line_buff_size, delim, f)) != -1) {
		if (atoi(line_buff) == 0) {
			*is_ascii = true;
		} else {
			*is_ascii = false;
		}
	}

	/* Field : data-size (int) */
	if ((ret = getline(&line_buff, &line_buff_size, f)) != -1) {
		*data_size = atoi(line_buff);
	}

	/*Check end of $MeshFormat block */
	if ((ret = getline(&line_buff, &line_buff_size, f)) != -1) {
		assert(strcmp(line_buff, "$EndMeshFormat\n") == 0);
	}

	free(line_buff);
	return 1;
}

int gmsh_v2_extract_physicalnames(gdn_mesh *m, FILE *f)
{
	char *line_buff = NULL;
	size_t line_buff_size = 0;
	ssize_t ret = 1;
	int delim = (int)' ';

	m->has_physical_groups = true;
	m->physical_group_is_alloc = false;

	int number_of_names = 0;
	int physical_dimension = 0;
	int physical_number = 0;
	char *physical_name = NULL;

	/* Field : number-of-names (int) */
	if ((ret = getline(&line_buff, &line_buff_size, f)) != -1) {
		number_of_names = atoi(line_buff);
	}

	for (int i = 0; i < number_of_names; i++) {
		/* Field : physical-dimension (int) */
		if ((ret = getdelim(&line_buff, &line_buff_size, delim, f)) != -1) {
			physical_dimension = atoi(line_buff);
		}

		/* Field : physical-number (int) */
		if ((ret = getdelim(&line_buff, &line_buff_size, delim, f)) != -1) {
			physical_number = atoi(line_buff);
		}

		/* Field : physical-name (str) */
		if ((ret = getline(&line_buff, &line_buff_size, f)) != -1) {
			/* We remove the first char '"' and last '"' + '\n' at the end of
					* line_buff */
			int str_len = strlen(line_buff) - 3;
			physical_name = (char *)calloc(str_len, sizeof(char));
			strncpy(physical_name, line_buff + 1, str_len);
			free(physical_name);
		}
	}
	free(line_buff);
	return 1;
}

int gmsh_v2_extract_nodes(gdn_mesh *m, FILE *f)
{
	char *line_buff = NULL;
	size_t line_buff_size = 0;
	ssize_t ret = 1;
	int delim = (int)' ';

	/* Field : number-of-nodes (int) */
	if ((ret = getline(&line_buff, &line_buff_size, f)) != -1) {
		m->nbnodes = atoi(line_buff);
	}

	/*TODO : Deal with different physical dim = 1, 2 ,3 of nodes ? */

	/* Allocate node array */
	if (m->node) {
		free(m->node);
	}

	m->node = (gdn_real *)calloc(3 * (m->nbnodes), sizeof(gdn_real));
	assert(m->node);

	int node_number = 0;
	for (int i = 0; i < m->nbnodes; i++) {
		/* Field : node-number (int) */
		if ((ret = getdelim(&line_buff, &line_buff_size, delim, f)) != -1) {
			node_number = atoi(line_buff);
		}
		assert(node_number - i == 1);

		/* Field : x-coord (gdn_real) */
		if ((ret = getdelim(&line_buff, &line_buff_size, delim, f)) != -1) {
			m->node[3 * i + 0] = atof(line_buff);
		}

		/* Field : y-coord (gdn_real) */
		if ((ret = getdelim(&line_buff, &line_buff_size, delim, f)) != -1) {
			m->node[3 * i + 1] = atof(line_buff);
		}

		/* Field : z-coord (gdn_real) (read from current fp position to last \n) */
		if ((ret = getline(&line_buff, &line_buff_size, f)) != -1) {
			m->node[3 * i + 2] = atof(line_buff);
		}
	}

	free(line_buff);
	return 1;
}

int gmsh_v2_extract_elements(gdn_mesh *m, FILE *f)
{
	char *line_buff = NULL;
	size_t line_buff_size = 0;
	ssize_t ret = 1;
	int delim = (int)' ';

	int nb_all_elem = 0;

	/* Field : number-of-elements (int) */
	if ((ret = getline(&line_buff, &line_buff_size, f)) != -1) {
		nb_all_elem = atoi(line_buff);
	}
	/* Allocate element2node array */
	if (m->elem2node) {
		free(m->elem2node);
	}

	/* Allocate elem2node using the value nb_elem read from the file.
	 * This value is the sum of all mixed element type present in the mesh.
	 * We will realloc elem2node later using the proper number of elements that we 
	 * need for computation. 
   */
	int nb_loc_nodes = NNO;
	m->elem2node = (int *)calloc(nb_loc_nodes * nb_all_elem, sizeof(int));
	assert(m->elem2node);

	int elem_number = 0;
	int elem_type = 0;
	int number_of_tags = 0;
	int node_number = 0;
	int countnode = 0;
	int countelem = 0;

	for (int i = 0; i < nb_all_elem; i++) {
		/* Field : elm-number */
		if ((ret = getdelim(&line_buff, &line_buff_size, delim, f)) != -1) {
			elem_number = atof(line_buff);
		}
		assert(elem_number - i == 1);

		/* Field : elm-type */
		if ((ret = getdelim(&line_buff, &line_buff_size, delim, f)) != -1) {
			elem_type = atof(line_buff);
		}

		/* Element is a T4 (GSMH_CODE = 4)*/
		if (elem_type == 4) {
			countelem++;

			/* Field number_of_tags */
			if ((ret = getdelim(&line_buff, &line_buff_size, delim, f)) != -1) {
				number_of_tags = atoi(line_buff);
			}

			/* First allocation of physical group */
			if ((m->has_physical_groups) && (!(m->physical_group_is_alloc))) {
				m->physical_group = (int *)calloc(nb_all_elem, sizeof(int));
				assert(m->physical_group);
				m->physical_group_is_alloc = true;
			}

			int *tags_val = (int *)calloc(number_of_tags, sizeof(int));

			for (int i_tag = 0; i_tag < number_of_tags; i_tag++) {
				if ((ret = getdelim(&line_buff, &line_buff_size, delim, f)) !=
					-1) {
					tags_val[i_tag] = atoi(line_buff);
				}
			}

			/* Store the physical group associated to the element */
			if (m->has_physical_groups) {
				m->physical_group[countelem - 1] = tags_val[0];
			}

			free(tags_val);
			tags_val = NULL;

			/* Field node-number-list ... */
			for (int j = 0; j < nb_loc_nodes - 1; j++) {
				/* Field : node-number (int) */
				if ((ret = getdelim(&line_buff, &line_buff_size, delim, f)) !=
					-1) {
					node_number = atoi(line_buff);
					m->elem2node[countnode] = node_number - 1;
					countnode++;
					assert(node_number - 1 <= m->nbnodes);
				}
			}

			/* Field : node-number (int) (Last one) */
			if ((ret = getline(&line_buff, &line_buff_size, f)) != -1) {
				node_number = atoi(line_buff);
				m->elem2node[countnode] = node_number - 1;
				countnode++;
				assert(node_number - 1 <= m->nbnodes);
			}
		} else {
			/* Skip non T4 elements lines */
			if ((ret = getline(&line_buff, &line_buff_size, f)) != -1) {
			}
		}
	}

	m->nbelems = countelem;
	m->elem2node = (int *)realloc(m->elem2node, countnode * sizeof(int));

	if ((m->has_physical_groups) && (m->physical_group_is_alloc)) {
		m->physical_group =
			(int *)realloc(m->physical_group, m->nbelems * sizeof(int));
	}

	free(line_buff);
	return 1;
}

void mesh_read_from_msh22_file(gdn_mesh *m, const char *filename)
{
	gdn_gdn_mesh_free(m);

	FILE *f = NULL;
	char *line_buff = NULL;
	size_t line_buff_size = 0;
	ssize_t ret = 1;
	printf("[gmsh_reader] : reading file       : %s\n", filename);

	/* Read the msh file */
	assert((f = fopen(filename, "r")));

	float version = 0;
	bool is_ascii = false;
	int data_size;

	/* Field : $MeshFormat */
	assert(gmsh_v2_extract_meshformat(&version, &is_ascii, &data_size, f));
	assert((version > 2) && (version < 3));
	assert(is_ascii);

	/* Field : $PhysicaNames (if exists) or $Nodes */
	if ((ret = getline(&line_buff, &line_buff_size, f)) != -1) {
		if (strcmp(line_buff, "$PhysicalNames\n") == 0) {
			m->has_physical_groups = true;
		} else if (strcmp(line_buff, "$Nodes\n") == 0) {
			m->has_physical_groups = false;
		} else {
			printf(
				"error: mesh_read_from_msh22_file: $PhysicalNames or $Nodes not located\n");
			exit(1);
		}
	}

	if (m->has_physical_groups) {
		assert(gmsh_v2_extract_physicalnames(m, f));
		/* Check end of $PhysicaNames block */
		if ((ret = getline(&line_buff, &line_buff_size, f)) != -1) {
			assert(strcmp(line_buff, "$EndPhysicalNames\n") == 0);
		}
		/* Move to the begining of  $Nodes block */
		if ((ret = getline(&line_buff, &line_buff_size, f)) != -1) {
			assert(strcmp(line_buff, "$Nodes\n") == 0);
		}
	}

	/* Field : $Nodes */
	assert(gmsh_v2_extract_nodes(m, f));
	/*Check end of $Nodes block */
	if ((ret = getline(&line_buff, &line_buff_size, f)) != -1) {
		assert(strcmp(line_buff, "$EndNodes\n") == 0);
	}

	/* Field : $Elements */
	if ((ret = getline(&line_buff, &line_buff_size, f)) != -1) {
		assert(strcmp(line_buff, "$Elements\n") == 0);
	}
	assert(gmsh_v2_extract_elements(m, f));

	/* Check end of $EndElements block */
	if ((ret = getline(&line_buff, &line_buff_size, f)) != -1) {
		assert(strcmp(line_buff, "$EndElements\n") == 0);
	}

	// 	/* NB : When physical entities are not explicitely defined
	//    * gmsh still puts a common identifier (1000) for volumic entities but does
	//    * not declare physical groups in the gmsh file header
	//    */
	// m->nb_phys_groups = 1;
	// m->phys_group_id = (int *)calloc(m->nb_phys_groups, sizeof(int));
	// m->phys_group_id[0] = m->physical_group[0];
	// m->phys_group_id_is_alloc = true;
	// for (int i = 1; i < m->nbelems; i++) {
	// 	int gid1 = m->physical_group[i];
	// 	bool is_new = true;

	// 	for (int k = 0; k < m->nb_phys_groups; k++) {
	// 		if (m->phys_group_id[k] == gid1) {
	// 			is_new = false;
	// 			break;
	// 		}
	// 	}
	// 	if (is_new) {
	// 		m->nb_phys_groups++;
	// 		m->phys_group_id =
	// 			(int *)realloc(m->phys_group_id, m->nb_phys_groups * sizeof(int));
	// 		m->phys_group_id[m->nb_phys_groups - 1] = m->physical_group[i];
	// 	}
	// }

	// printf("[Info - gmsh_reader()] : Found %d physical groups [",
	// 			 m->nb_phys_groups);
	// 	for (int k = 0; k < m->nb_phys_groups; k++) {
	// 		printf(" %d ", m->phys_group_id[k]);
	// 	}
	// 	printf("]\n");
	// }

	m->is_read = true;
	m->elem2elem = NULL;
	fclose(f);
	free(line_buff);
}

void gdn_gdn_mesh_print(gdn_mesh *m)
{
	int start = 0;
	printf("[gmsh_print]\n");
	printf("  #%-24s\n", "Nodes");
	printf("\t%-24s : %-12d\n", "m->nbnodes", m->nbnodes);
	printf("\t%-24s : %-12d\n", "m->nbqnodes_in", m->nbqnodes_in);

	if (m->qnode) {
		printf("\t%-12s %-12s %-12s %-12s\n", "index", "x", "y", "z");
		for (int i = 0; i < m->nbqnodes_in; i++) {
			printf("\t%-12d %-12f %-12f %-12f\n", i, m->qnode[3 * i + 0],
				   m->qnode[3 * i + 1], m->qnode[3 * i + 2]);
		}
	}

	int nb_loc_nodes = NQN;
	printf("\n  #%-24s\n", "Elements to qnodes");
	printf("\t%-24s : %-12d\n", "m->nbelems", m->nbelems);

	if (m->elem2qnode) {
		for (int i = 0; i < m->nbelems; i++) {
			printf("\t%-7d ", i + start);
			for (int j = 0; j < nb_loc_nodes; j++) {
				printf("%-5d", m->elem2qnode[nb_loc_nodes * i + j] + start);
			}
			printf("\n");
		}
	}

	printf("\n  #%-24s\n", "Elements to elements");
	if (m->elem2elem) {
		for (int i = 0; i < m->nbelems; i++) {
			printf("\t%-7d ", i + start);
			for (int j = 0; j < NFA; j++) {
				printf("%-5d ", m->elem2elem[NFA * i + j] + start);
			}
			printf("\n");
		}
	}

	printf("\n  #%-24s\n", "Faces to elements");
	printf("\t%-24s : %-12d\n", "m->nbfaces", m->nbfaces);

	char idL[6];
	char idR[6];

	if (m->face2elem) {
		for (int ifa = 0; ifa < m->nbfaces; ifa++) {
			sprintf(idL, "%s%d", "loc", m->face2elem[E2F * ifa + 1]);
			sprintf(idR, "%s%d", "loc", m->face2elem[E2F * ifa + 3]);
			printf("\t%-6d Left : %-6d %-6s Right : %-6d %-6s\n", ifa,
				   m->face2elem[E2F * ifa + 0] + start, idL,
				   m->face2elem[E2F * ifa + 2] + start, idR);
		}
	}
	printf("\n  #%-24s\n", "Faces to nodes");
	if (m->fnode2enode) {
		for (int ifa = 0; ifa < m->nbfaces; ifa++) {
			printf("\tifa %-6d\n", ifa);
			for (int ii = 0; ii < NQV; ii++) {
				int iL = m->fnode2enode[ifa * 2 * NQV + ii];
				int iR = m->fnode2enode[ifa * 2 * NQV + NQV + ii];
				printf("\t\t%-6d idL=%-6d iR=%-6d\n", ii, iL, iR);
			}
		}
	}

	/* List of nodes neighbours */
	printf("\n  #%-24s\n", "Node to elements");
	for (int ino = 0; ino < m->nbqnodes_in; ino++) {
		printf("\t%-7d ", ino);
		int ii = 0;
		while (m->node2elem[ii + m->max_node2elem * ino] != -1) {
			printf("%-6d ", m->node2elem[ii + m->max_node2elem * ino]);
			ii++;
			assert(ii < m->max_node2elem);
		}
		printf("\n");
	}
}

/* Fill array of faces of subcells */
void build_face(gdn_mesh *m, Face3 *face)
{
	assert(face);
	/* Loop over elements and their four faces */
	for (int ie = 0; ie < m->nbelems; ie++) {
		for (int ifa = 0; ifa < NFA; ifa++) {
			Face3 *f = face + NFA * ie + ifa;
			/* Using elem2node (read from .msh) we get the begining (NNO * ie) 
       * of the 4 nodes (NN0) compounding the current T4 element. 
       * We extract for each face (ifa) and in the proper 
       * order (t4_ref_face2node) the 3 nodes (ino) associated to the current 
       * face.
       */
			for (int ino = 0; ino < NVE; ino++)
				f->node[ino] =
					m->elem2node[NNO * ie +
								 t4_get_ref_face2node_entry(ifa, ino)];

			f->left = ie;
			f->locfaceleft = ifa;

			/* Sort the nodes (ascending order) of the current face */
			mesh_face_sort_node_face3(f);
		}
	}
	/* Sort the list of faces */
	qsort(face, NFA * m->nbelems, sizeof(Face3), mesh_face_compare_face3);
}

void build_elem2elem(gdn_mesh *m, Face3 *face)
{
	/* build_elem2elem : build the connectivity between elements.
   * One given element is connected to max 4 (NFA) other elements.
   */

	assert(m->elem2elem == NULL);
	m->elem2elem = (int *)calloc(NFA * m->nbelems, sizeof(int));
	assert(m->elem2elem);

	for (int i = 0; i < NFA * m->nbelems; i++)
		m->elem2elem[i] = -1;

	m->nbfaces = 0;
	int stop = NFA * m->nbelems;

	/* Loop over the sorted face */
	for (int ifa = 0; ifa < stop; ifa++) {
		Face3 *f1 = face + ifa;
		Face3 *f2 = face + ifa + 1;
		/* Two successive equal faces correspond to two neighbours elements */
		if (ifa != (stop - 1) && mesh_face_compare_face3(f1, f2) == 0) {
			ifa++;
			int ie1 = f1->left;
			int if1 = f1->locfaceleft;
			int ie2 = f2->left;
			int if2 = f2->locfaceleft;
			m->elem2elem[NFA * ie1 + if1] = ie2;
			m->elem2elem[NFA * ie2 + if2] = ie1;
		}
		m->nbfaces++;
	}
	// printf("nfaces=%d\n", m->nbfaces);
}

void build_boundaryinterface(gdn_mesh *m)
{
	/* build_boundaryinterface : extract boundary interface (face) id */
	m->nboundaryfaces = 0;
	m->boundaryface = NULL;
	/* Count faces lying on a boundary */
	for (int ifa = 0; ifa < m->nbfaces; ifa++) {
		if (m->face2elem[E2F * ifa + 2] == -1)
			m->nboundaryfaces += 1;
	}

	int idx = -1;
	if (m->nboundaryfaces > 0) {
		assert(m->boundaryface == NULL);
		m->boundaryface = (int *)calloc(m->nboundaryfaces, sizeof(int));
		assert(m->boundaryface);

		/* Boundary if element 2 id (neighboor) == -1 */
		for (int ifa = 0; ifa < m->nbfaces; ifa++) {
			if (m->face2elem[E2F * ifa + 2] == -1) {
				idx += 1;
				m->boundaryface[idx] = ifa;
			}
		}
	}
	// printf("m->nboundaryfaces: %d\n", m->nboundaryfaces);
}

void build_interface(gdn_mesh *m)
{
	/* build_interface : extract NON boundary interface (face) id */
	m->nmacrointerfaces = 0;
	m->macrointerface = NULL;
	/* Count faces NOT lying on a boundary */
	for (int ifa = 0; ifa < m->nbfaces; ifa++) {
		if (m->face2elem[E2F * ifa + 2] != -1)
			m->nmacrointerfaces += 1;
	}

	int idx = -1;
	if (m->nmacrointerfaces > 0) {
		assert(m->macrointerface == NULL);
		m->macrointerface = (int *)calloc(m->nmacrointerfaces, sizeof(int));
		assert(m->macrointerface);

		/* Interface if element 2 id (neighboor) is != -1 */
		for (int ifa = 0; ifa < m->nbfaces; ifa++) {
			if (m->face2elem[E2F * ifa + 2] != -1) {
				idx += 1;
				m->macrointerface[idx] = ifa;
			}
		}
	}
	// printf("m->nmacrointerfaces: %d\n", m->nmacrointerfaces);
}

void build_bounds(gdn_mesh *m)
{
	/* build_bounds : Find the coordinates of the minimal bounding box of the mesh */
	gdn_real xmin = m->node[0];
	gdn_real xmax = xmin;
	gdn_real ymin = m->node[1];
	gdn_real ymax = ymin;
	gdn_real zmin = m->node[2];
	gdn_real zmax = zmin;

	for (int i = 0; i < m->nbnodes; i++) {
		gdn_real x = m->node[3 * i];
		gdn_real y = m->node[3 * i + 1];
		gdn_real z = m->node[3 * i + 2];

		if (x < xmin)
			xmin = x;
		if (x > xmax)
			xmax = x;
		if (y < ymin)
			ymin = y;
		if (y > ymax)
			ymax = y;
		if (z < zmin)
			zmin = z;
		if (z > zmax)
			zmax = z;
	}

	m->xmin[0] = xmin;
	m->xmax[0] = xmax;
	m->xmin[1] = ymin;
	m->xmax[1] = ymax;
	m->xmin[2] = zmin;
	m->xmax[2] = zmax;

	// printf("bounds:  %.5e, %.5e, %.5e, %.5e, %.5e, %.5e\n",
	//        m->xmin[0], m->xmax[0], m->xmin[1], m->xmax[1], m->xmin[2], m->xmax[2]);
}

void build_face2elem(gdn_mesh *m)
{
	/* For each face we store : 
   * - Elem. 1 id 
   * - Face id of element 1 : (0, 1, 2, 3)
   * - Element 2 id  (-1 if boundary)
   * - Face id of element 2 : (0, 1, 2, 3)  (-1 if boundary)
   */
	assert(m->face2elem == NULL);
	m->face2elem = (int *)calloc(E2F * m->nbfaces, sizeof(int));
	assert(m->face2elem);

	int facecount = 0;
	for (int ie = 0; ie < m->nbelems; ie++) {
		for (int ifa = 0; ifa < NFA; ifa++) {
			int ie2 = m->elem2elem[NFA * ie + ifa];

			/* Deal only ONCE with the dependency ie->ieX in m->elem2elem and not ieX->ie */
			if (ie2 < ie) {
				m->face2elem[E2F * facecount + 0] =
					ie; /* Element 1 id          */
				m->face2elem[E2F * facecount + 1] =
					ifa; /* Face id of element 1  */
				m->face2elem[E2F * facecount + 2] =
					ie2; /* Element 2 id          */

				/* Search for the corresponding face id of element 2 */
				if (ie2 >= 0) {
					m->face2elem[E2F * facecount + 3] = -1;
					for (int ifa2 = 0; ifa2 < NFA; ifa2++) {
						if (m->elem2elem[NFA * ie2 + ifa2] == ie) {
							m->face2elem[E2F * facecount + 3] = ifa2;
							break;
						}
					}
					assert(m->face2elem[E2F * facecount + 3] != -1);
					/* Boundary face */
				} else {
					m->face2elem[E2F * facecount + 3] = -1;
				}
				facecount++;
			}
		}
	}
	assert(facecount == m->nbfaces);
}

void build_qnodes(gdn_mesh *m)
{
	int lst_size = NQN * m->nbelems;
	Node *node = malloc(lst_size * sizeof(Node));

	/* Get t10 nodes list (in element accessing order) and sort it (ascending) in 
   * lexicographic coordinates order
   */
	gdn_node_get_t10_elem_order(m, node);

	qsort(node, lst_size, sizeof(Node), mesh_node_compare_coordinate);

	/* Extract unique quadratic nodes in node_set and count and store the number 
  * of copy of each node 
  */

	int n_copy = 1;
	int idx_u = 0; /* Insertion index of the new unique quadratic node */
	int max_copy = -1; /* Max number of elem per node */
	Node *node_set = malloc(lst_size * sizeof(Node));

	node_set[0] = node[0];
	for (int k = 1; k < NQN * m->nbelems; k++) {
		node[k].rank_c = k;

		if (mesh_node_compare_coordinate(&node_set[idx_u], &node[k]) != 0) {
			idx_u++;
			node_set[idx_u] = node[k];
			/* Max number of element having the same node */
			if (n_copy > max_copy)
				max_copy = n_copy;
			n_copy = 1;
		} else {
			n_copy++;
			node[k].n_copy = n_copy;
		}
		node[k].rank_e = -1; /* Clean all rank_e to overwrite it later */
	}

	int nbqnode = idx_u + 1; /* To count the first insertion */

	/* Sort (ascending) unique quadratic node according to id_elem */
	qsort(node_set, nbqnode, sizeof(Node), mesh_node_compare_elem_id);

	/* Insert in nodes the final rank_e for the unique nodes  */
	for (int k = 0; k < nbqnode; k++)
		node[node_set[k].rank_c].rank_e = k;

	/* Propagate rank_e of unique node to duplicated nodes  */
	int id = 0;
	for (int k = 1; k < NQN * m->nbelems; k++) {
		if (mesh_node_compare_coordinate(&node[id], &node[k]) == 0) {
			node[k].rank_e = node[id].rank_e;
		} else {
			id = k;
		}
	}

	/* Reorder according to element */
	qsort(node, lst_size, sizeof(Node), mesh_node_compare_elem_id);

	m->nbqnodes_in = nbqnode;

	/* Remove T4 m->node */
	if (m->node) {
		free(m->node);
		m->node = NULL;
	}
	assert(m->node == NULL);
	m->node = (gdn_real *)calloc(NQN * 3 * m->nbqnodes_in, sizeof(gdn_real));
	assert(m->node);

	/* Fill m->node with T10 */
	for (int k = 0; k < m->nbqnodes_in; k++) {
		m->node[3 * k + 0] = node_set[k].x[0];
		m->node[3 * k + 1] = node_set[k].x[1];
		m->node[3 * k + 2] = node_set[k].x[2];
	}

	if (m->elem2qnode) {
		free(m->elem2node);
		m->elem2qnode = NULL;
	}
	assert(m->elem2qnode == NULL);
	m->elem2qnode = (int *)calloc(NQN * m->nbelems, sizeof(int));
	assert(m->elem2qnode);

	/* Fill m->elem2qnode with T10 */
	for (int ie = 0; ie < m->nbelems; ie++) {
		for (int iqn = 0; iqn < NQN; iqn++) {
			int id_node = node[NQN * ie + iqn].rank_e;
			m->elem2qnode[NQN * ie + iqn] = id_node;
		}
	}

	free(node);
	free(node_set);
}

/**Build face node to element node 
 * m->fnode2enode[k] = {L0, L1, L2, L3, L4, L5, R0, R1, R2, R3, R4, R5}
 * Interger Li (resp. Ri) store for i = 0..5 (nodes on the face) the local T10 
 * node id (0, 1, 2,...9) in the left (resp. right) element.
 */
void build_fnode2enode(gdn_mesh *m)
{
	assert(m->fnode2enode == NULL);
	m->fnode2enode = (int *)calloc(NQV * 2 * m->nbfaces, sizeof(int));
	assert(m->fnode2enode);

	/* Extend element2qnode with bd faces */
	m->elem2qnode = (int *)realloc(
		m->elem2qnode,
		(NQN * m->nbelems + NQV * m->nboundaryfaces) * sizeof(int));

	int facount = -1;

	for (int ifa = 0; ifa < m->nbfaces; ifa++) {
		int ieL = m->face2elem[E2F * ifa + 0];
		int ifL = m->face2elem[E2F * ifa + 1];
		int ieR = m->face2elem[E2F * ifa + 2];

		assert(ieR == m->elem2elem[NFA * ieL + ifL]);

		if (ieR < 0) {
			facount++;
		}

		for (int iloc = 0; iloc < NQV; iloc++) {
			/* Get the correspounding left local T10 node id */
			int inoe = t10_get_ref_face2node_entry(ifL, iloc);
			m->fnode2enode[ifa * NQV * 2 + iloc] = inoe;

			/* On a boundary face we store iloc as the right local T10 node id.
       * To get a corresponding global T10 node id for iloc, we append at the 
       * end of elem2qnode the global T10 node left element
       */
			if (ieR < 0) {
				m->fnode2enode[ifa * NQV * 2 + NQV + iloc] = iloc;
				m->elem2qnode[m->nbelems * NQN + facount * NQV + iloc] =
					m->elem2qnode[ieL * NQN + inoe];
			} else {
				/* Get the left global T10 node id */
				int ino = m->elem2qnode[ieL * NQN + inoe];
				int icount = 0;

				/* Search for the same global T10 node id in the right element */
				while ((icount < NQN) &&
					   (m->elem2qnode[ieR * NQN + icount] != ino)) {
					icount++;
				}
				assert(icount < NQN);
				m->fnode2enode[ifa * NQV * 2 + NQV + iloc] = icount;
			}
		}
	}
	assert(facount + 1 == m->nboundaryfaces);
}

// Build the node to elems connectivity
void build_node2elem(gdn_mesh *m)
{
	int *count = (int *)calloc(m->nbqnodes_in, sizeof(int));
	m->max_node2elem = 0;

	for (int ie = 0; ie < m->nbelems; ie++) {
		for (int iloc = 0; iloc < NQN; iloc++) {
			count[m->elem2qnode[iloc + NQN * ie]]++;
		}
	}

	for (int ino = 0; ino < m->nbqnodes_in; ino++) {
		m->max_node2elem =
			m->max_node2elem > count[ino] ? m->max_node2elem : count[ino];
	}
	free(count);

	// add a column of -1's for marking the ends of neighbours list
	// printf("max number of elems touching a node: %d\n", m->max_node2elem);
	(m->max_node2elem)++;

	m->node2elem =
		(int *)malloc((m->max_node2elem + 1) * m->nbqnodes_in * sizeof(int));
	assert(m->node2elem);

	// fill the array with -1's for marking the end of neighbours
	for (int i = 0; i < m->max_node2elem * m->nbqnodes_in; i++) {
		m->node2elem[i] = -1;
	}

	for (int ie = 0; ie < m->nbelems; ie++) {
		for (int iloc = 0; iloc < NQN; iloc++) {
			int ino = m->elem2qnode[iloc + NQN * ie];
			int ii = 0;

			/* Find index ii where to insert ie */
			while (m->node2elem[m->max_node2elem * ino + ii] != -1) {
				ii++;
			}
			assert(ii < m->max_node2elem);
			if (ii < m->max_node2elem - 1) {
				m->node2elem[m->max_node2elem * ino + ii] = ie;
			}
		}
	}

	// send to infinity nodes that do not belong to any element
	for (int ino = 0; ino < m->nbqnodes_in; ino++) {
		if (m->node2elem[0 + m->max_node2elem * ino] == -1) {
			m->qnode[0 + 3 * ino] = 1e10;
			// printf("WARNING : nodes at infinty\n");
		}
	}
}

// Build other connectivity arrays
void mesh_build_connectivity(gdn_mesh *m)
{
	assert(m->is_read);
	assert(!m->is_built);
	// clock_t tic;
	// clock_t toc;
	// double time;

	// printf("Build connectivity...\n");
	build_bounds(m);

	Face3 *face = malloc(NFA * m->nbelems * sizeof(Face3));
	build_face(m, face);
	build_elem2elem(m, face);
	free(face);

	build_face2elem(m);

	build_boundaryinterface(m);
	build_interface(m);

	// tic = clock();
	build_qnodes(m);
	// toc = clock();
	// time = (double)(toc - tic) / CLOCKS_PER_SEC;
	// printf("[Profiling time] =  build_qnode_2 : %f\n", time);

	build_fnode2enode(m);
	build_node2elem(m);

	if (m->elem2node) {
		free(m->elem2node);
		m->elem2node = NULL;
	}
	m->is_built = true;
}