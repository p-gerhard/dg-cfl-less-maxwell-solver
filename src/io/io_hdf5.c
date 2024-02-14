/*
 * Copyright (c) 
 * 2019-2021 Philippe Helluy <helluy@math.unistra.fr>
 * 2020-2021 Pierre Gerhard <pierre.gerhard@gmail.com>
 *
 * gdon3d is free software; you can redistribute it and/or modify
 * it under the terms of the GPLv3 license. See LICENSE for details.
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <string.h>

#include <hdf5.h>

#include <mesh/t10.h>

#include <gdon3d.h>
#include <simulation/simulation.h>
#include "io_hdf5.h"

void set_name_mesh(char *filename, char *basename)
{
	int len_base = strlen(basename);
	int len_ext = strlen("_mesh.h5");
	char *locfilename = (char *)malloc((len_base + len_ext + 1) * sizeof(char));

	sprintf(locfilename, "%s_mesh.h5", basename);
	strcpy(filename, locfilename);

	free(locfilename);
}

void set_name_data(char *filename, char *basename, int iter_id)
{
	int len_base = strlen(basename);
	int len_ext = strlen("_data_0000");
	char *locfilename = (char *)malloc((len_base + len_ext + 1) * sizeof(char));

	sprintf(locfilename, "%s_data_%04d", basename, iter_id);
	strcpy(filename, locfilename);

	free(locfilename);
}

void set_name_xmf(char *filename, char *basename, int iter_id)
{
	int len_base = strlen(basename);
	int len_ext = strlen("_simulation_0000.xmf");
	char *locfilename = (char *)malloc((len_base + len_ext + 1) * sizeof(char));

	sprintf(locfilename, "%s_simulation_%04d.xmf", basename, iter_id);
	strcpy(filename, locfilename);

	free(locfilename);
}

void io_hdf5_dump_mesh(gdn_simulation *simu, char *geo_filename)
{
	/* HDF5 handlers */
	hid_t file_id, dataset_id, dataspace_id; /* Identifiers */
	hsize_t dims1[1]; /* Small scalar data */
	hsize_t dims_conn[2]; /* Connectivity */
	hsize_t dims_nodes[2]; /* Nodes */
	herr_t status __attribute__((unused));

	/* Create file */
	file_id = H5Fcreate(geo_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	/* Small scalar data */
	dims1[0] = 1;
	int data_i[1];
	dataspace_id = H5Screate_simple(1, dims1, NULL);

	data_i[0] = simu->tt->nbelems;
	dataset_id = H5Dcreate2(file_id, "nb_elems", H5T_NATIVE_INT, dataspace_id,
							H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
					  data_i);

	status = H5Dclose(dataset_id);
	status = H5Sclose(dataspace_id);

	/* Connectivity */
	dims_conn[0] = simu->tt->nbelems;
	dims_conn[1] = 10;

	int *conn_buffer =
		(int *)malloc(dims_conn[1] * simu->tt->nbelems * sizeof(int));

	for (unsigned int ie = 0; ie < simu->tt->nbelems; ie++) {
		for (unsigned int in = 0; in < dims_conn[1]; in++) {
			conn_buffer[ie * dims_conn[1] + in] =
				simu->tt->elem2qnode[ie * 10 + in];
		}
	}

	char dspace_name1[sizeof("connectivity_000000")];
	sprintf(dspace_name1, "connectivity_%06d", 0);

	dataspace_id = H5Screate_simple(2, dims_conn, NULL);

	dataset_id = H5Dcreate2(file_id, dspace_name1, H5T_NATIVE_INT, dataspace_id,
							H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
					  conn_buffer);

	status = H5Dclose(dataset_id);
	status = H5Sclose(dataspace_id);

	free(conn_buffer);

	/* Nodes */
	int nb_points = simu->tt->nbqnodes_in;
	dims_nodes[0] = nb_points;
	dims_nodes[1] = 3;

	double *nodes_buffer = (gdn_real *)malloc(3 * nb_points * sizeof(gdn_real));

	int ipoint = 0;

	for (unsigned int in = 0; in < nb_points; in++) {
		for (unsigned int id = 0; id < 3; id++) {
			nodes_buffer[ipoint + id] = simu->tt->node[3 * in + id];
		}
		ipoint += 3;
	}

	char dspace_name2[sizeof("nodes_000000")];
	sprintf(dspace_name2, "nodes_%06d", 0);

	dataspace_id = H5Screate_simple(2, dims_nodes, NULL);

	dataset_id =
		H5Dcreate2(file_id, dspace_name2, H5T_NATIVE_DOUBLE, dataspace_id,
				   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
					  H5P_DEFAULT, nodes_buffer);
	status = H5Dclose(dataset_id);
	status = H5Sclose(dataspace_id);
	status = H5Fclose(file_id);
	free(nodes_buffer);
}

void io_hdf5_dump_data(gdn_simulation *simu, char *data_filename_prefix)
{
	/* HDF5 handlers */

	hid_t file_id, dataset_id, dataspace_id; /* Identifiers */
	hsize_t dims_data[1];
	herr_t status __attribute__((unused));

	int data_filename_len =
		strlen(data_filename_prefix) + 1;
	
	
	int ext_len = strlen(".hdf5");

	char *data_filename = (char *)malloc((data_filename_len + ext_len) * sizeof(char));

	snprintf(data_filename, data_filename_len + ext_len, "%s.hdf5", data_filename_prefix);
	printf("%s\n", data_filename);
	file_id = H5Fcreate(data_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	free(data_filename);

	/* Loop oover nodes */
	int nb_points = simu->tt->nbqnodes_in;
	int max_node2elem = simu->tt->max_node2elem;
	int m = simu->m;

	gdn_real *value = (gdn_real *)malloc(m * nb_points * sizeof(gdn_real));

	for (unsigned int in = 0; in < nb_points; in++) {
		int ie = simu->tt->node2elem[max_node2elem *
									 in]; // get the first element of the list
		int iloc = 0;
		//Get local node id
		while (simu->tt->elem2qnode[ie * NQN + iloc] != in) {
			iloc++;
		}

		unsigned int ipg = ie * NQN + iloc;
		if (simu->use_kinetic_scheme) {
			int neq = simu->neq;
			int nb_w = simu->mdl->nb_w;
			int nb_v = simu->mdl->nb_v;

			/* Method 1 */
			for (int im = 0; im < m; im++) {
				value[im * nb_points + in] = simu->wn[im * neq + ipg];
			}

			/* Methode 2 */
			// for (int iv = 0; iv < nb_v; iv++) {
			// 	for (int iw = 0; iw < nb_w; iw++) {
			// 		int offset_wn = iv * nb_w * neq + iw * neq;
			// 		value[iw * nb_v + iv * nb_points + in] = simu->wn[offset_wn + ipg];
			// 	}
			// }

		} else {
			for (unsigned int iv = 0; iv < m; iv++) {
				unsigned int imem = simulation_get_varindex(simu->m, ipg, iv);
				value[iv * nb_points + in] = simu->wn[imem];
			}
		}
	}

	dims_data[0] = nb_points;
	dataspace_id = H5Screate_simple(1, dims_data, NULL);

	int lset = sizeof("f00000");
	char *dset_name = (char *)malloc(lset * sizeof(char));

	for (int iv = 0; iv < m; iv++) {
		sprintf(dset_name, "f%05d", (unsigned short)iv);

		dataset_id =
			H5Dcreate2(file_id, dset_name, H5T_NATIVE_DOUBLE, dataspace_id,
					   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
						  H5P_DEFAULT, value + iv * nb_points);

		status = H5Dclose(dataset_id);
	}
	free(dset_name);

	status = H5Sclose(dataspace_id);

	free(value);

	status = H5Fclose(file_id);

}

void io_hdf5_write_xdmf(gdn_simulation *simu, char *filename,
						char *geo_filename, char *data_filename_prefix)
{
	int data_filename_len =
		strlen(data_filename_prefix) + 5;
	char *data_filename = (char *)malloc(data_filename_len * sizeof(char));

	FILE *xdmffile;
	xdmffile = fopen(filename, "w");

	fprintf(xdmffile, "<?xml version=\"1.0\" ?>\n");
	fprintf(xdmffile, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
	fprintf(
		xdmffile,
		"<Xdmf Version=\"2.0\" xmlns:xi=\"[http://www.w3.org/2001/XInclude]\">\n");
	fprintf(xdmffile, "<Domain>\n");
	fprintf(xdmffile, "<Grid Name=\"TheMesh\">\n");

	fprintf(xdmffile, "<Time Value=\"%f\"/>\n", simu->t);

	fprintf(xdmffile,
			"<Topology TopologyType = \"Tet_10\"  NumberOfElements =\"%d \">\n",
			simu->tt->nbelems);

	fprintf(xdmffile, "<DataItem Dimensions=\"%d 10\" Format=\"HDF\">\n",
			simu->tt->nbelems);
	fprintf(xdmffile, "%s", geo_filename);
	fprintf(xdmffile, "://connectivity_%06d\n", 0);
	fprintf(xdmffile, "</DataItem>\n");
	fprintf(xdmffile, "</Topology>\n");
	fprintf(xdmffile, "<Geometry Type=\"XYZ\" >\n");
	fprintf(
		xdmffile,
		"<DataItem Name =\"xyzvalues\" Format=\"hdf5\" NumberType=\"Float\" Dimensions=\"%d 3\">\n",
		simu->tt->nbqnodes_in);

	fprintf(xdmffile, "%s", geo_filename);
	fprintf(xdmffile, "://nodes_%06d\n", 0);
	fprintf(xdmffile, "</DataItem>\n");
	fprintf(xdmffile, "</Geometry>\n");

	/* Micro variables */

	if (simu->use_kinetic_scheme) {
		for (unsigned int iv = 0; iv < simu->m; iv++) {
			fprintf(xdmffile, "<Attribute Name=\"f%d\" center=\"node\" >\n",
					iv);
			fprintf(
				xdmffile,
				"<DataItem Format=\"hdf5\" NumberType=\"Float\" Dimensions=\"%d 1\">\n",
				simu->tt->nbqnodes_in);
			sprintf(data_filename, "%s.hdf5", data_filename_prefix);
			fprintf(xdmffile, "%s://f%05d\n", data_filename,
					(unsigned short)iv);
			fprintf(xdmffile, "</DataItem>\n");
			fprintf(xdmffile, "</Attribute>\n");
		}

		// computing macro variables
		for (int iw = 0; iw < simu->mdl->nb_w; iw++) {
			fprintf(xdmffile, "<Attribute Name=\"w%d\" center=\"node\" >\n",
					iw);
			fprintf(xdmffile, "<DataItem ItemType=\"Function\"\n");
			fprintf(xdmffile, "  Function=\"");
			int dim = simu->mdl->nb_v;
			for (int iv = 0; iv < dim; iv++) {
				if (iv > 0) {
					fprintf(xdmffile, "+");
				}
				fprintf(xdmffile, "$%d", iv);
			}
			fprintf(xdmffile, "\"\n");
			fprintf(xdmffile, "  Dimensions=\"%d 1\">\n",
					simu->tt->nbqnodes_in);

			for (int iv = 0; iv < dim; iv++) {
				fprintf(
					xdmffile,
					"  <DataItem Format=\"hdf5\" NumberType=\"Float\" Dimensions=\"%d 1\">\n",
					simu->tt->nbqnodes_in);
				sprintf(data_filename, "%s.hdf5", data_filename_prefix);

				fprintf(xdmffile, "%s://f%05d\n", data_filename, iv + iw * dim);
				fprintf(xdmffile, "  </DataItem>\n");
			}
			fprintf(xdmffile, "</DataItem>\n");
			fprintf(xdmffile, "</Attribute>\n");
		}
	} else {
		for (unsigned int iv = 0; iv < simu->m; iv++) {
			fprintf(xdmffile, "<Attribute Name=\"w%d\" center=\"node\" >\n",
					iv);
			fprintf(
				xdmffile,
				"<DataItem Format=\"hdf5\" NumberType=\"Float\" Dimensions=\"%d 1\">\n",
				simu->tt->nbqnodes_in);
			sprintf(data_filename, "%s.hdf5", data_filename_prefix);
			fprintf(xdmffile, "%s://f%05d\n", data_filename,
					(unsigned short)iv);
			fprintf(xdmffile, "</DataItem>\n");
			fprintf(xdmffile, "</Attribute>\n");
		}
	}

	fprintf(xdmffile, "</Grid>\n");
	fprintf(xdmffile, "</Domain>\n");
	fprintf(xdmffile, "</Xdmf>\n");
	free(data_filename);

	fclose(xdmffile);
}

void io_save_all(gdn_simulation *simulation, char filename[1024], int iternb)
{
	char mesh_filename[1024] = { 0 };
	char data_filename[1024] = { 0 };
	char xmf_filename[1024] = { 0 };
	set_name_mesh(mesh_filename, filename);
	set_name_data(data_filename, filename, iternb);
	set_name_xmf(xmf_filename, filename, iternb);

	/* Dump the mesh and free memory */
	io_hdf5_dump_mesh(simulation, mesh_filename);
	io_hdf5_dump_data(simulation, data_filename);
	io_hdf5_write_xdmf(simulation, xmf_filename, mesh_filename, data_filename);
}
