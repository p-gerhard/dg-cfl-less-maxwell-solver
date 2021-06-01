/*
 * Copyright (c) 
 * 2019-2021 Philippe Helluy <helluy@math.unistra.fr>
 * 2020-2021 Pierre Gerhard <pierre.gerhard@gmail.com>
 *
 * gdon3d is free software; you can redistribute it and/or modify
 * it under the terms of the GPLv3 license. See LICENSE for details.
 */
#ifndef IO_HDF5_H
#define IO_HDF5_H

#include <simulation/simulation.h>

//!\brief function that computes the name of the mesh file
//!\param[in] basename is the user name for current simulation
//! \returns the new name of mesh file
void set_name_mesh(char *filename, char *basename);

//!\brief function that computes the name of the data files
//!\param[in] basename is the user name for current simulation
//!\param[in] iternb is an index of the data file
//! \returns the new name of data file
void set_name_data(char *filename, char *basename, int iternb);

//!\brief function that computes the name of the xmf files
//!\param[in] basename is the user name for current simulation
//!\param[in] iternb is an index of the data file
//! \returns the new name of xmf file
void set_name_xmf(char *filename, char *basename, int iternb);

//!\brief write mesh static mesh data into hdf5 file (for use with  gdn_simulation_xdmf_plot_simulations_hdf5_structured_light)
//!\param[in] tf a tetrasimulation
//!\param[in] geo_filename : name of the geometry filename . Geometry/Topology has only have to be dumped once per simulation
//! as long all xmf files point to the same geo data file
void io_hdf5_dump_mesh(gdn_simulation *tf, char *geo_filename);

//!\brief write simulation hdf5 data for unstructured tetras mesh
//!\param[in] tf a tetrasimulation
//!\param[in] data_filename_prefix : self explanatory , the routine will append mpi node number and hdf5 extension to the prefix.

void io_hdf5_dump_data(gdn_simulation *tf, char *data_filename_prefix);


//!\brief write simulation hdf5 data for unstructured tetras mesh
//!\param[in] tf : a tetrasimulation
//!\param[in] geo_filename : name of the hdf5 file containing nodes/topology data
//!\param[in] data_filename_prefix : prefix used for hdf5 simulation data files
void io_hdf5_write_xdmf (gdn_simulation * tf, char *filename, char *geo_filename,
                                   char* data_filename_prefix);

void io_save_all(gdn_simulation *simulation, char filename[1024], int iternb);
#endif
