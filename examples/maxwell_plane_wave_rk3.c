#undef NDEBUG
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <time.h>

#include <dg/dg.h>
#include <models/model.h>
#include <models/maxwell.h>
#include <io/io_hdf5.h>
#include <simulation/simulation.h>
#include <utils/timing.h>
#include <utils/convergence.h>


static char base_name[256] = "maxwell_rk3_plane_wave";

/**
 * \fn imposed_macro_cos
 * \brief Imposed oscillatory solution of frequence nu for Maxwell's equations.
 *        This solution is constructed using the method of characteristics, 
 *        therby it is an exact solution for all (x, t) of Maxwell's equations.
 * \param[in] x    Position vector
 * \param[in] t    Time
 * \param[inout] w Macro state
 */
void static imposed_macro_cos(const gdn_real *x, const gdn_real t, gdn_real *w)
{
	// const gdn_real nu = 0.5;
	const gdn_real nu = 1;
	const gdn_real s = cos(nu * M_PI * (x[0] - t));
	// const gdn_real s = (x[0] - t) * (x[0] - t);
	w[0] = 0;
	w[1] = 0;
	w[2] = s;
	w[3] = 0;
	w[4] = -s;
	w[5] = 0;
}



int test_maxwell_rk3_plane_wave(void);

int test_maxwell_rk3_plane_wave(void)
{
	bool export_xdmf = true;

	gdn_mesh mesh = { 0 };
	gdn_simulation simu = { 0 };
	gdn_model model = { 0 };

	model.nb_w = 6;
	model.get_num_flux = maxwell_num_flux_upwind;
	model.get_num_flux_boundary = maxwell_num_flux_upwind;
	model.get_imposed_data = imposed_macro_cos;
	model.use_relax_scheme = false;
	model.lambda = 1.0;

	char *mesh_name = "../data/mesh/t4_cube_8.msh";

	/* Time settings */
	int mesh_raf = 9;
	gdn_real tmax = 1.0;
	gdn_real dt = 0.032 * (1. / mesh_raf);
	
	export_xdmf = true;

	utils_l2_time_convergence(&model, 4, mesh_name, dt, tmax, export_xdmf);
	
	return 1;
	
	
	struct timeval start;
	mesh_read_from_msh22_file(&mesh, mesh_name);
	mesh_build_connectivity(&mesh);
	simulation_init(&mesh, &model, &simu);
	simulation_set_time_parameters(&simu, tmax, dt, -1);
	simulation_display_info(&simu);
	tic(&start);
	dg_init_sol_macro(&simu);
	// dg_solve(&simu);
	toc(&start);

	/* Compute L2 error */
	gdn_real error_l2 = simulation_error_l2(&simu);
	simulation_dump_info(0, &simu, error_l2, 0);
  
  if (export_xdmf) {
	  mesh_raf = mesh_raf * 2;
	  int iternb = (int)floor(tmax / dt);
	  char filename[1024];
	  snprintf(filename, 1024, "%s", base_name);
	  io_save_all(&simu, filename, iternb);
  }
  simulation_free(&simu);
  return 1;
}

int main(void)
{
	int resu = test_maxwell_rk3_plane_wave();
	if (resu)
		printf("[Info] - test_maxwell_rk3_plane_wave : OK\n");
	else
		printf("[Info] -  test_maxwell_rk3_plane_wave : FAILED\n");
	return !resu;
}