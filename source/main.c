/*
 * main.c
 *
 * Main function
 *
 * Felipe Contreras
 * felipe.contrerass@postgrado.uv.cl
 */

/*
 * Copyright(c) 2022 Felipe Contreras
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/** ! testing
* Better comments API
* 
* Mymethod
*
*
* // Underline
* * Important information is highlighted
* ! Deprecated method, do not use
* ? Should this method be exposed in the public API?
* TODO: refactor this method so that it conforms to the API
* @param myParam The parameter for this method
 */

/** ! testign */




#include "main.h"

void grid_check()
{

	vtype aux_pot;
	vtype aux_dens;
	vtype aux_ax;
	vtype aux_ay;
	vtype aux_az;
	struct node *ptr_node;
	int no_pts;
	for (int i = 0; i < GL_tentacles_level_max + 1; i++)
	{
		no_pts = GL_tentacles_size[i];
		aux_pot = 0;
		aux_dens = 0;
		aux_ax = 0;
		aux_ay = 0;
		aux_az = 0;
		for (int j = 0; j < no_pts; j++)
		{
			ptr_node = GL_tentacles[i][j];
			for (int k = 0; k < ptr_node->grid_intr_size; k++)
			{
				aux_pot += ptr_node->ptr_pot[ptr_node->ptr_grid_intr[k]];
				aux_dens += ptr_node->ptr_d[ptr_node->ptr_grid_intr[k]];
				aux_ax += ptr_node->ptr_ax[ptr_node->ptr_grid_intr[k]];
				aux_ay += ptr_node->ptr_ay[ptr_node->ptr_grid_intr[k]];
				aux_az += ptr_node->ptr_az[ptr_node->ptr_grid_intr[k]];
			}
		}
		printf("\nInterior, lv = %d\n", i);
		//printf("pot = %.16Lf, d = %.16Lf, ax = %.16Lf, ay = %.16Lf, az = %.16Lf\n", aux_pot, aux_dens, aux_ax, aux_ay, aux_az);
		//printf("pot = %.16f, d = %.16f, ax = %.16f, ay = %.16f, az = %.16f\n", aux_pot, aux_dens, aux_ax, aux_ay, aux_az);
		printf("pot = %.10f, d = %.10f, ax = %.10f, ay = %.10f, az = %.10f\n", (double)aux_pot, (double)aux_dens, (double)aux_ax, (double)aux_ay, (double)aux_az);

		aux_pot = 0;
		aux_dens = 0;
		aux_ax = 0;
		aux_ay = 0;
		aux_az = 0;
		for (int j = 0; j < no_pts; j++)
		{
			ptr_node = GL_tentacles[i][j];
			for (int k = 0; k < ptr_node->grid_bder_size; k++)
			{
				aux_pot += ptr_node->ptr_pot[ptr_node->ptr_grid_bder[k]];
				aux_dens += ptr_node->ptr_d[ptr_node->ptr_grid_bder[k]];
				aux_ax += ptr_node->ptr_ax[ptr_node->ptr_grid_bder[k]];
				aux_ay += ptr_node->ptr_ay[ptr_node->ptr_grid_bder[k]];
				aux_az += ptr_node->ptr_az[ptr_node->ptr_grid_bder[k]];
			}
		}
		printf("\nBorder, lv = %d\n", i);
		//printf("pot = %.16Lf, d = %.16Lf, ax = %.16Lf, ay = %.16Lf, az = %.16Lf\n", aux_pot, aux_dens, aux_ax, aux_ay, aux_az);
		//printf("pot = %.16f, d = %.16f, ax = %.16f, ay = %.16f, az = %.16f\n", aux_pot, aux_dens, aux_ax, aux_ay, aux_az);
		printf("pot = %.10f, d = %.10f, ax = %.10f, ay = %.10f, az = %.10f\n", (double)aux_pot, (double)aux_dens, (double)aux_ax, (double)aux_ay, (double)aux_az);
	}
}


int main(int argc, char **argv) {

	printf("\n ********** RUNNING N-BODY CODE ********** \n\n");


	//printf(std::scientific);
	//printf(std::setprecision(digits_precision));

	//* DEFINING OTHER VARIABLES
	//** >> Defining timestep **/
	vtype dt = 0;
	vtype actualtime = 0;
	int Number_timesteps = 0;
	int Number_outputs = 0;
	/**
	 *! ******************************************************************************
	 *! ******************************************************************************
	 *! ******************************************************************************
	 */

	// OBSERVABLES:
	//   Energies[0] = Kinetic
	//   Energies[1] = Potential
	//   Energies[2] = Total = K + P
	vtype *energies;
	energies = (vtype *)calloc(3, sizeof(vtype));

	//** >> GLOBAL VARIABLES **/
	//printf("Global variables\n");
	GL_clock_begin = clock();
	global_variables();
    GL_times[0] += (double) (clock() - GL_clock_begin) / CLOCKS_PER_SEC;

	//** >> INPUT **/
	//printf("Input\n");
	GL_clock_begin = clock();
	if (input() == _FAILURE_){
		printf("\n\n Error running input() function\n\n");
		return _FAILURE_;
	}
    GL_times[1] += (double) (clock() - GL_clock_begin) / CLOCKS_PER_SEC;

	//** >> INITIALIZATION **/
	//printf("Initialization\n");
	GL_clock_begin = clock();
	if (initialization() == _FAILURE_)
	{
		printf("\n\n Error running initialization() function\n\n");
		return _FAILURE_;
	}
	GL_times[2] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

	//** >> BUILDING INITIAL TREE **/
	//printf("Building initial Tree\n");
	GL_clock_begin = clock();
	if (tree_construction() == _FAILURE_)
	{
		printf("\n\n Error running tree_construction() function\n\n");
		return _FAILURE_;
	}
	GL_times[3] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

	//** >> DENSITY COMPUTATION **/
	//printf("Initial density\n");
	GL_clock_begin = clock();
	grid_density();
	GL_times[4] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

	//** >> POTENTIAL COMPUTATION **/
	//printf("Initial Potential\n");
	GL_clock_begin = clock();
	if (potential() == _FAILURE_)
	{
		printf("\n\n Error running potential() function\n\n");
		return _FAILURE_;
	}
	GL_times[5] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

	//** >> GRID ACCELERATION **/
	//printf("Initial Grid acceleration\n");
	GL_clock_begin = clock();
	if (grid_acceleration() == _FAILURE_)
	{
		printf("\n\n Error running grid_acceleration() function\n\n");
		return _FAILURE_;
	}
	GL_times[6] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

	//** >> PARTICLE ACCELERATION **/
	//printf("Initial Particle Acceleration\n");
	GL_clock_begin = clock();
	if (particle_acceleration() == _FAILURE_)
	{
		printf("\n\n Error running particle_acceleration() function\n\n");
		return _FAILURE_;
	}
	GL_times[7] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

	//** >> OBSERVABLES **/
	//printf("Initial Observables\n");
	GL_clock_begin = clock();
	observables(energies);
	GL_times[13] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

	//** >> OUPUT SNAPSHOTS **/
	//printf("Initial Snapshot output\n");
	GL_clock_begin = clock();
	if (output_snapshots(energies, actualtime, Number_outputs) == _FAILURE_)
	{
		printf("\n\n Error running output_snapshots() function\n\n");
		return _FAILURE_;
	}
	GL_times[19] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

	while (actualtime < Maxdt && Number_timesteps < MaxIterations)
	{
		
		// printf("\n Inicio de while\n");
		// for (int i = 0; i < 10000; i++)
		// {
		// 	printf("%d, x = %.16f, y = %.16f, z = %.16f\n", i, (double)GL_ptcl_x[i], (double)GL_ptcl_y[i], (double)GL_ptcl_z[i]);
		// 	printf("%d, ax = %.16f, ay = %.16f, az = %.16f\n", i, (double)GL_ptcl_ax[i], (double)GL_ptcl_ay[i], (double)GL_ptcl_az[i]);
		// }
		// grid_check();

		//** >> TIMESTEP **/
		// printf("Time-step compute\n");
		GL_clock_begin = clock();
		if (timestep_2(&dt) == _FAILURE_)
		{
			printf("\n\n Error running timestep_compute() function\n\n");
			return _FAILURE_;
		}
		GL_times[8] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

		actualtime += dt; // Updating actual time

		//** >> PARTICLE UPDATING A **/
		// printf("Particles Updating A\n");
		GL_clock_begin = clock();
		if (particle_updating_A(dt) == _FAILURE_)
		{
			printf("\n\n Error running particle_updating() function\n\n");
			return _FAILURE_;
		}
		GL_times[9] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

		//printf("\n AFTER Particle Updating A\n");
		// for (int i = 0; i < 10000; i++)
		// {
		// 	printf("%d, x = %.16f, y = %.16f, z = %.16f\n", i, GL_ptcl_x[i], GL_ptcl_y[i], GL_ptcl_z[i]);
		// 	printf("%d, ax = %.16f, ay = %.16f, az = %.16f\n", i, GL_ptcl_ax[i], GL_ptcl_ay[i], GL_ptcl_az[i]);
		// }
		//grid_check();

		//** >> RESET **/
		// printf("Reset\n");
		GL_clock_begin = clock();
		if (reset() == _FAILURE_)
		{
			printf("\n\n Error running initialization() function\n\n");
			return _FAILURE_;
		}
		GL_times[11] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

		// printf("tree adaptation\n");
		//** >> TREE ADAPTATION **/
		GL_clock_begin = clock();
		if (tree_adaptation() == _FAILURE_)
		{
			printf("\n\n Error running tree_construction() function\n\n");
			return _FAILURE_;
		}
		GL_times[10] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

		//printf("\n AFTER Tree adaptation \n");
		// for (int i = 0; i < 10000; i++)
		// {
		// 	printf("%d, x = %.16f, y = %.16f, z = %.16f\n", i, GL_ptcl_x[i], GL_ptcl_y[i], GL_ptcl_z[i]);
		// 	printf("%d, ax = %.16f, ay = %.16f, az = %.16f\n", i, GL_ptcl_ax[i], GL_ptcl_ay[i], GL_ptcl_az[i]);
		// }
		//grid_check();

		// printf("density \n");
		//** >> DENSITY COMPUTATION **/
		GL_clock_begin = clock();
		grid_density();
		GL_times[4] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

		//printf("\n AFTER density \n");
		// for (int i = 0; i < 10000; i++)
		// {
		// 	printf("%d, x = %.16f, y = %.16f, z = %.16f\n", i, GL_ptcl_x[i], GL_ptcl_y[i], GL_ptcl_z[i]);
		// 	printf("%d, ax = %.16f, ay = %.16f, az = %.16f\n", i, GL_ptcl_ax[i], GL_ptcl_ay[i], GL_ptcl_az[i]);
		// }
		//grid_check();

		// printf("Potential\n");
		//** >> POTENTIAL COMPUTATION **/
		GL_clock_begin = clock();
		if (potential() == _FAILURE_)
		{
			printf("\n\n Error running potential() function\n\n");
			return _FAILURE_;
		}
		GL_times[5] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

		// printf("Grid acceleration\n");
		//** >> GRID ACCELERATION **/
		GL_clock_begin = clock();
		if (grid_acceleration() == _FAILURE_)
		{
			printf("\n\n Error running grid_acceleration() function\n\n");
			return _FAILURE_;
		}
		GL_times[6] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

		//printf("\n AFTER Grid acceleration \n");
		// for (int i = 0; i < 10000; i++)
		// {
		// 	printf("%d, x = %.16f, y = %.16f, z = %.16f\n", i, GL_ptcl_x[i], GL_ptcl_y[i], GL_ptcl_z[i]);
		// 	printf("%d, ax = %.16f, ay = %.16f, az = %.16f\n", i, GL_ptcl_ax[i], GL_ptcl_ay[i], GL_ptcl_az[i]);
		// }
		//grid_check();

		//** >> PARTICLE ACCELERATION **/
		// printf("Particle acceleration\n");
		GL_clock_begin = clock();
		if (particle_acceleration() == _FAILURE_)
		{
			printf("\n\n Error running particle_acceleration() function\n\n");
			return _FAILURE_;
		}
		GL_times[7] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

		//printf("\n AFTER Particle acceleration \n");
		// for (int i = 0; i < 10000; i++)
		// {
		// 	printf("%d, x = %.16f, y = %.16f, z = %.16f\n", i, GL_ptcl_x[i], GL_ptcl_y[i], GL_ptcl_z[i]);
		// 	printf("%d, ax = %.16f, ay = %.16f, az = %.16f\n", i, GL_ptcl_ax[i], GL_ptcl_ay[i], GL_ptcl_az[i]);
		// }
		//grid_check();

		//** >> PARTICLE UPDATING B **/
		// printf("Particle Updating B\n");
		GL_clock_begin = clock();
		if (particle_updating_B(dt) == _FAILURE_)
		{
			printf("\n\n Error running particle_updating_B() function\n\n");
			return _FAILURE_;
		}
		GL_times[12] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

		//** >> OUPUT SNAPSHOTS **/

		if (Number_timesteps % fr_output == 0 || actualtime >= Maxdt)
		{
			++Number_outputs;

			//** >> OBSERVABLES **/
			// printf("Observables\n");
			GL_clock_begin = clock();
			observables(energies);
			GL_times[13] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

			//** >> OUPUT SNAPSHOTS **/
			// printf("Output snapshots\n");
			GL_clock_begin = clock();
			if (output_snapshots(energies, actualtime, Number_outputs) == _FAILURE_)
			{
				printf("\n\n Error running output_snapshots() function\n\n");
				return _FAILURE_;
			}
			GL_times[19] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;
		}

		// BAR PROGRESS
		fflush(stdout);
		printf("\r[");
		for (int j = 0; j < 50; j++)
		{
			if (j < (int)(50 * actualtime / Maxdt))
				printf("=");
			else if (j == (int)(50 * actualtime / Maxdt))
				printf(">");
			else
				printf(" ");
		}
		printf("] %.2f %%", (double)(actualtime * 100.0 / Maxdt));

		// Particles printf
		// if (Number_timesteps % 500 == 0)
		// {
		// 	printf("\n Fin while\n");
		// 	for (int i = 0; i < 10000; i++)
		// 	{
		// 		printf("%d, x = %.10f, y = %.10f, z = %.10f\n", i, (double)GL_ptcl_x[i], (double)GL_ptcl_y[i], (double)GL_ptcl_z[i]);
		// 		printf("%d, ax = %.10f, ay = %.10f, az = %.10f\n", i, (double)GL_ptcl_ax[i], (double)GL_ptcl_ay[i], (double)GL_ptcl_az[i]);
		// 	}
		// 	grid_check();
		// }


		Number_timesteps++; // Increasing the number of time-steps
	}

	printf("\n\n%sMaxdt = %f, Final time = %f%s\n\n",KRED ,(double) Maxdt, (double) actualtime,KNRM);

	/**
	*! ******************************************************************************
	*! ******************************************************************************
	*! ****************************************************************************** 
	*/

	//** >> OUPUT MAIN PARAMETERS **/
	GL_clock_begin = clock();
	output_main_parameters(actualtime, Number_timesteps, Number_outputs);
	GL_times[18] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

	//** >> TERMINAL PRINT **/
	terminal_print();

	printf("Total number of time-steps = %d\n", Number_timesteps);

	printf("\n\n ******** FINALIZED N-BODY CODE ********** \n\n");

	return _SUCCESS_;
}
