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
	GL_clock_begin = clock();
	global_variables();
    GL_times[0] += (double) (clock() - GL_clock_begin) / CLOCKS_PER_SEC;

	//** >> INPUT **/
	GL_clock_begin = clock();
	if (input() == _FAILURE_){
		printf("\n\n Error running input() function\n\n");
		return _FAILURE_;
	}
    GL_times[1] += (double) (clock() - GL_clock_begin) / CLOCKS_PER_SEC;

	//** >> INITIALIZATION **/
	GL_clock_begin = clock();
	if (initialization() == _FAILURE_)
	{
		printf("\n\n Error running initialization() function\n\n");
		return _FAILURE_;
	}
	GL_times[2] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

	//** >> BUILDING INITIAL TREE **/
	GL_clock_begin = clock();
	if (tree_construction() == _FAILURE_)
	{
		printf("\n\n Error running tree_construction() function\n\n");
		return _FAILURE_;
	}
	GL_times[3] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

	//** >> DENSITY COMPUTATION **/
	GL_clock_begin = clock();
	grid_density();
	GL_times[4] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

	//** >> POTENTIAL COMPUTATION **/
	GL_clock_begin = clock();
	if (potential() == _FAILURE_)
	{
		printf("\n\n Error running potential() function\n\n");
		return _FAILURE_;
	}
	GL_times[5] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

	//** >> GRID ACCELERATION **/
	GL_clock_begin = clock();
	if (grid_acceleration() == _FAILURE_)
	{
		printf("\n\n Error running grid_acceleration() function\n\n");
		return _FAILURE_;
	}
	GL_times[6] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

	//** >> PARTICLE ACCELERATION **/
	GL_clock_begin = clock();
	if (particle_acceleration() == _FAILURE_)
	{
		printf("\n\n Error running particle_acceleration() function\n\n");
		return _FAILURE_;
	}
	GL_times[7] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

	//** >> OBSERVABLES **/
	GL_clock_begin = clock();
	observables(energies);
	GL_times[13] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

	//** >> OUPUT SNAPSHOTS **/
	GL_clock_begin = clock();
	if (output_snapshots(energies, actualtime, Number_outputs) == _FAILURE_)
	{
		printf("\n\n Error running output_snapshots() function\n\n");
		return _FAILURE_;
	}
	GL_times[19] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

	while (actualtime < Maxdt && Number_timesteps < MaxIterations)
	{
		//** >> TIMESTEP **/
		GL_clock_begin = clock();
		if (timestep_2(&dt) == _FAILURE_)
		{
			printf("\n\n Error running timestep_compute() function\n\n");
			return _FAILURE_;
		}
		GL_times[8] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

		actualtime += dt; // Updating actual time
		
		//** >> PARTICLE UPDATING A **/
		GL_clock_begin = clock();
		if (particle_updating_A(dt) == _FAILURE_)
		{
			printf("\n\n Error running particle_updating() function\n\n");
			return _FAILURE_;
		}
		GL_times[9] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

		//** >> RESET **/
		GL_clock_begin = clock();
		if (reset() == _FAILURE_)
		{
			printf("\n\n Error running initialization() function\n\n");
			return _FAILURE_;
		}
		GL_times[11] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

		printf("tree adaptation\n");
		//** >> TREE ADAPTATION **/
		GL_clock_begin = clock();
		if (tree_adaptation() == _FAILURE_)
		{
			printf("\n\n Error running tree_construction() function\n\n");
			return _FAILURE_;
		}
		GL_times[10] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

		printf("density \n");
		//** >> DENSITY COMPUTATION **/
		GL_clock_begin = clock();
		grid_density();
		GL_times[4] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

		// printf("Potential\n");
		// //** >> POTENTIAL COMPUTATION **/
		// GL_clock_begin = clock();
		// if (potential() == _FAILURE_)
		// {
		// 	printf("\n\n Error running potential() function\n\n");
		// 	return _FAILURE_;
		// }
		// GL_times[5] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

		printf("Grid acceleration\n");
		//** >> GRID ACCELERATION **/
		GL_clock_begin = clock();
		if (grid_acceleration() == _FAILURE_)
		{
			printf("\n\n Error running grid_acceleration() function\n\n");
			return _FAILURE_;
		}
		GL_times[6] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

		//** >> PARTICLE ACCELERATION **/
		GL_clock_begin = clock();
		if (particle_acceleration() == _FAILURE_)
		{
			printf("\n\n Error running particle_acceleration() function\n\n");
			return _FAILURE_;
		}
		GL_times[7] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;


		//** >> PARTICLE UPDATING B **/
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
			GL_clock_begin = clock();
			observables(energies);
			GL_times[13] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

			//** >> OUPUT SNAPSHOTS **/
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

		Number_timesteps++; // Increasing the number of time-steps
	}

	printf("\n\n%sMaxdt = %f, Final time = %f%s\n\n",KRED ,Maxdt,actualtime,KNRM);

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
