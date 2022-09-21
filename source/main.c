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

int main() {

	//printf(std::scientific);
	//printf(std::setprecision(digits_precision));

	printf("\n ********** RUNNING N-BODY CODE ********** \n\n");

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

	printf("INITIAL PARAMETERS ...\n\n");
	printf("\nSIMULATION PARAMETERS ...\n");
	printf("User Box size = %f\n", (double)_User_BoxSize_);
	printf("Number of Particles = %d\n", (int)GL_no_ptcl_initial);
	printf("Max time of the simulation = %f My\n", (double)Maxdt / _Mgyear_);
	printf("Min level = %d\n", (int)lmin);
	printf("Max level = %d\n", (int)lmax);
	printf("Boundary type = %d\n", (int)boundary_type);
	printf("Output frecuency = %d\n", (int)fr_output);
	printf("\n\nUNITS ...\n");
	printf("Length unit = kpc\n");
	printf("Mass unit = 1 Solar Mass\n");
	printf("Time unit = %f\n", (double)tt);
	printf("G = %f\n", (double)_G_);
	printf("\n\nREFINEMENT CRITERIA ...\n");
	printf("Refinement criterion = %d particles\n", (int)ref_criterion_ptcl);
	printf("n_expand = %d\n", (int)n_exp);
	printf("CFL criterion = %f\n", (double)_CFL_);
	printf("Max timestep allow = %f My\n", (double)_MAX_dt_ / _Mgyear_);
	printf("\n\nPOISSON PARAMETERS ...\n");
	printf("Threshold over ther density = %f\n", (double)_ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_);
	printf("Threshold over ther previous solution = %f\n", (double)_ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_2);
	printf("Poisson error checking method = %d\n", (int)check_poisson_error_method);
	printf("Multigrid cycle = %d\n", (int)multigrid_cycle);
	printf("Solver type in the PreSmoothing of multigrid = %d\n", (int)solverPreS);
	printf("Solver type in the PostSmoothing of multigrid = %d\n", (int)solverPostS);
	printf("Solver type in the multigrid coarsest level = %d\n", (int)solverfinddphic);
	printf("Number of solver iterations in the PreSmoothing of multigrid = %d\n", (int)_NiterPreS_);
	printf("Number of solver iterations in the PostSmoothing of multigrid = %d\n", (int)_NiterPostS_);
	printf("Number of solver iterations in the multigrid coarsest level = %d\n", (int)_Niterfinddphic_);
	printf("Number of solver iterations in refinement levels = %d\n", (int)_Iter_branches_solver_);
	printf("Parameter of Successive Over Relaxation in the coarsest level = %f\n", (double)_w_SOR_HEAD_);
	printf("Parameter of Successive Over Relaxation in finner level = %f\n", (double)_w_SOR_);
	printf("Solver type for the corasest level of refinement = %d\n", (int)head_pot_method);
	printf("\n\nOTHER PARAMETERS ...\n");
	printf("Force stencil type: 3-points (=0) or 5-points (=1) = %d\n", (int)force_stencil);
	printf("Method for the computation of the observables = %d\n", (int)potential_energy_type);
	printf("Number of iterations between each Garbage collector = %d\n", (int)Garbage_Collector_iter);
	printf("\n\nBEGINNING the SIMULATION ...\n\n");

	//** >> INPUT **/
	//printf("Input\n");
	GL_clock_begin = clock();
	if (input() == _FAILURE_){
		printf("\n\n Error running input() function\n\n");
		return _FAILURE_;
	}
    GL_times[1] += (double) (clock() - GL_clock_begin) / CLOCKS_PER_SEC;

	//** >> INITIALIZATION **/
	// printf("Initialization\n");
	GL_clock_begin = clock();
	if (initialization() == _FAILURE_)
	{
		printf("\n\n Error running initialization() function\n\n");
		return _FAILURE_;
	}
	GL_times[2] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

	//** >> BUILDING INITIAL TREE **/
	// printf("Building initial Tree\n");
	GL_clock_begin = clock();
	if (tree_construction() == _FAILURE_)
	{
		printf("\n\n Error running tree_construction() function\n\n");
		return _FAILURE_;
	}
	GL_times[3] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

	//** >> DENSITY COMPUTATION **/
	// printf("Initial density\n");
	GL_clock_begin = clock();
	grid_density();
	GL_times[4] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

	//** >> POTENTIAL COMPUTATION **/
	// printf("Initial Potential\n");
	GL_clock_begin = clock();
	if (potential() == _FAILURE_)
	{
		printf("\n\n Error running potential() function\n\n");
		return _FAILURE_;
	}
	GL_times[5] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

	//** >> GRID ACCELERATION **/
	// printf("Initial Grid acceleration\n");
	GL_clock_begin = clock();
	if (grid_acceleration() == _FAILURE_)
	{
		printf("\n\n Error running grid_acceleration() function\n\n");
		return _FAILURE_;
	}
	GL_times[6] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

	//** >> PARTICLE ACCELERATION **/
	// printf("Initial Particle Acceleration\n");
	GL_clock_begin = clock();
	if (particle_acceleration() == _FAILURE_)
	{
		printf("\n\n Error running particle_acceleration() function\n\n");
		return _FAILURE_;
	}
	GL_times[7] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

	//** >> OBSERVABLES **/
	// printf("Initial Observables\n");
	GL_clock_begin = clock();
	observables(energies);
	GL_times[13] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

	//** >> OUPUT SNAPSHOTS **/
	// printf("Initial Snapshot output\n");
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
		//printf("\nTime-step compute\n");

		GL_clock_begin = clock();
		if (timestep_2(&dt) == _FAILURE_)
		{
			printf("\n\n Error running timestep_compute() function\n\n");
			return _FAILURE_;
		}
		GL_times[8] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

		actualtime += dt; // Updating actual time

		//** >> PARTICLE UPDATING A **/
		//printf("Particles Updating A\n");
		GL_clock_begin = clock();
		if (particle_updating_A(dt) == _FAILURE_)
		{
			printf("\n\n Error running particle_updating() function\n\n");
			return _FAILURE_;
		}
		GL_times[9] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

		//** >> RESET **/
		//printf("Reset\n");
		GL_clock_begin = clock();
		if (reset() == _FAILURE_)
		{
			printf("\n\n Error running initialization() function\n\n");
			return _FAILURE_;
		}
		GL_times[11] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

		//** >> TREE ADAPTATION **/
		//printf("tree adaptation\n");
		GL_clock_begin = clock();
		if (tree_adaptation() == _FAILURE_)
		{
			printf("\n\n Error running tree_adaptation() function\n\n");
			return _FAILURE_;
		}
		GL_times[10] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

		//** >> DENSITY COMPUTATION **/
		//printf("Density \n");
		GL_clock_begin = clock();
		grid_density();
		GL_times[4] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

		
		//** >> POTENTIAL COMPUTATION **/
		//printf("Potential\n");
		GL_clock_begin = clock();
		if (potential() == _FAILURE_)
		{
			printf("\n\n Error running potential() function\n\n");
			return _FAILURE_;
		}
		GL_times[5] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

		//** >> GRID ACCELERATION **/
		//printf("Grid acceleration\n");
		GL_clock_begin = clock();
		if (grid_acceleration() == _FAILURE_)
		{
			printf("\n\n Error running grid_acceleration() function\n\n");
			return _FAILURE_;
		}
		GL_times[6] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

		//** >> PARTICLE ACCELERATION **/
		//printf("Particle acceleration\n");
		GL_clock_begin = clock();
		if (particle_acceleration() == _FAILURE_)
		{
			printf("\n\n Error running particle_acceleration() function\n\n");
			return _FAILURE_;
		}
		GL_times[7] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

		//** >> PARTICLE UPDATING B **/
		//printf("Particle Updating B\n");
		GL_clock_begin = clock();
		particle_updating_B(dt);
		GL_times[12] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

		//printf("Finalized part updating B\n");

		//** >> OUPUT SNAPSHOTS **/
		//printf("Output and Observables:\n");
		if (Number_timesteps % fr_output == 0 || actualtime >= Maxdt)
		{
			++Number_outputs;

			//** >> OBSERVABLES **/
			//printf("Observables\n");
			GL_clock_begin = clock();
			if (observables(energies) == _FAILURE_)
			{
				printf("\n\n Error running observables() function\n\n");
				return _FAILURE_;
			}
			
			GL_times[13] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

			//** >> OUPUT SNAPSHOTS **/
			//printf("Output snapshots\n");
			GL_clock_begin = clock();
			if (output_snapshots(energies, actualtime, Number_outputs) == _FAILURE_)
			{
				printf("\n\n Error running output_snapshots() function\n\n");
				return _FAILURE_;
			}
			GL_times[19] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;
		}

		
		//** >> GARBAGE COLLECTOR **/
		//printf("Garbage Collector\n");
		GL_clock_begin = clock();
		if ((Number_timesteps + 1) % Garbage_Collector_iter == 0 && lmin < lmax)
		{
			//printf("Garbage Collector\n");
			if(garbage_collector() == _FAILURE_)
			{
				printf("\n\n Error running garbage_collector() function\n\n");
				return _FAILURE_;
			}
		}
		GL_times[14] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

		// BAR PROGRESS
		//printf("Bar\n");
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

		//printf("end cycle\n");
	}

	printf("\n\n%sMaxdt = %1.3f Myr, Final time = %1.3f Myr %s\n\n", KRED, (double)(Maxdt / _Mgyear_), (double)(actualtime / _Mgyear_), KNRM);

	/**
	*! ******************************************************************************
	*! ******************************************************************************
	*! ****************************************************************************** 
	*/

	//** >> OUPUT MAIN PARAMETERS **/
	// printf("Output main parameters\n");
	GL_clock_begin = clock();
	output_main_parameters(actualtime, Number_timesteps, Number_outputs);
	GL_times[18] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;

	//** >> TERMINAL PRINT **/
	terminal_print();

	printf("Total number of time-steps = %d\n", Number_timesteps);

	printf("\n\n ******** FINALIZED N-BODY CODE ********** \n\n");

	return _SUCCESS_;
}
