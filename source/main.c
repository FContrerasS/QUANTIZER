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

/**
 * @file main.c
 *
 * \f[{\color{magenta} \mathbf{ DOCUMENTED\ ``main.c"\ MODULE}}\f]
 *
 * @brief Main function in charge of executing the modules.
 *
 * **VERSION INFORMATION**: Felipe Contreras, 2022-10-01, version 1.0.
 *
 * - [a] **Node**: The most basic structure of the tree which represents an
 *                isolated refinement zone at a given level of refinement.
 *
 * - [b] **Tentacles**:  Pointer structure which has access to every node in
 *                      every level of refinement in the tree.
 *
 * **SHORT DESCRIPTION**: The code Begins and ends in this scrip. The main.c
 * script executes all modules; those who are executed only (a) \b one \b time:
 * global_variables.c, input.c, initialization.c, tree_construction.c, and those
 * who are executed (b) \b cyclically: grid_density.c, potential.c,
 * grid_acceleration.c, particle_acceleration.c observables.c,
 * output_snapshots.c, timestep.c, particle_updating.c, and tree_adaptation.c
 *
 * **PREREQUISITES**: Always used.
 *
 * **RETURN**: No parameter is returned.
 *
 * **LONG DESCRIPTION**:
 *
 * The goal of the main.c script is to call every module in the code to perform
 * the whole N-body simulation.
 *
 * The flux of this function can be seen in the figure (work in progress), and
 * it is explained below:
 *
 * - [0]  <b> THE N-BODY SIMULATION STARTS....</b>
 *
 * - [1]  Reading global variables and the input parameters with the module
 *   global_variables.c.
 *
 * - [2]  Reading the user's initial particle distribution with the module
 *   input.c.
 *
 * - [3]  Using the module initialization.c, \b (a) some particles global
 *   parameters, \b (b) the coarsest level of refinement, also known as the
 *   header node, and \b (c) the tentacles which go through all the nodes in the
 *   tree are initialized.
 *
 * - [4]  If there is more than 1 level of refinement, the tree of nodes is
 *   built with the module tree_construction.c.
 *
 * - [5]  Using the initial particle distribution given by the user, the grid
 *   density in every node is computed through the module grid_density.c.
 *   Currently the only method to perform this is the Cloud In Cell (CIC)
 *   scheme.
 *
 * - [6]  Using the grid density as input in the Poisson equation, the module
 *   potential.c computes the grid potential for every node. Currently the
 *   methods implemented to perform this task are the method of Multigrid with
 *   Gauss-Saidel or Jacoby, and the Successive Over Relaxation method (SOR).
 *
 * - [7]  Using the density in every node the grid acceleration is computing
 *   through the module grid_acceleration.c. There are 2 methods currently
 *   implemented, 3-stencil and 5-stencil points. However, the boundary grid
 *   points only use the 3-stencil scheme.
 *
 * - [8]  In the module particle_acceleration.c, the acceleration of the
 *   particles is computed with the Cloud In Cell (CIC) scheme.
 *
 * - [9]  The observables are computed using the module observables.c. It is
 *   possible to compute the potential energy using the exact method of particle
 *   to particle interaction of order O(N^2) or with an approximation of order
 *   O(N) using the grid potential obtained through the Poisson equation.
 *
 * - [10] The particle properties and observables are exported using the module
 *   output_snapshots.c. At this point, the code has computed the 0 time-step.
 *
 * - [11] Now, a \c "while" cycle is performed running until reaching the final
 *   time or the maximum number of time-steps, both given by the user. The steps
 *   of the \c "while" cycle below:
 *
 * - [12] Using the Courant-Friedrichs-Levy number given by the user, the next
 *   time-step is computed through the module timestep.c. Currently there are
 *   two methods to perform this: using the velocity only or using the velocity
 *   and the acceleration.
 *
 * - [13] In the scheme of Leapfrog, the particles are updated by a half-step of
 *   of velocity (kick) and by a step of position (drift), also known as the \e
 *   "Predictor step", through the module particle_updating.c.
 *
 * - [14] The particle acceleration and the grid density of the coarsest level
 *   of refinement are reset in the module reset.c.
 *
 * - [15] The well known \e "Adaptive Mesh Refinement" (AMR) scheme is performed
 *   through the adaptation of the tree structure of nodes in the module
 *   tree_adaptation.c.
 *
 * - [16] Same as the step [5], the grid density is recomputed using the module
 *   grid_density.c
 *
 * - [17] Same as the step [6], the grid potential is recomputed using the
 *   module potential.c
 *
 * - [18] Same as the step [7], the grid acceleration is recomputed using the
 *   module grid_acceleration.c
 *
 * - [19] Same as the step [8], the the acceleration of the particles is
 *   recomputed using the module particle_acceleration.c
 *
 * - [20] Similar the step [13], the particle velocity is updated by the
 *   remaining half-step of the Leapfrog scheme (kick). This update is also
 *   known as the \e "Corrector step", and it is executed in the module
 *   particle_updating.c.
 *
 * - [21] Same as the step [9], the observables are recomputed using the module
 *   observables.c
 *
 * - [22] Same as the step [10], the particle properties and observables are
 *   recomputed using the module output_snapshots.c
 *
 * - [23] If the user deems it convenient, a Garbage Collector of RAM is
 *   performed through the function garbage_collector().
 *
 * - [24] At this point a time-step cycle has been performed. Steps [12] to step
 *   [23] are repeated until reaching the final time or the maximum number of
 *   time-steps, both given by the user.
 *
 * - [25] Finished the previous \c "while" cycle, a .txt file with the main
 *   parameters used in the simulation are exported through the function
 *   output_main_parameters().
 *
 * - [26] <b> THE N-BODY SIMULATION ENDS....</b>
 *
 * **RATIONALES**:
 *  - [a] Trivial.
 *
 * **NOTES**:
 *  - [a] ll
 */

#include "main.h"

/**
 * @brief Make the calls of all the modules for each N-body simulation.
 *
 * **SHORT DESCRIPTION**: Make the calls of all the modules for each simulation.
 *
 * <b> SEE main.c FOR A WHOLE DESCRIPTION.</b>
 */



int main()
{

  // initialize parallelization with a default of 8 threads
  

  // printf(std::scientific);
  // printf(std::setprecision(digits_precision));

  printf("\n ********** RUNNING QUANTIZER ********** \n\n");

  //* DEFINING OTHER VARIABLES
  //* >> Defining timestep **/
  vtype dt = 0;
  vtype actualtime = 0;
  int Number_timesteps = 0;
  int Number_outputs = 0;
  /*
   *! ******************************************************************************
   *! ******************************************************************************
   *! ******************************************************************************
   */

  

  // OBSERVABLES:
  //   GL_energies[0] = Kinetic
  //   GL_energies[1] = Potential
  //   GL_energies[2] = Total = K + P

  //* >> GLOBAL VARIABLES **/
  // printf("Global variables\n");
  //GL_clock_begin = clock();
  clock_gettime( CLOCK_REALTIME, &GL_start);
  global_variables();
  //GL_times[0] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;
  clock_gettime( CLOCK_REALTIME, &GL_finish);
  GL_times[0] += ( GL_finish.tv_sec - GL_start.tv_sec ) + ( GL_finish.tv_nsec - GL_start.tv_nsec )/ 1000000000.;


  //omp_set_num_threads(NUMBER_OF_THEADS);

  printf("INITIAL PARAMETERS ...\n\n");
  printf("\nSIMULATION PARAMETERS ...\n");
  printf("Number of Threads = %d\n", (int)NUMBER_OF_THEADS);
  printf("User Box size = %f\n", (double)_User_BoxSize_);
  printf("Number of Particles = %d\n", (int)GL_no_ptcl_initial);
  printf("Max time of the simulation = %e years\n", (double) (Maxdt / _conversion_time_));
  printf("Min level = %d\n", (int)lmin);
  printf("Max level = %d\n", (int)lmax);
  printf("Boundary type = %d\n", (int)bdry_cond_type);
  printf("Output frecuency = %d\n", (int)fr_output);
  printf("Max iterations = %d\n",(int)MaxIterations);
  printf("\n\nUNITS ...\n");
  printf("Length unit = kpc\n");
  printf("Mass unit = 1 Solar Mass\n");
  //printf("Time unit = %f\n", (double)tt);
  printf("G = %f\n", (double)_G_);
  printf("\n\nREFINEMENT CRITERIA ...\n");
  printf("Refinement criterion = %d particles\n", (int)ref_criterion_ptcl);
  printf("n_expand = %d\n", (int)n_exp);
  printf("CFL criterion = %f\n", (double)_CFL_);
  printf("Max timestep allow = %e My\n", (double) (_MAX_dt_ /_conversion_time_));
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

  //* >> INPUT **/
  // printf("Input\n");
  //GL_clock_begin = clock();
  clock_gettime( CLOCK_REALTIME, &GL_start);
  if (input() == _FAILURE_)
  {
    printf("\n\n Error running input() function\n\n");
    return _FAILURE_;
  }
  //GL_times[1] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;
  clock_gettime( CLOCK_REALTIME, &GL_finish);
  GL_times[1] += ( GL_finish.tv_sec - GL_start.tv_sec ) + ( GL_finish.tv_nsec - GL_start.tv_nsec )/ 1000000000.;

  //* >> INITIALIZATION **/
  // printf("Initialization\n");
  //GL_clock_begin = clock();
  clock_gettime( CLOCK_REALTIME, &GL_start);
  if (initialization() == _FAILURE_)
  {
    printf("\n\n Error running initialization() function\n\n");
    return _FAILURE_;
  }
  //GL_times[2] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;
  clock_gettime( CLOCK_REALTIME, &GL_finish);
  GL_times[2] += ( GL_finish.tv_sec - GL_start.tv_sec ) + ( GL_finish.tv_nsec - GL_start.tv_nsec )/ 1000000000.;

  //* >> BUILDING INITIAL TREE **/
  // printf("Building initial Tree\n");
  //GL_clock_begin = clock();
  clock_gettime( CLOCK_REALTIME, &GL_start);
  if (tree_construction() == _FAILURE_)
  {
    printf("\n\n Error running tree_construction() function\n\n");
    return _FAILURE_;
  }
  //GL_times[3] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;
  clock_gettime( CLOCK_REALTIME, &GL_finish);
  GL_times[3] += ( GL_finish.tv_sec - GL_start.tv_sec ) + ( GL_finish.tv_nsec - GL_start.tv_nsec )/ 1000000000.;

  //* >> DENSITY COMPUTATION **/
  //printf("Initial density\n");
  //GL_clock_begin = clock();
  clock_gettime( CLOCK_REALTIME, &GL_start);
  grid_density();
  //GL_times[4] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;
  clock_gettime( CLOCK_REALTIME, &GL_finish);
  GL_times[4] += ( GL_finish.tv_sec - GL_start.tv_sec ) + ( GL_finish.tv_nsec - GL_start.tv_nsec )/ 1000000000.;

  //* >> POTENTIAL COMPUTATION **/
  //printf("Initial Potential\n");
  //GL_clock_begin = clock();
  clock_gettime( CLOCK_REALTIME, &GL_start);
  if (potential() == _FAILURE_)
  {
    printf("\n\n Error running potential() function\n\n");
    return _FAILURE_;
  }
  //GL_times[5] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;
  clock_gettime( CLOCK_REALTIME, &GL_finish);
  GL_times[5] += ( GL_finish.tv_sec - GL_start.tv_sec ) + ( GL_finish.tv_nsec - GL_start.tv_nsec )/ 1000000000.;

  //* >> GRID ACCELERATION **/
  //printf("Initial Grid acceleration\n");
  //GL_clock_begin = clock();
  clock_gettime( CLOCK_REALTIME, &GL_start);
  if (grid_acceleration() == _FAILURE_)
  {
    printf("\n\n Error running grid_acceleration() function\n\n");
    return _FAILURE_;
  }
  //GL_times[6] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;
  clock_gettime( CLOCK_REALTIME, &GL_finish);
  GL_times[6] += ( GL_finish.tv_sec - GL_start.tv_sec ) + ( GL_finish.tv_nsec - GL_start.tv_nsec )/ 1000000000.;

  //* >> PARTICLE ACCELERATION **/
  //printf("Initial Particle Acceleration\n");
  //GL_clock_begin = clock();
  clock_gettime( CLOCK_REALTIME, &GL_start);
  if (particle_acceleration() == _FAILURE_)
  {
    printf("\n\n Error running particle_acceleration() function\n\n");
    return _FAILURE_;
  }
  //GL_times[7] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;
  clock_gettime( CLOCK_REALTIME, &GL_finish);
  GL_times[7] += ( GL_finish.tv_sec - GL_start.tv_sec ) + ( GL_finish.tv_nsec - GL_start.tv_nsec )/ 1000000000.;

  //* >> OBSERVABLES **/
  //printf("Initial Observables\n");
  //GL_clock_begin = clock();
  clock_gettime( CLOCK_REALTIME, &GL_start);
  if (compute_observables_FLAG)
  {
    if (observables() == _FAILURE_)
    {
      printf("\n\n Error running observables() function\n\n");
      return _FAILURE_;
    }
  }
  //GL_times[13] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;
  clock_gettime( CLOCK_REALTIME, &GL_finish);
  GL_times[13] += ( GL_finish.tv_sec - GL_start.tv_sec ) + ( GL_finish.tv_nsec - GL_start.tv_nsec )/ 1000000000.;

  //* >> OUPUT SNAPSHOTS **/
  //printf("Initial Snapshot output\n");
  //GL_clock_begin = clock();
  clock_gettime( CLOCK_REALTIME, &GL_start);
  if (output_snapshots(actualtime, Number_outputs) == _FAILURE_)
  {
    printf("\n\n Error running output_snapshots() function\n\n");
    return _FAILURE_;
  }
  //GL_times[19] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;
  clock_gettime( CLOCK_REALTIME, &GL_finish);
  GL_times[19] += ( GL_finish.tv_sec - GL_start.tv_sec ) + ( GL_finish.tv_nsec - GL_start.tv_nsec )/ 1000000000.;

  while (actualtime < Maxdt && Number_timesteps < MaxIterations)
  {
    //* >> TIMESTEP **/
    //printf("\nTime-step compute\n");
    //GL_clock_begin = clock();
    clock_gettime( CLOCK_REALTIME, &GL_start);
    if (timestep(&dt) == _FAILURE_)
    {
      printf("\n\n Error running timestep_compute() function\n\n");
      return _FAILURE_;
    }

    //GL_times[8] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;
    clock_gettime( CLOCK_REALTIME, &GL_finish);
    GL_times[8] += ( GL_finish.tv_sec - GL_start.tv_sec ) + ( GL_finish.tv_nsec - GL_start.tv_nsec )/ 1000000000.;

    actualtime += dt; // Updating actual time

    //* >> PARTICLE UPDATING A **/
    // printf("Particles Updating A\n");
    //GL_clock_begin = clock();
    clock_gettime( CLOCK_REALTIME, &GL_start);
    if (particle_updating_A(dt) == _FAILURE_)
    {
      printf("\n\n Error running particle_updating() function\n\n");
      return _FAILURE_;
    }
    //GL_times[9] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;
    clock_gettime( CLOCK_REALTIME, &GL_finish);
    GL_times[9] += ( GL_finish.tv_sec - GL_start.tv_sec ) + ( GL_finish.tv_nsec - GL_start.tv_nsec )/ 1000000000.;

    //* >> RESET **/
    // printf("Reset\n");
    //GL_clock_begin = clock();
    clock_gettime( CLOCK_REALTIME, &GL_start);
    if (reset() == _FAILURE_)
    {
      printf("\n\n Error running initialization() function\n\n");
      return _FAILURE_;
    }
    //GL_times[11] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;
    clock_gettime( CLOCK_REALTIME, &GL_finish);
    GL_times[11] += ( GL_finish.tv_sec - GL_start.tv_sec ) + ( GL_finish.tv_nsec - GL_start.tv_nsec )/ 1000000000.;

    //* >> TREE ADAPTATION **/
    // printf("tree adaptation\n");
    //GL_clock_begin = clock();
    clock_gettime( CLOCK_REALTIME, &GL_start);
    if (tree_adaptation() == _FAILURE_)
    {
      printf("\n\n Error running tree_adaptation() function\n\n");
      return _FAILURE_;
    }
    //GL_times[10] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;
    clock_gettime( CLOCK_REALTIME, &GL_finish);
    GL_times[10] += ( GL_finish.tv_sec - GL_start.tv_sec ) + ( GL_finish.tv_nsec - GL_start.tv_nsec )/ 1000000000.;

    //* >> DENSITY COMPUTATION **/
    // printf("Density \n");
    //GL_clock_begin = clock();
    clock_gettime( CLOCK_REALTIME, &GL_start);
    grid_density();
    //GL_times[4] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;
    clock_gettime( CLOCK_REALTIME, &GL_finish);
    GL_times[4] += ( GL_finish.tv_sec - GL_start.tv_sec ) + ( GL_finish.tv_nsec - GL_start.tv_nsec )/ 1000000000.;

    //* >> POTENTIAL COMPUTATION **/
    // printf("Potential\n");
    //GL_clock_begin = clock();
    clock_gettime( CLOCK_REALTIME, &GL_start);
    if (potential() == _FAILURE_)
    {
      printf("\n\n Error running potential() function\n\n");
      return _FAILURE_;
    }
    //GL_times[5] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;
    clock_gettime( CLOCK_REALTIME, &GL_finish);
    GL_times[5] += ( GL_finish.tv_sec - GL_start.tv_sec ) + ( GL_finish.tv_nsec - GL_start.tv_nsec )/ 1000000000.;

    //* >> GRID ACCELERATION **/
    // printf("Grid acceleration\n");
    //GL_clock_begin = clock();
    clock_gettime( CLOCK_REALTIME, &GL_start);
    if (grid_acceleration() == _FAILURE_)
    {
      printf("\n\n Error running grid_acceleration() function\n\n");
      return _FAILURE_;
    }
    //GL_times[6] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;
    clock_gettime( CLOCK_REALTIME, &GL_finish);
    GL_times[6] += ( GL_finish.tv_sec - GL_start.tv_sec ) + ( GL_finish.tv_nsec - GL_start.tv_nsec )/ 1000000000.;

    //* >> PARTICLE ACCELERATION **/
    // printf("Particle acceleration\n");
    //GL_clock_begin = clock();
    clock_gettime( CLOCK_REALTIME, &GL_start);
    if (particle_acceleration() == _FAILURE_)
    {
      printf("\n\n Error running particle_acceleration() function\n\n");
      return _FAILURE_;
    }
    //GL_times[7] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;
    clock_gettime( CLOCK_REALTIME, &GL_finish);
    GL_times[7] += ( GL_finish.tv_sec - GL_start.tv_sec ) + ( GL_finish.tv_nsec - GL_start.tv_nsec )/ 1000000000.;

    //* >> PARTICLE UPDATING B **/
    // printf("Particle Updating B\n");
    //GL_clock_begin = clock();
    clock_gettime( CLOCK_REALTIME, &GL_start);
    particle_updating_B(dt);
    //GL_times[12] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;
    clock_gettime( CLOCK_REALTIME, &GL_finish);
    GL_times[12] += ( GL_finish.tv_sec - GL_start.tv_sec ) + ( GL_finish.tv_nsec - GL_start.tv_nsec )/ 1000000000.;

    // printf("Finalized part updating B\n");

    //* >> OUPUT SNAPSHOTS **/
    // printf("Output and Observables:\n");
    if (Number_timesteps % fr_output == 0 || actualtime >= Maxdt)
    {
      ++Number_outputs;

      //* >> OBSERVABLES **/
      // printf("Observables\n");
      //GL_clock_begin = clock();
      clock_gettime( CLOCK_REALTIME, &GL_start);
      if (compute_observables_FLAG)
      {
        if (observables() == _FAILURE_)
        {
          printf("\n\n Error running observables() function\n\n");
          return _FAILURE_;
        }
      }
      //GL_times[13] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;
      clock_gettime( CLOCK_REALTIME, &GL_finish);
      GL_times[13] += ( GL_finish.tv_sec - GL_start.tv_sec ) + ( GL_finish.tv_nsec - GL_start.tv_nsec )/ 1000000000.;
      
      //* >> OUPUT SNAPSHOTS **/
      // printf("Output snapshots\n");
      //GL_clock_begin = clock();
      clock_gettime( CLOCK_REALTIME, &GL_start);
      if (output_snapshots(actualtime, Number_outputs) == _FAILURE_)
      {
        printf("\n\n Error running output_snapshots() function\n\n");
        return _FAILURE_;
      }
      //GL_times[19] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;
      clock_gettime( CLOCK_REALTIME, &GL_finish);
      GL_times[19] += ( GL_finish.tv_sec - GL_start.tv_sec ) + ( GL_finish.tv_nsec - GL_start.tv_nsec )/ 1000000000.;
    }

    //* >> GARBAGE COLLECTOR **/
    // printf("Garbage Collector\n");
    //GL_clock_begin = clock();
    clock_gettime( CLOCK_REALTIME, &GL_start);
    if ((Number_timesteps + 1) % Garbage_Collector_iter == 0 && lmin < lmax)
    {
      // printf("Garbage Collector\n");
      if (garbage_collector() == _FAILURE_)
      {
        printf("\n\n Error running garbage_collector() function\n\n");
        return _FAILURE_;
      }
    }
    //GL_times[14] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;
    clock_gettime( CLOCK_REALTIME, &GL_finish);
    GL_times[14] += ( GL_finish.tv_sec - GL_start.tv_sec ) + ( GL_finish.tv_nsec - GL_start.tv_nsec )/ 1000000000.;

    // BAR PROGRESS
    // printf("Bar\n");
    fflush(stdout);
    printf("\r[");
    for (int j = 0; j < 50; j++)
    {
      if (j < (int)(50.0 * actualtime / Maxdt))
        printf("=");
      else if (j == (int)(50.0 * actualtime / Maxdt))
        printf(">");
      else
        printf(" ");
    }
    printf("] %.7f %%", (double)(actualtime * 100.0 / Maxdt));

    Number_timesteps++; // Increasing the number of time-steps


    printf("Max lv reached = %d\n", GL_tentacles_level_max + lmin);

    // printf("end cycle\n");
  }

  printf("\n\n%sMaxdt = %1.3f Myr, Final time = %1.3f Myr %s\n\n", KRED, (double)(Maxdt / _Myear_), (double)(actualtime / _Myear_), KNRM);

  /*
   *! ******************************************************************************
   *! ******************************************************************************
   *! ******************************************************************************
   */

  //* >> OUPUT MAIN PARAMETERS **/
  // printf("Output main parameters\n");
  //GL_clock_begin = clock();
  clock_gettime( CLOCK_REALTIME, &GL_start);
  output_main_parameters(actualtime, Number_timesteps, Number_outputs);
  //GL_times[18] += (double)(clock() - GL_clock_begin) / CLOCKS_PER_SEC;
  clock_gettime( CLOCK_REALTIME, &GL_finish);
  GL_times[18] += ( GL_finish.tv_sec - GL_start.tv_sec ) + ( GL_finish.tv_nsec - GL_start.tv_nsec )/ 1000000000.;

  //* >> TERMINAL PRINT **/
  terminal_print();

  printf("Total number of time-steps = %d\n", Number_timesteps);

  printf("\n\n ******** FINALIZED QUANTIZER ********** \n\n");


  


  return _SUCCESS_;
}
