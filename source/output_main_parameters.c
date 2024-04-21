/*
 * output_main_parameters.c
 *
 * Export a .txt file with the main parameters used in the simulation
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

#include "output_main_parameters.h"

static int computing_memory(void);

static int computing_memory(void)
{

  int no_pts;

  struct node *ptr_node;

  //* >> Tentacles *//
  TOTAL_MEMORY_TENTACLES += (lmax - lmin + 1) * sizeof(struct node **);

  for (int lv = 0; lv < lmax - lmin + 1; lv++)
  {
    TOTAL_MEMORY_TENTACLES += GL_tentacles_cap[lv] * sizeof(struct node *);
  }

  //* >> Global particles *//
  TOTAL_MEMORY_PARTICLES += GL_no_ptcl_initial * (sizeof(bool) + sizeof(int) + 10.0 * sizeof(vtype)); // Global particles

  //* >> Node properties *//
  for (int lv = 0; lv < GL_tentacles_level_max + 1; lv++)
  {
    no_pts = GL_tentacles_size[lv];

    for (int i = 0; i < no_pts; i++)
    {
      ptr_node = GL_tentacles[lv][i];
      TOTAL_MEMORY_CELDAS += 4.0 * ptr_node->cell_cap * sizeof(int);

      if (lmin < lmax)
      {
        TOTAL_MEMORY_CELL_STRUCT += ptr_node->box_cap * sizeof(struct cell_struct);
        TOTAL_MEMORY_CELL_STRUCT += ptr_node->cell_struct_old_cap * sizeof(struct cell_struct);
        if (ptr_node->ptr_cell_struct != NULL)
        {
          for (int j = 0; j < ptr_node->box_cap; j++)
          {
            TOTAL_MEMORY_CELL_STRUCT += ptr_node->ptr_cell_struct[j].ptcl_cap * sizeof(int);
          }
        }

        if (ptr_node->ptr_cell_struct_old != NULL)
        {
          for (int j = 0; j < ptr_node->cell_struct_old_cap; j++)
          {
            TOTAL_MEMORY_CELL_STRUCT += ptr_node->ptr_cell_struct_old[j].ptcl_cap * sizeof(int);
          }
        }
      }
      TOTAL_MEMORY_CAJAS += 2.0 * ptr_node->box_cap * (sizeof(int) + sizeof(vtype)); // Boxes and mass boxes

      TOTAL_MEMORY_GRID_POINTS += 4.0 * (ptr_node->grid_bdry_cap + ptr_node->grid_intr_cap + ptr_node->grid_sim_bdry_cap) * sizeof(int); // Grid interior and border points

      TOTAL_MEMORY_GRID_PROPERTIES += 6.0 * ptr_node->grid_properties_cap * sizeof(vtype); // Grid properties, accelerations, potential and density
      TOTAL_MEMORY_AUX += ptr_node->zones_cap * sizeof(int *);
      for (int j = 0; j < ptr_node->zones_cap; j++)
      {
        TOTAL_MEMORY_AUX += ptr_node->ptr_zone_cap[j] * sizeof(int);
      }
      TOTAL_MEMORY_AUX += ptr_node->aux_idx_cap * sizeof(int);
      TOTAL_MEMORY_AUX += ptr_node->pbc_bool_bdry_anomalies_cap * 3.0 * sizeof(bool);

      TOTAL_MEMORY_AUX += 4.0 * ptr_node->links_cap * sizeof(int);

      //* Subzones
      TOTAL_MEMORY_AUX += 6.0 * ptr_node->pbc_subzones_cap * sizeof(int);
      // TOTAL_MEMORY_AUX += ptr_node->pbc_subzones_cap * sizeof(int *);
      // for(int j=0; j < ptr_node->pbc_subzones_cap; j++)
      // {
      //     TOTAL_MEMORY_AUX += ptr_node->ptr_subzone_cap[j] * sizeof(int);
      // }
    }
    TOTAL_MEMORY_NODES += no_pts * sizeof(struct node);
  }

  // Stack of memory pool
  int cntr_nodes_memory_pool = 0;
  ptr_node = GL_pool_node_start;
  while (ptr_node != NULL)
  {
    cntr_nodes_memory_pool++;
    TOTAL_MEMORY_STACK += 4.0 * ptr_node->cell_cap * sizeof(int);
    TOTAL_MEMORY_STACK += ptr_node->box_cap * sizeof(struct cell_struct);

    TOTAL_MEMORY_STACK += ptr_node->box_cap * sizeof(struct cell_struct);
    TOTAL_MEMORY_STACK += ptr_node->cell_struct_old_cap * sizeof(struct cell_struct);
    if (ptr_node->ptr_cell_struct != NULL)
    {
      for (int j = 0; j < ptr_node->box_cap; j++)
      {
        TOTAL_MEMORY_STACK += ptr_node->ptr_cell_struct[j].ptcl_cap * sizeof(int);
      }
    }

    if (ptr_node->ptr_cell_struct_old != NULL)
    {
      for (int j = 0; j < ptr_node->cell_struct_old_cap; j++)
      {
        TOTAL_MEMORY_STACK += ptr_node->ptr_cell_struct_old[j].ptcl_cap * sizeof(int);
      }
    }

    TOTAL_MEMORY_STACK += 2 * ptr_node->box_cap * (sizeof(int) + sizeof(vtype));                                                          // Boxes and mass boxes
    TOTAL_MEMORY_STACK += 4 * (ptr_node->grid_bdry_cap + ptr_node->grid_intr_cap + ptr_node->grid_sim_bdry_cap) * sizeof(int);     // Grid interior, border and simulation boundary grid points
    TOTAL_MEMORY_STACK += 6 * ptr_node->grid_properties_cap * sizeof(vtype);                                                              // Grid properties, accelerations, potential and density
    TOTAL_MEMORY_STACK += ptr_node->zones_cap * sizeof(int *);
    for (int j = 0; j < ptr_node->zones_cap; j++)
    {
      TOTAL_MEMORY_STACK += ptr_node->ptr_zone_cap[j] * sizeof(int);
    }
    TOTAL_MEMORY_STACK += ptr_node->aux_idx_cap * sizeof(int);
    TOTAL_MEMORY_STACK += ptr_node->pbc_bool_bdry_anomalies_cap * 3.0 * sizeof(bool);
    TOTAL_MEMORY_STACK += sizeof(struct node);

    TOTAL_MEMORY_STACK += 4.0 * ptr_node->links_cap * sizeof(int);

    //* Subzones
    TOTAL_MEMORY_STACK += 6.0 * ptr_node->pbc_subzones_cap * sizeof(int);
    // TOTAL_MEMORY_STACK += ptr_node->pbc_subzones_cap * sizeof(int *);
    // for(int j=0; j < ptr_node->pbc_subzones_cap; j++)
    // {
    //     TOTAL_MEMORY_STACK += ptr_node->ptr_subzone_cap[j] * sizeof(int);
    // }

    if (ptr_node == GL_pool_node_end)
    {
      ptr_node = NULL;
    }
    else
    {
      ptr_node = ptr_node->ptr_pt;
    }
  }

  return cntr_nodes_memory_pool;
}

void output_main_parameters(vtype final_time, int Number_timesteps, int Number_outputs)
{
  char parameters_name[1100];
  sprintf(parameters_name, "%s/Parameters.dat", folder_name);

  FILE *file = NULL;
  file = fopen(parameters_name, "w");
  if (file == NULL)
  {
    printf("Error opening file!");
    exit(EXIT_FAILURE);
  }

  fprintf(file, "lmin %d\n", (int)lmin);
  fprintf(file, "lmax %d\n", (int)lmax);
  fprintf(file, "_MAX_NUMBER_OF_ITERATIONS_IN_POISSON_EQUATION_ %d\n", (int)_MAX_NUMBER_OF_ITERATIONS_IN_POISSON_EQUATION_);
  fprintf(file, "_ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_ %1.5e\n", (double)_ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_);
  fprintf(file, "_ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_2 %1.5e\n", (double)_ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_2);
  fprintf(file, "check_poisson_error_method %d\n", (int)check_poisson_error_method);
  fprintf(file, "multigrid_cycle %d\n", (int)multigrid_cycle);
  fprintf(file, "solverPreS %d\n", (int)solverPreS);
  fprintf(file, "solverfinddphic %d\n", (int)solverfinddphic);
  fprintf(file, "solverPostS %d\n", (int)solverPostS);
  fprintf(file, "_NiterPreS_ %d\n", (int)_NiterPreS_);
  fprintf(file, "_NiterPostS_ %d\n", (int)_NiterPostS_);
  fprintf(file, "_Niterfinddphic_ %d\n", _Niterfinddphic_);
  fprintf(file, "_Iter_branches_solver_ %d\n", (int)_Iter_branches_solver_);
  fprintf(file, "_w_SOR_ %1.5e\n", (double)_w_SOR_);
  fprintf(file, "_w_SOR_HEAD_ %1.5e\n", (double)_w_SOR_HEAD_);
  fprintf(file, "head_pot_method %d\n", (int)head_pot_method);
  fprintf(file, "_CFL_ %1.5e\n", (double)_CFL_);
  fprintf(file, "_MAX_dt_ %1.5e\n", (double)_MAX_dt_);
  fprintf(file, "ref_criterion_ptcl %d\n", (int)ref_criterion_ptcl);
  fprintf(file, "n_exp %d\n", (int)n_exp);
  fprintf(file, "force_stencil %d\n", (int)force_stencil);
  fprintf(file, "potential_energy_type %d\n", (int)potential_energy_type);
  fprintf(file, "GC_iter %d\n", (int)Garbage_Collector_iter);

  fprintf(file, "Ngrid %d\n", (int)(no_lmin_cell + 1));
  fprintf(file, "Nparticles_initial %d\n", (int)GL_no_ptcl_initial);
  fprintf(file, "Nparticles_final %d\n", (int)GL_no_ptcl_final);
  fprintf(file, "Number_timesteps	%d\n", (int)Number_timesteps);
  fprintf(file, "Number_outputs %d\n", (int)Number_outputs);
  fprintf(file, "fr_output %d\n", (int)fr_output);
  fprintf(file, "MaxIterations %d\n", (int)MaxIterations);
  fprintf(file, "BoxSize %1.5e\n", (double)_User_BoxSize_);
  fprintf(file, "vtype %d\n", _VTYPE_);
  fprintf(file, "Maxdt %1.5e years\n", (double)(final_time / _conversion_time_)); // In Mega years
  
  // fprintf(file, "Length_unit kpc\n");
  // fprintf(file, "Mass_unit Solar_mass\n");
  //fprintf(file, "Time_unit %1.5e\n", (double)tt);
  fprintf(file, "G %1.5e\n", (double)_G_);
  fprintf(file, "boundary_type %d\n", (int)bdry_cond_type);
  fprintf(file, "GL_total_mass_initial %f\n", (double)GL_total_mass_initial);
  fprintf(file, "GL_total_mass_final %f\n", (double)GL_total_mass_final);

  //* >> TIMES
  char Time_names[50][100] = {
      "Global variables",
      "Input",
      "Initialization",
      "Tree Construction",
      "Grid Density",
      "Total potential",
      "Grid Acceleration",
      "Particle acceleration",
      "Time-step computing",
      "Particle Updating A",
      "Tree Adaptation",
      "Reset",
      "Particle Updating B",
      "Observables",
      "Garbage Collector",
      "Output Main Parameters",
      "Output Snapshots",
      ""};

  char Time_names_potential[50][100] = {
      "Head Potential",
      "Error Check of Poisson equation in the Head",
      "Potential transfer in branches",
      "Filling Red and Black in branches",
      "Potential in branches",
      "Error Check of Poisson equation in branches",
      ""};

  double TOTAL_TIME = 0;
  for (int i = 0; i < 20; i++)
  {
    TOTAL_TIME += GL_times[i];
  }

  fprintf(file, "\n\nTIMES in seconds\n\n");

  fprintf(file, "TOTAL %.3e\n", TOTAL_TIME);

  for (int i = 0; i < 15; i++)
  {
    fprintf(file, "%s = %1.2e ~ %.1f %%\n", Time_names[i], GL_times[i], GL_times[i] * 100.0 / TOTAL_TIME);
  }
  fprintf(file, "Output Snapshots = %.2e ~ %.1f %%\n", GL_times[19], GL_times[19] * 100.0 / TOTAL_TIME);

  fprintf(file, "\n\nPOTENTIAL TIME [s], Percentage over total potential time\n\n");
  for (int i = 20; i < 26; i++)
  {
    fprintf(file, "%s = %1.2e ~ %.1f %%\n", Time_names_potential[i - 20], GL_times[i], GL_times[i] * 100.0 / GL_times[5]);
  }

  //* >> TREE ADAPTATION TIME *//
  char Time_names_tree_adaptation[50][100] = {
      "updating_cell_struct",
      "check_error",
      "initialization_node_boxes",
      "fill_cell_ref",
      "fill_zones_ref",
      "create_links",
      "remov_cells_nolonger_require_refinement",
      "adapt_child_nodes",
      "create_new_child_nodes",
      "moving_old_child_to_new_child",
      "moving_new_zones_to_new_child",
      "adding_boundary_simulation_box_status_to_children_nodes",
      "reorganization_child_node",
      "reorganization_grandchild_node",
      "updating_ref_zones_children",
      "filling_child_grid_point_arrays",
      "computing_no_ptcl_outside_refinement_zones",
      "moved_unused_child_node_to_memory_pool",
      "tentacles_updating",
      "updating_tentacles_max_lv",
      ""};

  fprintf(file, "\n\nTREE ADAPTATION TIME [s], Percentage over tree adaptation time\n\n");
  for (int i = 30; i < 50; i++)
  {
    fprintf(file, "%d: %s = %1.2e ~ %.1f %%\n", i - 30, Time_names_tree_adaptation[i - 30], GL_times[i], GL_times[i] * 100.0 / GL_times[10]);
  }





  //* >> MEMORY *//
  int cntr_nodes_memory_pool;

  TOTAL_MEMORY_NODES = 0.0;
  TOTAL_MEMORY_CELDAS = 0.0;
  TOTAL_MEMORY_PARTICLES = 0.0;
  TOTAL_MEMORY_CELL_STRUCT = 0.0;
  TOTAL_MEMORY_CAJAS = 0.0;
  TOTAL_MEMORY_GRID_POINTS = 0.0;
  TOTAL_MEMORY_GRID_PROPERTIES = 0.0;
  TOTAL_MEMORY_AUX = 0.0;
  TOTAL_MEMORY_TENTACLES = 0.0;
  TOTAL_MEMORY_STACK = 0.0;

  cntr_nodes_memory_pool = computing_memory();

  double sum = TOTAL_MEMORY_NODES + TOTAL_MEMORY_CELDAS + TOTAL_MEMORY_PARTICLES + TOTAL_MEMORY_CELL_STRUCT + TOTAL_MEMORY_CAJAS + TOTAL_MEMORY_GRID_POINTS + TOTAL_MEMORY_GRID_PROPERTIES + TOTAL_MEMORY_AUX + TOTAL_MEMORY_TENTACLES + TOTAL_MEMORY_STACK;
  int max_memory_secction = 0;

  double Total_memory[10] = {TOTAL_MEMORY_NODES,
                             TOTAL_MEMORY_CELDAS,
                             TOTAL_MEMORY_PARTICLES,
                             TOTAL_MEMORY_CELL_STRUCT,
                             TOTAL_MEMORY_CAJAS,
                             TOTAL_MEMORY_GRID_POINTS,
                             TOTAL_MEMORY_GRID_PROPERTIES,
                             TOTAL_MEMORY_AUX,
                             TOTAL_MEMORY_TENTACLES,
                             TOTAL_MEMORY_STACK};

  for (int i = 0; i < 10; i++)
  {
    max_memory_secction = Total_memory[max_memory_secction] > Total_memory[i] ? max_memory_secction : i;
  }

  fprintf(file,"\n\nMEMORY [MB]:\n\n");

  char Memory_names[50][100] = {
      "Nodes",
      "Celdas",
      "Particulas",
      "Cell struct",
      "Cajas",
      "Grid Points",
      "Grid Properties",
      "Auxiliary",
      "Tentacles",
      "Memory Pool",
  };



  for (int i = 0; i < 10; i++)
  {
    fprintf(file,"%s = %f ~ %.1f %%\n", Memory_names[i], Total_memory[i] / 1000000.0, Total_memory[i] * 100.0 / sum);
  }
  fprintf(file,"Nodes in memory pool = %d\n", cntr_nodes_memory_pool);

  fprintf(file,"\nTOTAL = %f MB\n\n", sum / 1000000.0);












  fclose(file);
}