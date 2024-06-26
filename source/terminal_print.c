/*
 * terminal_print.c
 *
 * Print variables in terminal
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

#include "terminal_print.h"

//* >> Local Functions
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

void terminal_print(void)
{
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

  printf("\n\n%sMEMORY [MB]:%s\n\n", KRED, KNRM);

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
    if (i == max_memory_secction)
    {
      printf("%s%s = %f ~ %.1f %% %s\n", KCYN, Memory_names[i], Total_memory[i] / 1000000.0, Total_memory[i] * 100.0 / sum, KNRM);
    }
    else
    {
      printf("%s = %f ~ %.1f %%\n", Memory_names[i], Total_memory[i] / 1000000.0, Total_memory[i] * 100.0 / sum);
    }
  }
  printf("Nodes in memory pool = %d\n", cntr_nodes_memory_pool);

  printf("\n%sTOTAL = %f MB%s\n\n", KMAG, sum / 1000000.0, KNRM);

  printf("\n\n%sTIEMPOS [s]:%s\n\n", KRED, KNRM);

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
      "Error Check in the Head",
      "Potential transfer in branches",
      "Filling Red and Black in branches",
      "Potential in branches",
      "Error Check in branches",
      ""};

  double TOTAL_TIME = 0.0;
  int MAX_TIME_SECCTION = 0;
  for (int i = 0; i < 20; i++)
  {
    TOTAL_TIME += GL_times[i];

    MAX_TIME_SECCTION = GL_times[MAX_TIME_SECCTION] > GL_times[i] ? MAX_TIME_SECCTION : i;
  }

  int MAX_TIME_SECCTION_POTENTIAL = 20;
  for (int i = 20; i < 26; i++)
  {
    MAX_TIME_SECCTION_POTENTIAL = GL_times[MAX_TIME_SECCTION_POTENTIAL] > GL_times[i] ? MAX_TIME_SECCTION_POTENTIAL : i;
  }

  //* >> Input time *//
  //* >> Initialization Time
  //* >> Tree Construction time *//
  //* >> Grid Density time *//
  for (int i = 0; i < 5; i++)
  {
    if (i == MAX_TIME_SECCTION)
    {
      printf("%s%s = %1.2e ~ %.1f %% %s\n", KCYN, Time_names[i], GL_times[i], GL_times[i] * 100.0 / TOTAL_TIME, KNRM);
    }
    else
    {
      printf("%s = %1.2e ~ %.1f %%\n", Time_names[i], GL_times[i], GL_times[i] * 100.0 / TOTAL_TIME);
    }
  }
  //* Potential times *//

  if (MAX_TIME_SECCTION == 5)
  {
    for (int i = 20; i < 26; i++)
    {
      if (MAX_TIME_SECCTION_POTENTIAL == i)
      {
        printf("%s%s = %1.2e ~ %.1f %% %s\n", KYEL, Time_names_potential[i - 20], GL_times[i], GL_times[i] * 100.0 / TOTAL_TIME, KNRM);
      }
      else
      {
        printf("%s = %1.2e ~ %.1f %%\n", Time_names_potential[i - 20], GL_times[i], GL_times[i] * 100.0 / TOTAL_TIME);
      }
    }
  }
  else
  {
    for (int i = 20; i < 26; i++)
    {
      printf("%s = %1.2e ~ %.1f %%\n", Time_names_potential[i - 20], GL_times[i], GL_times[i] * 100.0 / TOTAL_TIME);
    }
  }

  //* Total potential time *//
  //* >> Grid Acceleration time *//
  //* >> Particle acceleration time *//
  //* >> Time-step time *//
  //* >> Particle Updating time *//
  for (int i = 5; i < 15; i++)
  {
    if (i == MAX_TIME_SECCTION)
    {
      printf("%s%s = %1.2e ~ %.1f %% %s\n", KCYN, Time_names[i], GL_times[i], GL_times[i] * 100.0 / TOTAL_TIME, KNRM);
    }
    else
    {
      printf("%s = %1.2e ~ %.1f %%\n", Time_names[i], GL_times[i], GL_times[i] * 100.0 / TOTAL_TIME);
    }
  }

  //* >> Output Time *//
  if (18 == MAX_TIME_SECCTION)
  {
    printf("%sOutput Main Parameters  = %1.2e ~ %.1f %% %s\n", KCYN, GL_times[18], GL_times[18] * 100.0 / TOTAL_TIME, KNRM);
  }
  else
  {
    printf("Output Main Parameters = %.2e ~ %.1f %%\n", GL_times[18], GL_times[18] * 100.0 / TOTAL_TIME);
  }

  //* >> Output Time *//
  if (19 == MAX_TIME_SECCTION)
  {
    printf("%sOutput Snapshots = %.2e ~ %.1f %% %s\n", KCYN, GL_times[19], GL_times[19] * 100.0 / TOTAL_TIME, KNRM);
  }
  else
  {
    printf("Output Snapshots = %.2e ~ %.1f %%\n", GL_times[19], GL_times[19] * 100.0 / TOTAL_TIME);
  }

  //* >> Total time *//
  printf("\n%sTOTAL = %.3e s%s\n", KMAG, TOTAL_TIME, KNRM);

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
      "error",

      ""};

  printf("\n\n");
  //* >> TREE ADAPTATION TIME *//
  printf("\n\nTREE ADAPTATION TIME [s], Percentage over tree adaptation time\n\n");
  for (int i = 30; i < 51; i++)
  {
    printf("%d: %s = %1.2e ~ %.1f %%\n", i - 30, Time_names_tree_adaptation[i - 30], GL_times[i], GL_times[i] * 100.0 / GL_times[10]);
  }
  printf("\n\n");


  printf("total - observabables = %1.3f\n",TOTAL_TIME-GL_times[13]);


}
