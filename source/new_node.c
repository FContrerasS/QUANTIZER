/*
 * new_node.c
 *
 * Return a pointer to a node using the stack of the memory pool
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

#include "new_node.h"

//* >> Local Functions
static void re_initialize_node(struct node *ptr_node);

static void re_initialize_node(struct node *ptr_node)
{
  int box_idx_node;

  //* >> Boxes *//
  ptr_node->box_min_x = INT_MAX;  // Already minimal box value index in the real local space at the dimension X
  ptr_node->box_min_y = INT_MAX;  // Already minimal box value index in the real local space at the dimension Y
  ptr_node->box_min_z = INT_MAX;  // Already minimal box value index in the real local space at the dimension Z
  ptr_node->box_max_x = INT_MIN;  // Already maximum box value index in the real local space at the dimension X
  ptr_node->box_max_y = INT_MIN;  // Already maximum box value index in the real local space at the dimension Y
  ptr_node->box_max_z = INT_MIN;  // Already maximum box value index in the real local space at the dimension Z
  ptr_node->box_check_fit = true; // Check if the new box will fit in the old one

  //* >> Struct of cells (Particles and cell mass)
  for (int i = 0; i < ptr_node->cell_size; i++)
  {
    box_idx_node = ptr_node->ptr_box_idx[i];
    ptr_node->ptr_cell_struct[box_idx_node].cell_mass = 0.0;
    ptr_node->ptr_cell_struct[box_idx_node].ptcl_size = 0;
  }

  //* >> Cells in the node *//
  ptr_node->cell_size = 0; // Number of existing cells in the node

  //* >> Total mass in the node *//
  ptr_node->node_mass = 0.0;
  //* >> Total number of particles in the node
  ptr_node->no_ptcl_full_node = 0;
  ptr_node->no_ptcl_outs_ref_zones = 0;

  //* >> Grid points *//
  ptr_node->grid_intr_size = 0;                // Number of existing grid interior points in the block
  ptr_node->grid_bdry_size = 0;                // Number of existing grid border points in the block
  ptr_node->grid_sim_bdry_size = 0;     // Number of the boundary simulation grid points array of the block

  //* >> Potential, acceleration and density of the grid *//

  for (int i = 0; i < ptr_node->grid_properties_cap; i++)
  {
    ptr_node->ptr_d[i] = 0.0;
  }

  //* >> Tree structure *//
  ptr_node->chn_size = 0; // Number of children of the node

  //* >> Auxililary arrays to go from old box to new box *//
  // ptr_node->cell_ref_size = 0;	// Number of cells to be refined

  // Zones of refinement
  //  for (int j = 0; j < ptr_node->zones_cap; j++)
  //  {
  //    ptr_node->ptr_zone_size[j] = 0;
  //  }
  // ptr_node->zones_size = 0;		// Number of refined zones in the node

  // auxiliary booleans
  for (int j = 0; j < ptr_node->pbc_bool_bdry_anomalies_cap; j++)
  {
    ptr_node->ptr_pbc_bool_bdry_anomalies_x[j] = false;
    ptr_node->ptr_pbc_bool_bdry_anomalies_y[j] = false;
    ptr_node->ptr_pbc_bool_bdry_anomalies_z[j] = false;
  }

  // Subzones
  //  for (int j = 0; j < ptr_node->subzones_cap; j++)
  //  {
  //  	ptr_node->ptr_subzone_size[j] = 0;
  //  }
  // ptr_node->subzones_size = 0;	//NUmber of refined subzones in the node

  //* >> Boundary of the simulation box *//
  ptr_node->sim_bdry_contact = false;
  ptr_node->sim_bdry_contact_x = false;
  ptr_node->sim_bdry_contact_y = false;
  ptr_node->sim_bdry_contact_z = false;

  ptr_node->pbc_crosses_sim_box_bdry = false;
  ptr_node->pbc_crosses_sim_box_bdry_x = false; // when one node croses the simulation box at X direcction
  ptr_node->pbc_crosses_sim_box_bdry_y = false;
  ptr_node->pbc_crosses_sim_box_bdry_z = false;

  ptr_node->pbc_crosses_whole_sim_box = false;
  ptr_node->pbc_crosses_whole_sim_box_x = false;
  ptr_node->pbc_crosses_whole_sim_box_y = false;
  ptr_node->pbc_crosses_whole_sim_box_z = false;

  ptr_node->pbc_correction_due_pbc_flag = false;
}

struct node *new_node(void)
{
  struct node *ptr_node;

  if (GL_pool_node_start == NULL)
  {
    ptr_node = (struct node *)malloc(sizeof(struct node));
    initialize_node(ptr_node);
  }
  else
  {
    ptr_node = GL_pool_node_start;
    re_initialize_node(ptr_node);
    if (GL_pool_node_start == GL_pool_node_end)
    {
      GL_pool_node_start = NULL;
      GL_pool_node_end = NULL;
    }
    else
    {
      GL_pool_node_start = ptr_node->ptr_pt;
    }
  }

  return ptr_node;
}
