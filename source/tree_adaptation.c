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
 * @file tree_adaptation.c ******************** Documented \e "tree_adaptation.c" module ******************** \n
 * @brief Adaptation of the tree to the new particle configuration.
 *
 * **VERSION INFORMATION**: Felipe Contreras, 2022-10-01, version 1.0.
 *
 * **Reminder**:
 * - [a] **Node**: The most basic structure of the tree which represents an
 * isolated refinement zone at a given level of refinement.
 *
 * - \anchor tree_adaptation__REMINDER__Tentacles [b] **Tentacles**: Pointer
 * structure which has access to every node in every level of refinement in the
 * tree.
 *
 * **SHORT DESCRIPTION**: Adaptation of the tree of refinement zones to the new
 * particle configuration. The nodes are created, removed, or adapted according
 * to the new refinement zones. The structure of tentacles is modified to
 * pointer to the new configurations of nodes of the tree.
 *
 * **PREREQUISITES**: Always called, but only used if the user uses more than
 * one level of refinement, i.e. \f$ l_{min} < l_{max} \f$.
 *
 * **RETURN**: The error status.
 *
 * **LONG DESCRIPTION**:
 *
 * The goal of the tree_adaptation.c module is to perform the adaptation of the
 * tree of refinements, creating, removing, and adapting the nodes. In this
 * module, the tentacles structure is also adapted to the new node
 * configuration. The modifications go from the finest level of refinement to
 * the coarsest one, running over every node of those levels. How this process
 * should be understood is by considering the parent node which is the current
 * node of the loop, its child nodes, and its grandchild nodes. In every loop
 * over a node, the parent node performs modifications over the spatial geometry
 * of its child nodes without modifying its geometry or the grandchild geometry,
 * but only modifying his refinement zones and who will be the new parents of
 * its grandchild nodes.
 *
 * The flux of this function can be seen in the figure (work in progress), and
 * it is explained below:
 *
 * - [0]  <b> THE tree_adaptation.c MODULE STARTS....</b>
 *
 * - [1]  The tree_adaptation() begins to work.
 *
 * - [2]  Some useful parameters are defined.
 *
 * - [3]  Run a "for" loop over the refinement levels but avoiding the finest 
 * refinement level. The loop goes from a finer refinement level to a coarser 
 * one.
 *
 * - [4]  The number of tentacles in the next level of refinement and the number
 * of parent nodes in the actual level of refinement are reinitialized.
 *
 * - [5]  Run a "for" loop over all the nodes in the level of refinement
 * begins. From now on we will call these nodes parent nodes.
 *
 * - [6]  The cell structure of the parent node is updated using the child nodes
 * through the function updating_cell_struct(). Here the cell information of the
 *  parent node for every cell in a refinement zone is updated.
 *
 * - [7]  Using the function initialization_node_boxes(), the box associated
 * with the parent node is reinitialized, and the old box is stored in an
 * auxiliary box.
 *
 * - [8]  With the new configuration of the particles obtained in the
 * particle_updating.c module, cells in the parent node are labeled as
 * refinement cells (-1 in the code language) in the reinitialized box through
 * the function fill_cell_ref().
 *
 * - [9]  Using the new map of refinement cells, the different refinement zones
 * (blocks in code language) in the parent node are built through the function
 * fill_zones_ref(). Every one of these zones is independent of the other.
 *
 * - [10] Using the old map of refinement and the new map of refinement of the
 * parent node (old and new boxes), a link between the old zones and new zones
 * of refinement is created through the function create_links(). Every old zone
 * of refinement is associated with a single non-associated new zone of
 * refinement. In the following steps the old refinement zones associated with a
 * new one will be adapted, and for the new zones without association new child
 * nodes will be created to match with these zones. Remember that refinement
 * zones correspond to child nodes of the parent node.
 *
 * - [11] Cells that previously required refinement but are no anymore removed
 * from the child node associated with that refinement zone of the parent node
 * using the function remov_cells_nolonger_require_refinement().
 *
 * - [12] Using the function adapt_child_nodes(), the already child nodes
 * corresponding to the old refinement zones linked with a new refinement zone
 * are adapted in their parameter to fit with the new zone of refinement.
 *
 * - [13] Using the function create_new_child_nodes(), the new refinement zones
 * without a link with any old zone are fitted in a new child node obtained
 * using the memory pool.
 *
 * - [14] Every new refinement zone is associated with an old child node or with
 * a new child node, but the cell structure of these child nodes needs to be
 * inserted. Using the function transferring_old_child_to_new_child(), the old 
 * child cells which belong to old refinement zones are transferred to their new 
 * child node which corresponds to the linked new zone of refinement. Another 
 * way to understand the goal of this function is that it moves cells from one 
 * child node to another child node.
 *
 * - [15] Some cells require refinement in the new map but didn't on the old
 * map. These cells are put in the corresponding new child nodes through the
 * function transferring_new_zones_to_new_child().
 *
 * - [16] As the new child nodes can have different IDs of identification, they
 * require to be sorted properly. The function reorganization_child_node() do
 * that task.
 *
 * - [17] Similar to the step before, the grandchild nodes need to be sorted
 * properly. The function reorganization_grandchild_node() connect every
 * grandchild node with its right child node.
 *
 * - [18] Thanks to the previous step, the grandchild nodes are put in the right
 * child node. However, they are not necessarily sorted according to the child
 * box status of the refinement zone. So, the function
 * updating_ref_zones_children() is in charge of modifying the child box
 * status to add the right grandchild node identification.
 *
 * - [19] Besides to contain the status of every cell in the node as No-Exist
 * (-4), exist (-3), refinement (\f$ \geq \f$ 0), the boxes of the nodes
 * contains the information if the node is in the boundary of the simulation
 * box. To add this information, the function
 * adding_boundary_simulation_box_status_to_children_nodes() is called.
 *
 * - [20] Using the function filling_child_grid_point_arrays() the grid points 
 * of the child nodes are updated.
 *
 * - [21] In a "for" loop, using the function
 * computing_no_ptcl_outside_refinement_zones() the number of particles
 * outside of the refinement zones in every child node is computed.
 *
 * - [22] Using the function transferring_unused_child_node_to_memory_pool() the 
 * child nodes which was removed from the parent node are put in the memory pool 
 * to be used in the future for any other parent node as its child node.
 *
 * - [23] Using the function tentacles_updating() the tentacles structure is
 * updated with the new chlild node distribution.
 *
 * - [24] At this point the loop over the level of refinement and the loop over
 * the parent nodes end. The loops went from step [3] to step [23].
 *
 * - [25] The particles in the Head node, which represent the whole coarsest
 * level of refinement, are computed as in the step [21] using the function
 * computing_no_ptcl_outside_refinement_zones().
 *
 * - [26] The new maximum number of level of refinement is updated through the
 * function updating_tentacles_max_lv().
 *
 * - [27] <b> THE tree_adaptation.c MODULE ENDS....</b>
 *
 * **ILLUSTRATIVE EXAMPLES**:
 * - [a]  Trivial.
 *
 * **RATIONALES**:
 * - [a]  It is possible to join several functions in only one function, for
 * example, the functions fill_cell_ref() and fill_zones_ref(), which are in
 * charge of labeling the new box of the parent node with the cells to be
 * refined, and create the refinement zones using this new map of refinement
 * respectively. However, we decide to separate the functions as much as
 * possible to increase the modularity and the cleanness of the code.
 *
 * **NOTES**:
 * - [a]  The order in which some functions are called by the code is important.
 * For example, the function which removes cells from old child nodes
 * remov_cells_nolonger_require_refinement(), must be called before the function
 * which adapts the old child nodes to fit the new ones adapt_child_nodes().
 * This is because the adaptation of the old child nodes can rebuild the old box
 * to a new one erasing the old information.
 *
 */

#include "tree_adaptation.h"

//* >> Local Functions
void check_error(struct node *ptr_node, int type);
static int find_min_max_subzones_ref_PERIODIC_BOUNDARY(struct node *ptr_node, int zone_idx);
static int updating_cell_struct(struct node *ptr_node);
static void initialization_node_boxes(struct node *ptr_node);
static int fill_cell_ref(struct node *ptr_node);
static int fill_zones_ref(struct node *ptr_node);
static int create_links(struct node *ptr_node);
static void remov_cells_nolonger_require_refinement(struct node *ptr_node);
static int adapt_child_nodes(struct node *ptr_node);
static int create_new_child_nodes(struct node *ptr_node);
static int transferring_old_child_to_new_child(struct node *ptr_node);
static int transferring_new_zones_to_new_child(struct node *ptr_node);
static void reorganization_child_node(struct node *ptr_node);
static int reorganization_grandchild_node(struct node *ptr_node);
static void updating_ref_zones_children(struct node *ptr_node);
static void adding_boundary_simulation_box_status_to_children_nodes(struct node *ptr_node);
static int filling_child_grid_point_arrays(struct node *ptr_node);
static void computing_no_ptcl_outside_refinement_zones(struct node *ptr_node);
static void transferring_unused_child_node_to_memory_pool(struct node *ptr_node);
static int tentacles_updating(const struct node *ptr_node, int tentacle_lv);
static void updating_tentacles_max_lv(void);

void check_error(struct node *ptr_node, int type)
{
  int box_idx_x_ch;
  int box_idx_y_ch;
  int box_idx_z_ch;
  int box_idx_ch;

  int box_idx_x_node;
  int box_idx_y_node;
  int box_idx_z_node;
  int box_idx_node;

  struct node *ptr_ch;

  int box_real_dim_X_node = ptr_node->box_real_dim_x;
  int box_real_dim_X_times_Y_node = ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;

  int box_real_dim_X_ch;
  int box_real_dim_X_times_Y_ch;

  int ptcl_idx;

  int lv = ptr_node->lv;

  int aux_idx_x;
  int aux_idx_y;
  int aux_idx_z;

  // Checking ID of children
  for (int i = 0; i < ptr_node->chn_size; i++)
  {
    ptr_ch = ptr_node->pptr_chn[i];
    if (ptr_ch->ID != i)
    {
      printf("\n error, child ID = %d, its different to his position as a child = %d\n", ptr_ch->ID, i);
    }
  }

  // Check ptr_box_idx
  // Children:
  for (int i = 0; i < ptr_node->chn_size; i++)
  {
    ptr_ch = ptr_node->pptr_chn[i];

    box_real_dim_X_ch = ptr_ch->box_real_dim_x;
    box_real_dim_X_times_Y_ch = ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;

    for (int cell_idx = 0; cell_idx < ptr_ch->cell_size; cell_idx++)
    {

      aux_idx_x = ptr_ch->ptr_cell_idx_x[cell_idx];
      aux_idx_y = ptr_ch->ptr_cell_idx_y[cell_idx];
      aux_idx_z = ptr_ch->ptr_cell_idx_z[cell_idx];

      box_idx_x_ch = aux_idx_x - ptr_ch->box_ts_x;
      box_idx_y_ch = aux_idx_y - ptr_ch->box_ts_y;
      box_idx_z_ch = aux_idx_z - ptr_ch->box_ts_z;

      if (ptr_ch->pbc_crosses_sim_box_bdry == true)
      {
        if (aux_idx_x > ptr_ch->box_max_x)
        {
          box_idx_x_ch -= (1 << (lv + 1));
        }

        if (aux_idx_y > ptr_ch->box_max_y)
        {
          box_idx_y_ch -= (1 << (lv + 1));
        }

        if (aux_idx_z > ptr_ch->box_max_z)
        {
          box_idx_z_ch -= (1 << (lv + 1));
        }
      }

      box_idx_ch = box_idx_x_ch + box_idx_y_ch * box_real_dim_X_ch + box_idx_z_ch * box_real_dim_X_times_Y_ch;
      if (box_idx_ch != ptr_ch->ptr_box_idx[cell_idx])
      {

        if (type == 1)
        {
          printf("\npre check\n");
        }

        else
        {
          printf("\npost check\n");
        }

        printf("Child node at lv = %d, ID = %d\n", ptr_ch->lv, ptr_ch->ID);
        printf("\nError, child, ptr_box_idx no es igual \n");
        printf("box_idx_ch from cell idx = %d, ptr_ch box_idx = %d\n", box_idx_ch, ptr_ch->ptr_box_idx[cell_idx]);
      }
    }
  }

  // Parent node
  for (int cell_idx = 0; cell_idx < ptr_node->cell_size; cell_idx++)
  {

    aux_idx_x = ptr_node->ptr_cell_idx_x[cell_idx];
    aux_idx_y = ptr_node->ptr_cell_idx_y[cell_idx];
    aux_idx_z = ptr_node->ptr_cell_idx_z[cell_idx];

    box_idx_x_node = aux_idx_x - ptr_node->box_ts_x;
    box_idx_y_node = aux_idx_y - ptr_node->box_ts_y;
    box_idx_z_node = aux_idx_z - ptr_node->box_ts_z;

    if (ptr_node->pbc_crosses_sim_box_bdry == true)
    {
      if (aux_idx_x > ptr_node->box_max_x)
      {
        box_idx_x_node -= (1 << lv);
      }

      if (aux_idx_y > ptr_node->box_max_y)
      {
        box_idx_y_node -= (1 << lv);
      }

      if (aux_idx_z > ptr_node->box_max_z)
      {
        box_idx_z_node -= (1 << lv);
      }
    }

    box_idx_node = box_idx_x_node + box_idx_y_node * box_real_dim_X_node + box_idx_z_node * box_real_dim_X_times_Y_node;

    if (box_idx_node != ptr_node->ptr_box_idx[cell_idx])
    {
      if (type == 1)
      {
        printf("\npre check\n");
      }

      else
      {
        printf("\npost check\n");
      }

      printf("Parent node at lv = %d, ID = %d\n", ptr_node->lv, ptr_node->ID);
      printf("Error, parent ptr_box_idx no es igual \n");
    }
  }

  // Checking box value is in agreement with cell in the child node
  int cntr;
  int cntr2;

  for (int i = 0; i < ptr_node->chn_size; i++)
  {
    cntr = 0;
    ptr_ch = ptr_node->pptr_chn[i];

    for (int j = 0; j < ptr_node->box_real_dim_x * ptr_node->box_real_dim_y * ptr_node->box_real_dim_z; j++)
    {
      if (ptr_node->ptr_box[j] == i)
      {
        cntr++;
      }
    }
    if (ptr_ch->cell_size != 8 * cntr)
    {
      printf("Error, children size is different to parent box values * 8 (zone size)\n");
    }

    for (int j = 0; j < ptr_node->box_real_dim_x * ptr_node->box_real_dim_y * ptr_node->box_real_dim_z; j++)
    {
      if (ptr_node->ptr_box[j] == i)
      {
        cntr2 = 0;
        for (int k = 0; k < ptr_ch->cell_size; k++)
        {

          aux_idx_x = ptr_ch->ptr_cell_idx_x[k] >> 1;
          aux_idx_y = ptr_ch->ptr_cell_idx_y[k] >> 1;
          aux_idx_z = ptr_ch->ptr_cell_idx_z[k] >> 1;

          box_idx_x_node = aux_idx_x - ptr_node->box_ts_x;
          box_idx_y_node = aux_idx_y - ptr_node->box_ts_y;
          box_idx_z_node = aux_idx_z - ptr_node->box_ts_z;

          if (ptr_node->pbc_crosses_sim_box_bdry == true)
          {
            if (aux_idx_x > ptr_node->box_max_x)
            {
              box_idx_x_node -= (1 << lv);
            }

            if (aux_idx_y > ptr_node->box_max_y)
            {
              box_idx_y_node -= (1 << lv);
            }

            if (aux_idx_z > ptr_node->box_max_z)
            {
              box_idx_z_node -= (1 << lv);
            }
          }

          box_idx_node = box_idx_x_node + box_idx_y_node * box_real_dim_X_node + box_idx_z_node * box_real_dim_X_times_Y_node;

          if (box_idx_node == j)
          {
            cntr2++;
          }
        }
        if (cntr2 != 8)
        {
          printf("Child node at lv = %d, ID = %d\n", ptr_ch->lv, ptr_ch->ID);
          printf("Error, En el nodo hijo no se encuentra la celda del padre\n");
        }
      }
    }
  }

  // Repetitions of cells
  //  int aux_int;

  // //Children
  // for (int i = 0; i < ptr_node->chn_size; i++)
  // {
  //     ptr_ch = ptr_node->pptr_chn[i];
  //     for (int j = 0; j < ptr_ch->cell_size; j++)
  //     {
  //         aux_int = ptr_ch->ptr_box_idx[j];
  //         for (int k = j + 1; k < ptr_ch->cell_size; k++)
  //         {
  //             if (aux_int == ptr_ch->ptr_box_idx[k])
  //             {
  //                 printf("error, child repetition of cell\n");
  //             }
  //         }
  //     }
  // }

  // //Parent
  // for (int j = 0; j < ptr_node->cell_size; j++)
  // {
  //     aux_int = ptr_node->ptr_box_idx[j];
  //     for (int k = j + 1; k < ptr_node->cell_size; k++)
  //     {
  //         if (aux_int == ptr_node->ptr_box_idx[k])
  //         {
  //             printf("error, parent repetition of cell\n");
  //         }
  //     }
  // }

  // Checking box value with status No-exist with particles or mass > 0
  // Children
  for (int i = 0; i < ptr_node->chn_size; i++)
  {
    ptr_ch = ptr_node->pptr_chn[i];

    for (int j = 0; j < ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y * ptr_ch->box_real_dim_z; j++)
    {
      if (ptr_ch->ptr_box[j] < -3)
      {
        if (ptr_ch->ptr_cell_struct[j].cell_mass > 0 || ptr_ch->ptr_cell_struct[j].ptcl_size > 0)
        {
          printf("Child node at lv = %d, ID = %d\n", ptr_ch->lv, ptr_ch->ID);
          printf("error, child, has a box cell with status of No-Exist but with particles inside\n");
          printf("ptcl size = %d, cell mass = %f\n", ptr_ch->ptr_cell_struct[j].ptcl_size, ptr_ch->ptr_cell_struct[j].cell_mass);
          ptcl_idx = ptr_ch->ptr_cell_struct[j].ptr_ptcl[0];
          printf("ptcl id = %d, pos_x = %f, pos_y = %f, pos_z = %f\n", GL_ptcl_ID[ptcl_idx], (double)GL_ptcl_x[ptcl_idx], (double)GL_ptcl_y[ptcl_idx], (double)GL_ptcl_z[ptcl_idx]);
        }
      }
    }
  }

  // Checking cell mass vs cell particle
  //  Children
  for (int i = 0; i < ptr_node->chn_size; i++)
  {
    ptr_ch = ptr_node->pptr_chn[i];

    box_real_dim_X_ch = ptr_ch->box_real_dim_x;
    box_real_dim_X_times_Y_ch = ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;

    for (int cell_idx = 0; cell_idx < ptr_ch->cell_size; cell_idx++)
    {
      aux_idx_x = ptr_ch->ptr_cell_idx_x[cell_idx];
      aux_idx_y = ptr_ch->ptr_cell_idx_y[cell_idx];
      aux_idx_z = ptr_ch->ptr_cell_idx_z[cell_idx];

      box_idx_x_ch = aux_idx_x - ptr_ch->box_ts_x;
      box_idx_y_ch = aux_idx_y - ptr_ch->box_ts_y;
      box_idx_z_ch = aux_idx_z - ptr_ch->box_ts_z;

      if (ptr_ch->pbc_crosses_sim_box_bdry == true)
      {
        if (aux_idx_x > ptr_ch->box_max_x)
        {
          box_idx_x_ch -= (1 << (lv + 1));
        }

        if (aux_idx_y > ptr_ch->box_max_y)
        {
          box_idx_y_ch -= (1 << (lv + 1));
        }

        if (aux_idx_z > ptr_ch->box_max_z)
        {
          box_idx_z_ch -= (1 << (lv + 1));
        }
      }

      box_idx_ch = box_idx_x_ch + box_idx_y_ch * box_real_dim_X_ch + box_idx_z_ch * box_real_dim_X_times_Y_ch;

      if (ptr_ch->ptr_cell_struct[box_idx_ch].cell_mass != ptr_ch->ptr_cell_struct[box_idx_ch].ptcl_size * 100.0 || ptr_ch->ptr_cell_struct[box_idx_ch].ptcl_size < 0)
      {
        printf("Child node at lv = %d, ID = %d\n", ptr_ch->lv, ptr_ch->ID);
        printf("Error, child, box_idx = %d\n", box_idx_ch);
        printf("ptcl_size = %d, cell_mass = %f\n", ptr_ch->ptr_cell_struct[box_idx_ch].ptcl_size, (double)ptr_ch->ptr_cell_struct[box_idx_ch].cell_mass);
      }
    }

    if (ptr_ch->no_ptcl_full_node < ref_criterion_ptcl && type == 2)
    {
      printf("\nchild at lv = %d, ID = %d\n", ptr_ch->lv, ptr_ch->ID);
      printf("error, child no ptcl = %d is less than ref_criterion_ptcl = %d\n", ptr_ch->no_ptcl_full_node, ref_criterion_ptcl);
    }
  }

  // Parent node
  for (int cell_idx = 0; cell_idx < ptr_node->cell_size; cell_idx++)
  {
    aux_idx_x = ptr_node->ptr_cell_idx_x[cell_idx];
    aux_idx_y = ptr_node->ptr_cell_idx_y[cell_idx];
    aux_idx_z = ptr_node->ptr_cell_idx_z[cell_idx];

    box_idx_x_node = aux_idx_x - ptr_node->box_ts_x;
    box_idx_y_node = aux_idx_y - ptr_node->box_ts_y;
    box_idx_z_node = aux_idx_z - ptr_node->box_ts_z;

    if (ptr_node->pbc_crosses_sim_box_bdry == true)
    {
      if (aux_idx_x > ptr_node->box_max_x)
      {
        box_idx_x_node -= (1 << lv);
      }

      if (aux_idx_y > ptr_node->box_max_y)
      {
        box_idx_y_node -= (1 << lv);
      }

      if (aux_idx_z > ptr_node->box_max_z)
      {
        box_idx_z_node -= (1 << lv);
      }
    }

    box_idx_node = box_idx_x_node + box_idx_y_node * box_real_dim_X_node + box_idx_z_node * box_real_dim_X_times_Y_node;

    if (ptr_node->ptr_cell_struct[box_idx_node].cell_mass != ptr_node->ptr_cell_struct[box_idx_node].ptcl_size * 100.0 || ptr_node->ptr_cell_struct[box_idx_node].ptcl_size < 0)
    {
      printf("Parent node at lv = %d, ID = %d\n", ptr_node->lv, ptr_node->ID);
      printf("Error, parent, box_idx = %d\n", box_idx_node);
      printf("ptcl_size = %d, cell_mass = %f\n", ptr_node->ptr_cell_struct[box_idx_node].ptcl_size, (double)ptr_node->ptr_cell_struct[box_idx_node].cell_mass);
    }
  }
  // Notes: the parent node can have less particles than the ref_criterion becasue
  // when the parent node do not had children, then it is in fact a "child node"
  // if (ptr_node->no_ptcl_full_node < ref_criterion_ptcl && type == 2)
  // {
  //     printf("Parent at lv = %d, ID = %d\n", ptr_node->lv, ptr_node->ID);
  //     printf("error, Parent no ptcl = %d is less than ref_criterion_ptcl = %d\n", ptr_node->no_ptcl_full_node, ref_criterion_ptcl);
  // }

  // Checking maximum size of the box allow
  // Children
  for (int i = 0; i < ptr_node->chn_size; i++)
  {
    ptr_ch = ptr_node->pptr_chn[i];
    if (ptr_ch->box_max_x - ptr_ch->box_min_x >= (1 << (lv + 1)) ||
        ptr_ch->box_max_y - ptr_ch->box_min_y >= (1 << (lv + 1)) ||
        ptr_ch->box_max_z - ptr_ch->box_min_z >= (1 << (lv + 1)))
    {
      printf("Child node at lv = %d, ID = %d\n", ptr_ch->lv, ptr_ch->ID);
      printf("error, child, node has a dimension bigger than the simulation box\n");
    }
  }
  // Parent
  if (ptr_node->box_max_x - ptr_node->box_min_x >= (1 << lv) ||
      ptr_node->box_max_y - ptr_node->box_min_y >= (1 << lv) ||
      ptr_node->box_max_z - ptr_node->box_min_z >= (1 << lv))
  {
    printf("Parent node at lv = %d, ID = %d\n", ptr_node->lv, ptr_node->ID);
    printf("error, parent node has a dimension bigger than the simulation box\n");
  }

  // Checking maximum box dimensions:
  // Children
  for (int i = 0; i < ptr_node->chn_size; i++)
  {
    ptr_ch = ptr_node->pptr_chn[i];
    if (ptr_ch->box_max_x - ptr_ch->box_min_x + 1 != ptr_ch->box_dim_x ||
        ptr_ch->box_max_y - ptr_ch->box_min_y + 1 != ptr_ch->box_dim_y ||
        ptr_ch->box_max_z - ptr_ch->box_min_z + 1 != ptr_ch->box_dim_z)
    {
      printf("Child node at lv = %d, ID = %d\n", ptr_ch->lv, ptr_ch->ID);
      printf("error, child node has a box_dim different to max - min + 1 \n");
    }

    if (ptr_ch->box_max_x >= (1 << (lv + 1)) || ptr_ch->box_max_y >= (1 << (lv + 1)) || ptr_ch->box_max_z >= (1 << (lv + 1)))
    {
      printf("Child node at lv = %d, ID = %d\n", ptr_ch->lv, ptr_ch->ID);
      printf("error, child node box_max > max allow of %d \n", (1 << (lv + 1)));
    }

    if (ptr_ch->pbc_crosses_sim_box_bdry_x == false && ptr_ch->box_min_x < 0)
    {
      printf("Child node at lv = %d, ID = %d\n", ptr_ch->lv, ptr_ch->ID);
      printf("error, child node box_min_x < 0, but flag of crosees sim == false\n");
    }
    if (ptr_ch->pbc_crosses_sim_box_bdry_y == false && ptr_ch->box_min_y < 0)
    {
      printf("Child node at lv = %d, ID = %d\n", ptr_ch->lv, ptr_ch->ID);
      printf("error, child node box_min_y < 0, but flag of crosees sim == false\n");
    }
    if (ptr_ch->pbc_crosses_sim_box_bdry_z == false && ptr_ch->box_min_z < 0)
    {
      printf("Child node at lv = %d, ID = %d\n", ptr_ch->lv, ptr_ch->ID);
      printf("error, child node box_min_z < 0, but flag of crosees sim == false\n");
    }

    if (ptr_ch->sim_bdry_contact != (ptr_ch->sim_bdry_contact_x || ptr_ch->sim_bdry_contact_y || ptr_ch->sim_bdry_contact_z))
    {
      printf("Child node at lv = %d, ID = %d\n", ptr_ch->lv, ptr_ch->ID);
      printf("error, child node simulation contact flag dont match with simulation contact flag by direcction\n");
    }

    if (ptr_ch->pbc_crosses_sim_box_bdry != (ptr_ch->pbc_crosses_sim_box_bdry_x || ptr_ch->pbc_crosses_sim_box_bdry_y || ptr_ch->pbc_crosses_sim_box_bdry_z))
    {
      printf("Child node at lv = %d, ID = %d\n", ptr_ch->lv, ptr_ch->ID);
      printf("error, child node crosses the boundary simulation box flag dont match with simulation contact flag by direcction\n");
    }

    if (ptr_ch->pbc_crosses_whole_sim_box != (ptr_ch->pbc_crosses_whole_sim_box_x || ptr_ch->pbc_crosses_whole_sim_box_y || ptr_ch->pbc_crosses_whole_sim_box_z))
    {
      printf("Child node at lv = %d, ID = %d\n", ptr_ch->lv, ptr_ch->ID);
      printf("error, child node crosses the whole simulation box flag dont match with simulation contact flag by direcction\n");
    }

    if (ptr_ch->pbc_crosses_whole_sim_box == true)
    {
      if (ptr_ch->pbc_crosses_whole_sim_box_x == true)
      {
        if (ptr_ch->box_min_x != 0 && ptr_ch->box_max_x != (1 << (lv + 1)) - 1)
        {
          printf("Child node at lv = %d, ID = %d\n", ptr_ch->lv, ptr_ch->ID);
          printf("error, child node, the x direcction crosses the whole simulation box, but it has wrond box dim:\n");
          printf("box_dim_x = %d, box_min_x = %d, box_max_x = %d\n", ptr_ch->box_dim_x, ptr_ch->box_min_x, ptr_ch->box_max_x);
        }
      }

      if (ptr_ch->pbc_crosses_whole_sim_box_y == true)
      {
        if (ptr_ch->box_min_y != 0 && ptr_ch->box_max_y != (1 << (lv + 1)) - 1)
        {
          printf("Child node at lv = %d, ID = %d\n", ptr_ch->lv, ptr_ch->ID);
          printf("error, child node, the y direcction crosses the whole simulation box, but it has wrond box dim:\n");
          printf("box_dim_y = %d, box_min_y = %d, box_max_y = %d\n", ptr_ch->box_dim_y, ptr_ch->box_min_y, ptr_ch->box_max_y);
        }
      }

      if (ptr_ch->pbc_crosses_whole_sim_box_z == true)
      {
        if (ptr_ch->box_dim_z != 0 && ptr_ch->box_max_z != (1 << (lv + 1)) - 1)
        {
          printf("Child node at lv = %d, ID = %d\n", ptr_ch->lv, ptr_ch->ID);
          printf("error, child node, the z direcction crosses the whole simulation box, but it has wrond box dim:\n");
          printf("box_dim_z = %d, box_min_z = %d, box_max_z = %d\n", ptr_ch->box_dim_z, ptr_ch->box_min_z, ptr_ch->box_max_z);
        }
      }
    }
  }

  // Parent
  if (ptr_node->box_max_x - ptr_node->box_min_x + 1 != ptr_node->box_dim_x ||
      ptr_node->box_max_y - ptr_node->box_min_y + 1 != ptr_node->box_dim_y ||
      ptr_node->box_max_z - ptr_node->box_min_z + 1 != ptr_node->box_dim_z)
  {
    printf("Parent node at lv = %d, ID = %d\n", ptr_node->lv, ptr_node->ID);
    printf("error, Parent node has a box_dim different to max - min + 1 \n");
  }

  if (ptr_node->box_max_x >= (1 << lv) || ptr_node->box_max_y >= (1 << lv) || ptr_node->box_max_z >= (1 << lv))
  {
    printf("Parent node at lv = %d, ID = %d\n", ptr_node->lv, ptr_node->ID);
    printf("error, Parent node box_max > max allow of %d \n", 1 << lv);
  }

  if (ptr_node->pbc_crosses_sim_box_bdry_x == false && ptr_node->box_min_x < 0)
  {
    printf("Parent node at lv = %d, ID = %d\n", ptr_node->lv, ptr_node->ID);
    printf("error, Parent node box_min_x < 0, but flag of crosees sim == false\n");
  }
  if (ptr_node->pbc_crosses_sim_box_bdry_y == false && ptr_node->box_min_y < 0)
  {
    printf("Parent node at lv = %d, ID = %d\n", ptr_node->lv, ptr_node->ID);
    printf("error, Parent node box_min_y < 0, but flag of crosees sim == false\n");
  }
  if (ptr_node->pbc_crosses_sim_box_bdry_z == false && ptr_node->box_min_z < 0)
  {
    printf("Parent node at lv = %d, ID = %d\n", ptr_node->lv, ptr_node->ID);
    printf("error, Parent node box_min_z < 0, but flag of crosees sim == false\n");
  }

  if (ptr_node->sim_bdry_contact != (ptr_node->sim_bdry_contact_x || ptr_node->sim_bdry_contact_y || ptr_node->sim_bdry_contact_z))
  {
    printf("Parent node at lv = %d, ID = %d\n", ptr_node->lv, ptr_node->ID);
    printf("error, Parent node simulation contact flag dont match with simulation contact flag by direcction\n");
  }

  if (ptr_node->pbc_crosses_sim_box_bdry != (ptr_node->pbc_crosses_sim_box_bdry_x || ptr_node->pbc_crosses_sim_box_bdry_y || ptr_node->pbc_crosses_sim_box_bdry_z))
  {
    printf("Parent node at lv = %d, ID = %d\n", ptr_node->lv, ptr_node->ID);
    printf("error, Parent node crosses the boundary simulation box flag dont match with simulation contact flag by direcction\n");
  }

  if (ptr_node->pbc_crosses_whole_sim_box != (ptr_node->pbc_crosses_whole_sim_box_x || ptr_node->pbc_crosses_whole_sim_box_y || ptr_node->pbc_crosses_whole_sim_box_z))
  {
    printf("Parent node at lv = %d, ID = %d\n", ptr_node->lv, ptr_node->ID);
    printf("error, Parent node crosses the whole simulation box flag dont match with simulation contact flag by direcction\n");
  }

  if (ptr_node->pbc_crosses_whole_sim_box == true)
  {
    if (ptr_node->pbc_crosses_whole_sim_box_x == true)
    {
      if (ptr_node->box_min_x != 0 && ptr_node->box_max_x != (1 << lv) - 1)
      {
        printf("Parent node at lv = %d, ID = %d\n", ptr_node->lv, ptr_node->ID);
        printf("error, Parent node, the x direcction crosses the whole simulation box, but it has wrond box dim:\n");
        printf("box_dim_x = %d, box_min_x = %d, box_max_x = %d\n", ptr_node->box_dim_x, ptr_node->box_min_x, ptr_node->box_max_x);
      }
    }

    if (ptr_node->pbc_crosses_whole_sim_box_y == true)
    {
      if (ptr_node->box_max_y != 0 && ptr_node->box_max_y != (1 << lv) - 1)
      {
        printf("Parent node at lv = %d, ID = %d\n", ptr_node->lv, ptr_node->ID);
        printf("error, Parent node, the y direcction crosses the whole simulation box, but it has wrond box dim:\n");
        printf("box_dim_y = %d, box_min_y = %d, box_max_y = %d\n", ptr_node->box_dim_y, ptr_node->box_min_y, ptr_node->box_max_y);
      }
    }

    if (ptr_node->pbc_crosses_whole_sim_box_z == true)
    {
      if (ptr_node->box_max_z != 0 && ptr_node->box_max_z != (1 << lv) - 1)
      {
        printf("Parent node at lv = %d, ID = %d\n", ptr_node->lv, ptr_node->ID);
        printf("error, Parent node, the z direcction crosses the whole simulation box, but it has wrond box dim:\n");
        printf("box_dim_z = %d, box_min_z = %d, box_max_z = %d\n", ptr_node->box_dim_z, ptr_node->box_min_z, ptr_node->box_max_z);
      }
    }
  }

  // Particles are in the right cell
  //  Children
  for (int i = 0; i < ptr_node->chn_size; i++)
  {
    ptr_ch = ptr_node->pptr_chn[i];

    // box_real_dim_X_ch = ptr_ch->box_real_dim_x;
    // box_real_dim_X_times_Y_ch = ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;

    for (int cell_idx = 0; cell_idx < ptr_ch->cell_size; cell_idx++)
    {
      box_idx_ch = ptr_ch->ptr_box_idx[cell_idx];
      for (int j = 0; j < ptr_ch->ptr_cell_struct[box_idx_ch].ptcl_size; j++)
      {
        ptcl_idx = ptr_ch->ptr_cell_struct[box_idx_ch].ptr_ptcl[j];
        if (box_idx_ch != ptcl_idx_to_box_idx(ptr_ch, ptcl_idx))
        {
          if (type == 1)
          {
            printf("\npre check\n");
          }
          else
          {
            printf("\npost check\n");
          }
          printf("Error, child, particle is out of the cell\n");
          printf("parent: lv = %d, ID = %d\n", ptr_node->lv, ptr_node->ID);
          printf("child : lv = %d, ID = %d\n", ptr_ch->lv, ptr_ch->ID);
          printf("children size = %d\n", ptr_node->chn_size);
          printf("ptcl index = %d\n", ptcl_idx);
          printf("ptcl pos_x = %f\n", GL_ptcl_x[ptcl_idx]);
          printf("ptcl pos_y = %f\n", GL_ptcl_y[ptcl_idx]);
          printf("ptcl pos_z = %f\n", GL_ptcl_z[ptcl_idx]);
          printf("box_dim_x = %d\n", ptr_ch->box_dim_x);
          printf("box_dim_y = %d\n", ptr_ch->box_dim_y);
          printf("box_dim_z = %d\n", ptr_ch->box_dim_z);
          printf("box_real_dim_x = %d\n", ptr_ch->box_real_dim_x);
          printf("box_real_dim_y = %d\n", ptr_ch->box_real_dim_y);
          printf("box_real_dim_z = %d\n", ptr_ch->box_real_dim_z);
          printf("box_ts_x = %d\n", ptr_ch->box_ts_x);
          printf("box_ts_y = %d\n", ptr_ch->box_ts_y);
          printf("box_ts_z = %d\n", ptr_ch->box_ts_z);
          printf("box_min_x = %d, box_max_x = %d\n", ptr_ch->box_min_x, ptr_ch->box_max_x);
          printf("box_min_y = %d, box_max_y = %d\n", ptr_ch->box_min_y, ptr_ch->box_max_y);
          printf("box_min_z = %d, box_max_z = %d\n", ptr_ch->box_min_z, ptr_ch->box_max_z);

          printf("true = %d, false = %d\n", true, false);

          printf("ptr_ch->box_check_fit = %d\n", ptr_ch->box_check_fit);
          printf("ptr_ch->sim_bdry_contact = %d\n", ptr_ch->sim_bdry_contact);
          printf("ptr_ch->sim_bdry_contact_x = %d\n", ptr_ch->sim_bdry_contact_x);
          printf("ptr_ch->sim_bdry_contact_y = %d\n", ptr_ch->sim_bdry_contact_y);
          printf("ptr_ch->sim_bdry_contact_z = %d\n", ptr_ch->sim_bdry_contact_z);
          printf("ptr_ch->pbc_crosses_sim_box_bdry = %d\n", ptr_ch->pbc_crosses_sim_box_bdry);
          printf("ptr_ch->pbc_crosses_sim_box_bdry_x = %d\n", ptr_ch->pbc_crosses_sim_box_bdry_x);
          printf("ptr_ch->pbc_crosses_sim_box_bdry_y = %d\n", ptr_ch->pbc_crosses_sim_box_bdry_y);
          printf("ptr_ch->pbc_crosses_sim_box_bdry_z = %d\n", ptr_ch->pbc_crosses_sim_box_bdry_z);
          printf("ptr_ch->pbc_crosses_whole_sim_box = %d\n", ptr_ch->pbc_crosses_whole_sim_box);
          printf("ptr_ch->pbc_crosses_whole_sim_box_x = %d\n", ptr_ch->pbc_crosses_whole_sim_box_x);
          printf("ptr_ch->pbc_crosses_whole_sim_box_y = %d\n", ptr_ch->pbc_crosses_whole_sim_box_y);
          printf("ptr_ch->pbc_crosses_whole_sim_box_z = %d\n", ptr_ch->pbc_crosses_whole_sim_box_z);
        }
      }
    }
  }

  // Parent
  for (int cell_idx = 0; cell_idx < ptr_node->cell_size; cell_idx++)
  {
    box_idx_node = ptr_node->ptr_box_idx[cell_idx];
    for (int j = 0; j < ptr_node->ptr_cell_struct[box_idx_node].ptcl_size; j++)
    {
      ptcl_idx = ptr_node->ptr_cell_struct[box_idx_node].ptr_ptcl[j];
      if (box_idx_node != ptcl_idx_to_box_idx(ptr_node, ptcl_idx))
      {
        if (type == 1)
        {
          printf("\npre check\n");
        }

        else
        {
          printf("\npost check\n");
        }
        printf("Parent node at lv = %d, ID = %d\n", ptr_node->lv, ptr_node->ID);
        printf("Error, Parent, particle is out of the cell\n");
      }
    }
  }

  // Local mass and loca no particles:
  //  Children
  vtype aux_mass1;
  vtype aux_mass2;
  int aux_no_ptcl1;
  int aux_no_ptcl2;
  for (int i = 0; i < ptr_node->chn_size; i++)
  {
    aux_mass1 = 0;
    aux_mass2 = 0;
    aux_no_ptcl1 = 0;
    aux_no_ptcl2 = 0;
    ptr_ch = ptr_node->pptr_chn[i];
    for (int cell_idx = 0; cell_idx < ptr_ch->cell_size; cell_idx++)
    {
      box_idx_ch = ptr_ch->ptr_box_idx[cell_idx];

      aux_mass1 += ptr_ch->ptr_cell_struct[box_idx_ch].cell_mass;
      aux_no_ptcl1 += ptr_ch->ptr_cell_struct[box_idx_ch].ptcl_size;
    }
    for (int cell_idx = 0; cell_idx < ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y * ptr_ch->box_real_dim_z; cell_idx++)
    {
      aux_mass2 += ptr_ch->ptr_cell_struct[cell_idx].cell_mass;
      aux_no_ptcl2 += ptr_ch->ptr_cell_struct[cell_idx].ptcl_size;
    }
    if (aux_mass1 != aux_mass2 || aux_mass1 != ptr_ch->node_mass)
    {
      if (type == 1)
      {
        printf("pre check\n");
      }

      else
      {
        printf("post check\n");
      }

      printf("child : lv = %d, ID = %d\n", ptr_ch->lv, ptr_ch->ID);
      printf("\nError, Child local mass diferente a mass en caja\n ");
      printf("cell mass = %f, box mass = %f, local mass = %f\n", (double)aux_mass1, (double)aux_mass2, (double)ptr_ch->node_mass);

      printf("cell size = %d, zone size * 8 = %d\n", ptr_ch->cell_size, 8 * ptr_node->ptr_zone_size[i]);

      printf("child lv = %d, child ID = %d\n", ptr_ch->lv, ptr_ch->ID);
      printf("local ptcl = %d, local mass = %f\n", ptr_ch->no_ptcl_full_node, (double)ptr_ch->node_mass);

      printf("box_min_x = %d, box_max_x = %d\n", ptr_ch->box_min_x, ptr_ch->box_max_x);
      printf("box_min_y = %d, box_max_y = %d\n", ptr_ch->box_min_y, ptr_ch->box_max_y);
      printf("box_min_z = %d, box_max_z = %d\n", ptr_ch->box_min_z, ptr_ch->box_max_z);

      printf("box_real_dim_x = %d\n", ptr_ch->box_real_dim_x);
      printf("box_real_dim_y = %d\n", ptr_ch->box_real_dim_y);
      printf("box_real_dim_z = %d\n", ptr_ch->box_real_dim_z);

      printf("box_ts_x = %d\n", ptr_ch->box_ts_x);
      printf("box_ts_y = %d\n", ptr_ch->box_ts_y);
      printf("box_ts_z = %d\n", ptr_ch->box_ts_z);
    }

    if (aux_no_ptcl1 != aux_no_ptcl2 || aux_no_ptcl1 != ptr_ch->no_ptcl_full_node)
    {
      printf("\nError, Child local number of particles diferente a local no ptcl en caja\n ");
      printf("cell no particles = %d, box no ptcl = %d, local no ptcl = %d\n", (int)aux_no_ptcl1, (int)aux_no_ptcl2, (int)ptr_ch->no_ptcl_full_node);
    }
  }

  // Parent
  aux_mass1 = 0;
  aux_mass2 = 0;
  aux_no_ptcl1 = 0;
  aux_no_ptcl2 = 0;
  for (int cell_idx = 0; cell_idx < ptr_node->cell_size; cell_idx++)
  {
    box_idx_node = ptr_node->ptr_box_idx[cell_idx];

    aux_mass1 += ptr_node->ptr_cell_struct[box_idx_node].cell_mass;
    aux_no_ptcl1 += ptr_node->ptr_cell_struct[box_idx_node].ptcl_size;
  }
  for (int cell_idx = 0; cell_idx < ptr_node->box_real_dim_x * ptr_node->box_real_dim_y * ptr_node->box_real_dim_z; cell_idx++)
  {
    aux_mass2 += ptr_node->ptr_cell_struct[cell_idx].cell_mass;
    aux_no_ptcl2 += ptr_node->ptr_cell_struct[cell_idx].ptcl_size;
  }
  if (aux_mass1 != aux_mass2 || aux_mass1 != ptr_node->node_mass)
  {
    printf("\nError, Parent local mass diferente a mass en caja\n ");
    printf("cell mass = %f, box mass = %f, local mass = %f\n", (double)aux_mass1, (double)aux_mass2, (double)ptr_node->node_mass);
  }
  if (aux_no_ptcl1 != aux_no_ptcl2 || aux_no_ptcl1 != ptr_node->no_ptcl_full_node)
  {
    printf("\nError, Parent local number of particles diferente a local no ptcl en caja\n ");
    printf("cell no particles = %d, box no ptcl = %d, local no ptcl = %d\n", (int)aux_no_ptcl1, (int)aux_no_ptcl2, (int)ptr_node->no_ptcl_full_node);
  }

  // Boxes:
  //   Children
  int cntr_Periodic_boundary;
  int cntr_bdry_simulation;
  int cntr_No_exist;
  int cntr_Exist;
  int cntr_chn;

  for (int i = 0; i < ptr_node->chn_size; i++)
  {
    ptr_ch = ptr_node->pptr_chn[i];
    cntr_Periodic_boundary = 0;
    cntr_bdry_simulation = 0;
    cntr_No_exist = 0;
    cntr_Exist = 0;
    cntr_chn = 0;

    for (int cell_idx = 0; cell_idx < ptr_ch->cell_size; cell_idx++)
    {
      box_idx_ch = ptr_ch->ptr_box_idx[cell_idx];

      if (ptr_ch->ptr_box[box_idx_ch] < -3 || ptr_ch->ptr_box[box_idx_ch] == -1 || ptr_ch->ptr_box[box_idx_ch] == -2 || ptr_ch->ptr_box[box_idx_ch] > ptr_ch->chn_size - 1)
      {
        printf("\nError, hijo, las celdas de la caja no concuerdan\n ");
        printf("children size = %d, box_idx_ch = %d, box value = %d\n", ptr_ch->chn_size, box_idx_ch, ptr_ch->ptr_box[box_idx_ch]);
      }
    }

    for (int cell_idx = 0; cell_idx < ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y * ptr_ch->box_real_dim_z; cell_idx++)
    {
      if (ptr_ch->ptr_box[cell_idx] < -6)
      {
        printf("Error, celda < -6\n");
      }
      else if (ptr_ch->ptr_box[cell_idx] == -6)
      {
        cntr_Periodic_boundary++;
      }
      else if (ptr_ch->ptr_box[cell_idx] == -5)
      {
        cntr_bdry_simulation++;
      }
      else if (ptr_ch->ptr_box[cell_idx] == -4)
      {
        cntr_No_exist++;
      }
      else if (ptr_ch->ptr_box[cell_idx] == -3)
      {
        cntr_Exist++;
      }
      else if (ptr_ch->ptr_box[cell_idx] == -2 || ptr_ch->ptr_box[cell_idx] == -1)
      {
        printf("Error, Hijo, celda = %d\n", ptr_ch->ptr_box[cell_idx]);
      }
      else
      {
        cntr_chn++;
      }
    }

    if (ptr_ch->cell_size != cntr_Exist + cntr_chn)
    {
      printf("child lv = %d, child ID = %d\n", ptr_ch->lv, ptr_ch->ID);
      printf("local ptcl = %d, local mass = %f\n", ptr_ch->no_ptcl_full_node, (double)ptr_ch->node_mass);
      printf("Error, Hijo, cell size diferente la valor \n");
      printf("cell size = %d, exist + cntr_chn = %d\n", ptr_ch->cell_size, cntr_Exist + cntr_chn);
      printf("status -6 = %d\n", cntr_Periodic_boundary);
      printf("status -5 = %d\n", cntr_bdry_simulation);
      printf("status -4 = %d\n", cntr_No_exist);
      printf("status -3 = %d\n", cntr_Exist);
      printf("status >=0 = %d\n", cntr_chn);

      printf("box_min_x = %d, box_max_x = %d\n", ptr_ch->box_min_x, ptr_ch->box_max_x);
      printf("box_min_y = %d, box_max_y = %d\n", ptr_ch->box_min_y, ptr_ch->box_max_y);
      printf("box_min_z = %d, box_max_z = %d\n", ptr_ch->box_min_z, ptr_ch->box_max_z);

      printf("box_real_dim_x = %d\n", ptr_ch->box_real_dim_x);
      printf("box_real_dim_y = %d\n", ptr_ch->box_real_dim_y);
      printf("box_real_dim_z = %d\n", ptr_ch->box_real_dim_z);

      printf("box_ts_x = %d\n", ptr_ch->box_ts_x);
      printf("box_ts_y = %d\n", ptr_ch->box_ts_y);
      printf("box_ts_z = %d\n", ptr_ch->box_ts_z);

      if (type == 1)
        printf("pre check\n");
      else
        printf("post check\n");
    }
  }

  //   Parent
  cntr_Periodic_boundary = 0;
  cntr_bdry_simulation = 0;
  cntr_No_exist = 0;
  cntr_Exist = 0;
  cntr_chn = 0;
  for (int cell_idx = 0; cell_idx < ptr_node->cell_size; cell_idx++)
  {
    box_idx_node = ptr_node->ptr_box_idx[cell_idx];

    if (ptr_node->ptr_box[box_idx_node] < -3 || ptr_node->ptr_box[box_idx_node] == -1 || ptr_node->ptr_box[box_idx_node] == -2 || ptr_node->ptr_box[box_idx_node] > ptr_node->chn_size - 1)
    {
      printf("Error, Parent, las celdas de la caja no concuerdan\n ");
    }
  }

  for (int cell_idx = 0; cell_idx < ptr_node->box_real_dim_x * ptr_node->box_real_dim_y * ptr_node->box_real_dim_z; cell_idx++)
  {
    if (ptr_node->ptr_box[cell_idx] < -6)
    {
      printf("Error, celda < -6\n");
    }
    else if (ptr_node->ptr_box[cell_idx] == -6)
    {
      cntr_Periodic_boundary++;
    }
    else if (ptr_node->ptr_box[cell_idx] == -5)
    {
      cntr_bdry_simulation++;
    }
    else if (ptr_node->ptr_box[cell_idx] == -4)
    {
      cntr_No_exist++;
    }
    else if (ptr_node->ptr_box[cell_idx] == -3)
    {
      cntr_Exist++;
    }
    else if (ptr_node->ptr_box[cell_idx] == -2 || ptr_node->ptr_box[cell_idx] == -1)
    {
      printf("Error, celda = %d\n", ptr_node->ptr_box[cell_idx]);
    }
    else
    {
      cntr_chn++;
    }
  }
  if (ptr_node->cell_size != cntr_Exist + cntr_chn)
  {
    printf("Error, parent, cell size diferente la valor \n");
    if (type == 1)
      printf("pre check\n");
    else
      printf("post check\n");
  }

  // Minimum amount of cell size
  // if (ptr_node->cell_size < 64)
  // {
  //     printf("Error, Parent, cell size minimum amount of 64 broken = %d \n", ptr_node->cell_size);
  // }

  for (int i = 0; i < ptr_node->chn_size; i++)
  {
    ptr_ch = ptr_node->pptr_chn[i];
    cntr_chn = 0;

    for (int cell_idx = 0; cell_idx < ptr_node->box_real_dim_x * ptr_node->box_real_dim_y * ptr_node->box_real_dim_z; cell_idx++)
    {
      if (ptr_node->ptr_box[cell_idx] == i)
      {
        cntr_chn++;
      }
    }

    if (cntr_chn * 8 != ptr_ch->cell_size)
    {
      printf("Error, Parent, box hijo diferente cell size hijo \n");
    }
  }

  if (ptr_node->lv == lmin)
  {
    if (ptr_node->node_mass != GL_no_ptcl_final * 100.0)
    {
      printf("Error, Head lv, diferent mass\n");
      printf("local mass = %f\n", (double)ptr_node->node_mass);
    }

    if (ptr_node->no_ptcl_full_node != GL_no_ptcl_final)
    {
      printf("Error, Head lv, diferent local no ptcl\n");
      printf("local no_ptcl = %f\n", (double)ptr_node->no_ptcl_full_node);
    }

    if (GL_tentacles_size[0] != 1)
    {
      printf("Error, tentacles size != 1 at Head lv\n");
    }
  }

  // Checking grandchildren
  if (ptr_node->lv == lmax - 1)
  {
    for (int i = 0; i < ptr_node->chn_size; i++)
    {
      if (ptr_node->pptr_chn[i]->chn_size > 0)
      {
        printf("Error, there are grandchildren at level lmax -1\n");
      }
    }
  }

  // Checking head node
  if (ptr_node->lv == lmin)
  {
    if (ptr_node->cell_size != 1 << (3 * lmin))
    {
      printf("Error, head node cell size different to initial cell size\n");
    }
  }
}

/**
 * @brief Find the minimum and maximum cell indices of the subzone \e zone_idx.
 *
 * **SHORT DESCRIPTION**: Find the minimum and maximum cell indices in the 
 * *Code Space* (see Key Concepts \ref Key_Concepts_Code_Space "Code Space") 
 * that defines the *Smallest Box* (see Key Concepts \ref 
 * Key_Concepts_Smallest_Box "Smallest Box") associated with the corresponding 
 * subzone \e zone_idx.
 *
 * **PREREQUISITES**: Only used if periodic boundary conditions (pbc) are used.
 * Only called if the new refinement map of the parent node crosses the boundary
 * of the simulation.
 *
 * @param[in,out] ptr_node Pointer to node structure
 *
 * @param[in] zone_idx Idientification of some refinement zone
 *
 * **RETURN**: The error status.
 *
 * **LONG DESCRIPTION**:
 *
 * Find the minimum and maximum cell indices in the *Code Space* (see Key 
 * Concepts \ref Key_Concepts_Code_Space "Code Space") that defines the 
 * *Smallest Box* (see Key Concepts \ref Key_Concepts_Smallest_Box 
 * "Smallest Box") ssociated with the corresponding subzone \e zone_idx., 
 * storing them in the parent node structure parameters \ref
 * node::ptr_pbc_min_subzones_x "ptr_pbc_min_subzones_x" (\ref
 * node::ptr_pbc_min_subzones_y "y", \ref node::ptr_pbc_min_subzones_z "z"), and
 * \ref node::ptr_pbc_max_subzones_x "ptr_pbc_max_subzones_x" (\ref
 * node::ptr_pbc_max_subzones_y "y", \ref node::ptr_pbc_max_subzones_z "z").
 *
 * The flux of this function can be seen in the figure (work in progress), and
 * it is explained below:
 *
 * - [0]  <b> THE find_min_max_subzones_ref_PERIODIC_BOUNDARY() FUNCTION STARTS....</b>
 *
 * - [1]  Defining some internal useful parameters.
 *
 * - [2]  The current status of the new box contains the new map of refinement,
 * which will correspond to the final one. In this step, to find the subzones,
 * the cells of the zone of refinement with ID \e zone_idx are returned to the
 * previous status to \e refinement_required (-1) in the new box.
 *
 * - [3]  How we are going to create the subzones is very similar to the
 * creation of zones of refinement. However, there are some differences, for
 * example, there is a difference in the criterion of the definition of zone and
 * subzone of refinement. The first one accepts that a refinement zone crosses
 * the box of the simulation, while the second one does not.
 *
 * - [4]  A "while" loop runs over all cells in the zone of refinement of ID
 * \e zone_idx that have not been tagged to any subzone. A counter is used to
 * perform this task.
 *
 * - [5]  In every step of the "while" loop over the cells, we ask if the cell
 * belongs to a subzone of refinement (box value \f$ \geq 0 \f$). If it belongs,
 * we continue asking to the next cell, if it does not belong (box value
 * \f$ = -1 \f$), we change its box satus to "zone_idx" and pass to the next
 * step.
 *
 * - [6]  Having found a cell with no subzone, we are going to create an entire
 * new subzone starting with this cell as the foundation stone. To perform this,
 * a new "while" loop is executed running until there are no more cells without
 * analysis in the subzone, i.e. the final state of the block found is a solid
 * piece completely isolated from the other subzones, labeled with the
 * "subzone_idx" value in the parent node box.
 *
 * - [7]  At this point both "while" loops of the steps [4] and [6] end. Now,
 * using the box with this new information the minimum and maximum of every
 * subzone of the zone \e zone_idx are found.
 *
 * - [8]  Finally, the box status of the parent node is returned to its initial
 * state, putting the value of \e zone_idx in the corresponding cell of the new
 * box.
 *
 * - [9]  <b> THE find_min_max_subzones_ref_PERIODIC_BOUNDARY() FUNCTION ENDS....</b>
 *
 * **ILLUSTRATIVE EXAMPLES**:
 * - [a]  Trivial.
 *
 * **RATIONALES**:
 * - [a]  In step [2], we decide to use the same box to store the
 * information about the new refinement zones. It requires overwriting the
 * current information and then returning the old one. We choose this option
 * instead of creating a copy of the box because we are avoiding to use malloc,
 * memcpy, and free functions.
 *
 * **NOTES**:
 * - [a]  If the number of subzones found is equal to 1, then it is not
 * necessary to find the minimum and maximum cell indices, because they  are
 * equal to the minimum and maximum of the refinement level. It is  also valid
 * when more than one dimension crosses the box of the simulation.
 */

int static find_min_max_subzones_ref_PERIODIC_BOUNDARY(struct node *ptr_node, int zone_idx)
{
  //* >> Filling the auxiliary diferent refinement subzones *//
  int cntr_cell_add_all_subzones = 0; // Counter Number of cells added to any refinement subzone
  int cntr_cell_add;                  // Counter number of cells added per inspection
  int cntr_insp;                      // Numer of inspected cells in the subzone

  int subzone_size;    // Number of cells in the subzone
  int subzone_idx = 0; // Index of the subzone. Initialized at subzone 0

  int box_idx_node; // Box index

  int box_idxNbr_x_plus;  // Box index in the neigborhood on the right
  int box_idxNbr_x_minus; // Box index in the neigborhood on the left
  int box_idxNbr_y_plus;  // Box index in the neigborhood behind
  int box_idxNbr_y_minus; // Box index in the neigborhood in front
  int box_idxNbr_z_plus;  // Box index in the neigborhood up
  int box_idxNbr_z_minus; // Box index in the neigborhood down

  int cell_idx;         // The cell index is simply i of the for loop
  int cell_ref_idx = 0; // The index in the cell refined array pptr_zones[zone_idx]

  int box_real_dim_X_node = ptr_node->box_real_dim_x;
  int box_real_dim_X_times_Y_node = ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;

  // The auxiliary array ptr_aux_idx will be used to store the box indices of the subzone

  // Initializing the box in every zone_idx cell at value of REFINEMENT REQUIRED (-1)
  for (int i = 0; i < ptr_node->ptr_zone_size[zone_idx]; i++)
  {
    cell_idx = ptr_node->pptr_zones[zone_idx][i];
    box_idx_node = ptr_node->ptr_box_idx[cell_idx];
    ptr_node->ptr_box[box_idx_node] = -1;
  }

  //* >>  Changing the box status from REFINEMENT REQUIRED (-1) to the refinement subzone ID (>= 0) *//
  while (ptr_node->ptr_zone_size[zone_idx] > cntr_cell_add_all_subzones) // The loop while works as long as the number of cell addeed is less than the total refined cells
  {
    // printf("Inside of the "while" loop\n");
    //  Notes that the initiality we inspect the elements in the refined cell array until an element has been found that is not found in any of the current refinement subzones
    cell_idx = ptr_node->pptr_zones[zone_idx][cell_ref_idx]; // Index of the cells array in the node
    box_idx_node = ptr_node->ptr_box_idx[cell_idx];

    if (ptr_node->ptr_box[box_idx_node] == -1) // A cell without subzone has been founded
    {
      // printf("inside of the if ptr_box [] == -1\n");
      subzone_size = 0; // Initial number of element in the subzone

      //* >> Including the first element of the box to the auxiliary array ptr_aux_idx *//
      ptr_node->ptr_aux_idx[0] = box_idx_node;

      //* >>  Changing the box status from zone ID (>= 0) to the REFINEMENT REQUIRED (-1) in order not to repeat elements *//
      ptr_node->ptr_box[box_idx_node] = subzone_idx;

      subzone_size++;               // +1 to the number of cells in the zone
      cntr_cell_add_all_subzones++; // +1 to the number of cells added in total

      //* >> Building of the subzone of ID = subzone_idx *//
      cntr_insp = 0;                                                                           // Counter for the number of elements inspectioned in the current subzone array
      while (subzone_size > cntr_insp && ptr_node->cell_ref_size > cntr_cell_add_all_subzones) // The loop ends when all the elements of the subzone have been inspected or when the number of cells in the zone is equal to the total number of cells added in all the subzones.
      {
        // Note that the number of elements in the subzone increases if the neighbors of the inspected cell must be added to the subzone

        cntr_cell_add = 0; // Counter number of cells added per cell inspected
        box_idx_node = ptr_node->ptr_aux_idx[cntr_insp];

        box_idxNbr_x_plus = box_idx_node + 1;
        box_idxNbr_x_minus = box_idx_node - 1;
        box_idxNbr_y_plus = box_idx_node + box_real_dim_X_node;
        box_idxNbr_y_minus = box_idx_node - box_real_dim_X_node;
        box_idxNbr_z_plus = box_idx_node + box_real_dim_X_times_Y_node;
        box_idxNbr_z_minus = box_idx_node - box_real_dim_X_times_Y_node;

        //* Checking the nearest 6 neighbors of face
        // First neighbor
        if (ptr_node->ptr_box[box_idxNbr_x_plus] == -1)
        {
          ptr_node->ptr_aux_idx[subzone_size + cntr_cell_add] = box_idxNbr_x_plus; // Including the neighboring element of the box to the auxiliary array
          ptr_node->ptr_box[box_idxNbr_x_plus] = subzone_idx;                      // Changing the box status from refinement zone ID (>= 0) to the REFINEMENT REQUIRED (-1)
          cntr_cell_add++;                                                         // +1 to the number of cells added in the current inspection
        }
        // Second neighbor
        if (ptr_node->ptr_box[box_idxNbr_x_minus] == -1)
        {
          ptr_node->ptr_aux_idx[subzone_size + cntr_cell_add] = box_idxNbr_x_minus; // Including the neighboring element of the box to the auxiliary array
          ptr_node->ptr_box[box_idxNbr_x_minus] = subzone_idx;                      // Changing the box status from refinement zone ID (>= 0) to the REFINEMENT REQUIRED (-1)
          cntr_cell_add++;                                                          // +1 to the number of cells added in the current inspection
        }
        // Third neighbor
        if (ptr_node->ptr_box[box_idxNbr_y_plus] == -1)
        {
          ptr_node->ptr_aux_idx[subzone_size + cntr_cell_add] = box_idxNbr_y_plus; // Including the neighboring element of the box to the auxiliary array
          ptr_node->ptr_box[box_idxNbr_y_plus] = subzone_idx;                      // Changing the box status from refinement zone ID (>= 0) to the REFINEMENT REQUIRED (-1)
          cntr_cell_add++;                                                         // +1 to the number of cells added in the current inspection
        }
        // Fourth neighbor
        if (ptr_node->ptr_box[box_idxNbr_y_minus] == -1)
        {
          ptr_node->ptr_aux_idx[subzone_size + cntr_cell_add] = box_idxNbr_y_minus; // Including the neighboring element of the box to the auxiliary array
          ptr_node->ptr_box[box_idxNbr_y_minus] = subzone_idx;                      // Changing the box status from refinement zone ID (>= 0) to the REFINEMENT REQUIRED (-1)
          cntr_cell_add++;                                                          // +1 to the number of cells added in the current inspection
        }
        // Fifth neighbor
        if (ptr_node->ptr_box[box_idxNbr_z_plus] == -1)
        {
          ptr_node->ptr_aux_idx[subzone_size + cntr_cell_add] = box_idxNbr_z_plus; // Including the neighboring element of the box to the auxiliary array
          ptr_node->ptr_box[box_idxNbr_z_plus] = subzone_idx;                      // Changing the box status from refinement zone ID (>= 0) to the REFINEMENT REQUIRED (-1)
          cntr_cell_add++;                                                         // +1 to the number of cells added in the current inspection
        }
        // Sixth neighbor
        if (ptr_node->ptr_box[box_idxNbr_z_minus] == -1)
        {
          ptr_node->ptr_aux_idx[subzone_size + cntr_cell_add] = box_idxNbr_z_minus; // Including the neighboring element of the box to the auxiliary array
          ptr_node->ptr_box[box_idxNbr_z_minus] = subzone_idx;                      // Changing the box status from refinement zone ID (>= 0) to the REFINEMENT REQUIRED (-1)
          cntr_cell_add++;                                                          // +1 to the number of cells added in the current inspection
        }
        subzone_size += cntr_cell_add;               // Increasing the number of cells in the subzone
        cntr_cell_add_all_subzones += cntr_cell_add; // Increasing the number of cells added to the subzones
        cntr_insp++;                                 // Increasing the number of inspections
      }                                              // End "while" loop, now the box contains the the information about all cells of the subzone "subzone_idx"
      subzone_idx++;                                 // Increasing the subzone number
    }                                                // subZone defined in the box
    cell_ref_idx++;
  } // At this point the box contains the information of all refinement subzones

  // ptr_node->pbc_subzones_cap = subzone_idx_max; // Maximum amount of subzones
  ptr_node->pbc_subzones_size = subzone_idx; // Total amount of subzones

  if (ptr_node->pbc_subzones_size > 1)
  {
    //* >> Space checking of refinement subzones min and max arrays  *//
    if (space_check(&(ptr_node->pbc_subzones_cap), subzone_idx, 2.0f, "p6i1i1i1i1i1i1", &(ptr_node->ptr_pbc_min_subzones_x), &(ptr_node->ptr_pbc_max_subzones_x), &(ptr_node->ptr_pbc_min_subzones_y), &(ptr_node->ptr_pbc_max_subzones_y), &(ptr_node->ptr_pbc_min_subzones_z), &(ptr_node->ptr_pbc_max_subzones_z)) == _FAILURE_)
    {
      printf("Error, in space_check function\n");
      return _FAILURE_;
    }

    //* >> Initializing suzbones min and max
    if (ptr_node->ptr_pbc_bool_bdry_anomalies_x[zone_idx] == true)
    {
      for (int i = 0; i < subzone_idx; i++)
      {
        ptr_node->ptr_pbc_min_subzones_x[i] = INT_MAX;
        ptr_node->ptr_pbc_max_subzones_x[i] = INT_MIN;
      }
    }

    if (ptr_node->ptr_pbc_bool_bdry_anomalies_y[zone_idx] == true)
    {
      for (int i = 0; i < subzone_idx; i++)
      {
        ptr_node->ptr_pbc_min_subzones_y[i] = INT_MAX;
        ptr_node->ptr_pbc_max_subzones_y[i] = INT_MIN;
      }
    }

    if (ptr_node->ptr_pbc_bool_bdry_anomalies_z[zone_idx] == true)
    {
      for (int i = 0; i < subzone_idx; i++)
      {
        ptr_node->ptr_pbc_min_subzones_z[i] = INT_MAX;
        ptr_node->ptr_pbc_max_subzones_z[i] = INT_MIN;
      }
    }

    //* >> Adding the cells to the zone array pptr_zones *//
    for (int i = 0; i < ptr_node->ptr_zone_size[zone_idx]; i++)
    {
      cell_idx = ptr_node->pptr_zones[zone_idx][i];
      box_idx_node = ptr_node->ptr_box_idx[cell_idx];
      subzone_idx = ptr_node->ptr_box[box_idx_node];

      // Min and Max
      if (ptr_node->ptr_pbc_bool_bdry_anomalies_x[zone_idx] == true)
      {
        if (ptr_node->ptr_pbc_min_subzones_x[subzone_idx] > ptr_node->ptr_cell_idx_x[cell_idx])
        {
          ptr_node->ptr_pbc_min_subzones_x[subzone_idx] = ptr_node->ptr_cell_idx_x[cell_idx];
        }
        if (ptr_node->ptr_pbc_max_subzones_x[subzone_idx] < ptr_node->ptr_cell_idx_x[cell_idx])
        {
          ptr_node->ptr_pbc_max_subzones_x[subzone_idx] = ptr_node->ptr_cell_idx_x[cell_idx];
        }
      }

      if (ptr_node->ptr_pbc_bool_bdry_anomalies_y[zone_idx] == true)
      {
        if (ptr_node->ptr_pbc_min_subzones_y[subzone_idx] > ptr_node->ptr_cell_idx_y[cell_idx])
        {
          ptr_node->ptr_pbc_min_subzones_y[subzone_idx] = ptr_node->ptr_cell_idx_y[cell_idx];
        }
        if (ptr_node->ptr_pbc_max_subzones_y[subzone_idx] < ptr_node->ptr_cell_idx_y[cell_idx])
        {
          ptr_node->ptr_pbc_max_subzones_y[subzone_idx] = ptr_node->ptr_cell_idx_y[cell_idx];
        }
      }

      if (ptr_node->ptr_pbc_bool_bdry_anomalies_z[zone_idx] == true)
      {
        if (ptr_node->ptr_pbc_min_subzones_z[subzone_idx] > ptr_node->ptr_cell_idx_z[cell_idx])
        {
          ptr_node->ptr_pbc_min_subzones_z[subzone_idx] = ptr_node->ptr_cell_idx_z[cell_idx];
        }
        if (ptr_node->ptr_pbc_max_subzones_z[subzone_idx] < ptr_node->ptr_cell_idx_z[cell_idx])
        {
          ptr_node->ptr_pbc_max_subzones_z[subzone_idx] = ptr_node->ptr_cell_idx_z[cell_idx];
        }
      }
      ptr_node->ptr_box[box_idx_node] = zone_idx; // Returning the value zone_idx at the box in the zone partitioned
    }
  }
  else
  {
    // Returning the value zone_idx at the box in the zone partitioned
    for (int i = 0; i < ptr_node->ptr_zone_size[zone_idx]; i++)
    {
      cell_idx = ptr_node->pptr_zones[zone_idx][i];
      box_idx_node = ptr_node->ptr_box_idx[cell_idx];
      ptr_node->ptr_box[box_idx_node] = zone_idx;
    }
  }

  return _SUCCESS_;
}

/**
 * @brief Update the cell structure of the parent node.
 *
 * **SHORT DESCRIPTION**: Update the cells which currently belong to any
 * refinement zone of the parent node using the information of its child nodes.
 *
 * **PREREQUISITES**: Always used.
 *
 * @param[in,out] ptr_node Pointer to node structure
 *
 * **RETURN**: The error status.
 *
 * **LONG DESCRIPTION**:
 *
 * At the beginning of the tree_adaptaion.c module, the information of the cells
 * of the parent node which do not belong to any refinement zone have the right
 * information. But, the cells of the parent node which are in a refinement
 * zone (so they are in a child node), can have incorrect information in the
 * cell structure of the parent node. This function is in charge of update the
 * information of the cells which belong to any refinement zone using the
 * child nodes of the parent node. 
 * 
 * The flux of this function can be seen in the figure (work in progress), and 
 * it is explained below:
 *
 * - [0]  <b> THE updating_cell_struct() FUNCTION STARTS....</b>
 *
 * - [1]  Defining some internal useful parameters.
 *
 * - [2]  Run a "for" loop over the child nodes.
 *
 * - [3]  Run a "for" loop over the cells in every child node. The loop runs
 * each 8 cells which correspond to one cell in the parent node.
 *
 * - [4]  The parent box index associated to this package of 8 child cells is
 * computed.
 *
 * - [5]  The total mass and particles are copied from the 8 child cells to the
 * respective parent node cell.
 *
 * - [6]  The "for" loops ends.
 *
 * - [7]  <b> THE updating_cell_struct() FUNCTION ENDS....</b>
 *
 * **ILLUSTRATIVE EXAMPLES**:
 * - [a]  Trivial.
 *
 * **RATIONALES**:
 * - [a]  a
 *
 * **NOTES**:
 * - [a]  In this module, this function is one of the most expensive. We are 
 * working to find a more optimized way to perform this function.
 */

static int updating_cell_struct(struct node *ptr_node)
{
  int no_chn;          // Number of child nodes
  struct node *ptr_ch; // child node
  vtype aux_cell_mass; // Total mass in the 8 child cells
  int aux_no_ptcl;     // Total mass in the 8 child cells

  int box_idx_ch; // Box index of the child cell

  int box_idx_x_node; // Box index in X direcction of the node cell
  int box_idx_y_node; // Box index in Y direcction of the node cell
  int box_idx_z_node; // Box index in Z direcction of the node cell
  int box_idx_node;   // Box index of the node cell

  int *ptr_int_aux;

  no_chn = ptr_node->chn_size;

  int box_real_dim_X_node = ptr_node->box_real_dim_x;
  int box_real_dim_X_times_Y_node = ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;

  int aux_idx_x;
  int aux_idx_y;
  int aux_idx_z;

  int lv = ptr_node->lv;

  for (int i = 0; i < no_chn; i++) // Loop over children
  {
    ptr_ch = ptr_node->pptr_chn[i];
    for (int j = 0; j < ptr_ch->cell_size; j += 8) // Loop over packeges of 8 cells
    {
      //* >> Computing the mass of the 8 child cells *//

      aux_idx_x = (ptr_ch->ptr_cell_idx_x[j] >> 1);
      aux_idx_y = (ptr_ch->ptr_cell_idx_y[j] >> 1);
      aux_idx_z = (ptr_ch->ptr_cell_idx_z[j] >> 1);

      box_idx_x_node = aux_idx_x - ptr_node->box_ts_x;
      box_idx_y_node = aux_idx_y - ptr_node->box_ts_y;
      box_idx_z_node = aux_idx_z - ptr_node->box_ts_z;

      if (ptr_node->pbc_crosses_sim_box_bdry == true)
      {
        if (aux_idx_x > ptr_node->box_max_x)
        {
          box_idx_x_node -= (1 << lv);
        }

        if (aux_idx_y > ptr_node->box_max_y)
        {
          box_idx_y_node -= (1 << lv);
        }

        if (aux_idx_z > ptr_node->box_max_z)
        {
          box_idx_z_node -= (1 << lv);
        }
      }

      box_idx_node = box_idx_x_node + box_idx_y_node * box_real_dim_X_node + box_idx_z_node * box_real_dim_X_times_Y_node;

      aux_cell_mass = 0;
      aux_no_ptcl = 0;

      //* >> Counting mass and number of particles *//
      for (int k = j; k < j + 8; k++)
      {
        box_idx_ch = ptr_ch->ptr_box_idx[k];
        // box_idx_ch = (box_idx_x_ch + ii) + (box_idx_y_ch + jj) * ptr_ch->box_real_dim_x + (box_idx_z_ch + kk) * ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;
        aux_cell_mass += ptr_ch->ptr_cell_struct[box_idx_ch].cell_mass;
        aux_no_ptcl += ptr_ch->ptr_cell_struct[box_idx_ch].ptcl_size;
      }

      ptr_node->ptr_cell_struct[box_idx_node].cell_mass = aux_cell_mass;
      ptr_node->ptr_cell_struct[box_idx_node].ptcl_size = aux_no_ptcl;

      //* >> Space checking of the capacity of the number of particles in the cell *//
      if (space_check(&(ptr_node->ptr_cell_struct[box_idx_node].ptcl_cap), aux_no_ptcl, 1.0f, "p1i1", &(ptr_node->ptr_cell_struct[box_idx_node].ptr_ptcl)) == _FAILURE_)
      {
        printf("Error, in space_check function\n");
        return _FAILURE_;
      }

      ptr_int_aux = ptr_node->ptr_cell_struct[box_idx_node].ptr_ptcl;
      //* >> Transferring particles *//
      for (int k = j; k < j + 8; k++)
      {
        box_idx_ch = ptr_ch->ptr_box_idx[k];
        if (ptr_ch->ptr_cell_struct[box_idx_ch].ptcl_size > 0)
        {
          memcpy(ptr_int_aux, ptr_ch->ptr_cell_struct[box_idx_ch].ptr_ptcl, ptr_ch->ptr_cell_struct[box_idx_ch].ptcl_size * sizeof(int));
          ptr_int_aux += ptr_ch->ptr_cell_struct[box_idx_ch].ptcl_size;
        }
      }
    }
  }

  return _SUCCESS_;
}

/**
 * @brief Store the box array \ref node::ptr_box "ptr_box" in the old_box array
 * \ref node::ptr_box_old "ptr_box_old", and reset the box array \ref
 * node::ptr_box "ptr_box".
 *
 * **SHORT DESCRIPTION**: Store the box array \ref node::ptr_box "ptr_box" in
 * the old_box array \ref node::ptr_box_old "ptr_box_old", and reset the cells
 * that require refinement \f$ \geq 0 \f$ to existing cells (=-3) in the box
 * array \ref node::ptr_box "ptr_box".
 *
 * **PREREQUISITES**: Always used.
 *
 * @param[in,out] ptr_node Pointer to node structure
 *
 * **RETURN**: No parameter is returned.
 *
 * **LONG DESCRIPTION**:
 *
 * Store the box array \ref node::ptr_box "ptr_box" in the old_box array \ref
 * node::ptr_box_old "ptr_box_old", and reset the cells that require refinement
 * \f$ \geq 0 \f$ to existing cells (=-3) in the box array \ref node::ptr_box
 * "ptr_box".
 *
 * The flux of this function can be seen in the figure (work in progress), and
 * it is explained below:
 *
 * - [0]  <b> THE initialization_node_boxes() FUNCTION STARTS....</b>
 *
 * - [1]  Defining some internal useful parameters.
 *
 * - [2]  Using the function memcpy() the box array \ref node::ptr_box "ptr_box"
 * is copied to the old box array \ref node::ptr_box_old "ptr_box_old".
 *
 * - [3]  Run a "for" loop over the child nodes.
 *
 * - [4]  Run a "for" loop over the cells in every child node. The loop runs
 * each 8 cells which correspond to one cell in the parent node.
 *
 * - [5]  The parent box index associated to this package of 8 child cells is
 *        computed.
 *
 * - [6]  The box array \ref node::ptr_box "ptr_box" is reset to the Exist
 * status (=-3).
 *
 * - [7]  The "for" loops ends.
 *
 * - [8]  <b> THE initialization_node_boxes() FUNCTION ENDS....</b>
 *
 * **ILLUSTRATIVE EXAMPLES**:
 * - [a]  Trivial.
 *
 * **RATIONALES**:
 * - [a]  a
 *
 * **NOTES**:
 * - [a]  a
 */

static void initialization_node_boxes(struct node *ptr_node)
{
  struct node *ptr_ch;

  int cap; // Maximum amount of cells in the box

  int box_idx_x_node; // Box index in X direcction
  int box_idx_y_node; // Box index in Y direcction
  int box_idx_z_node; // Box index in Z direcction
  int box_idx_node;   // Box index

  int aux_idx_x;
  int aux_idx_y;
  int aux_idx_z;

  int lv = ptr_node->lv;

  int box_real_dim_X_node = ptr_node->box_real_dim_x;
  int box_real_dim_X_times_Y_node = ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;

  cap = ptr_node->box_real_dim_x * ptr_node->box_real_dim_y * ptr_node->box_real_dim_z; // In general, the size of each side must be 3 times bigger than the same side of the "Smallest Box"

  memcpy(ptr_node->ptr_box_old, ptr_node->ptr_box, cap * sizeof(int));

  for (int i = 0; i < ptr_node->chn_size; i++)
  {
    ptr_ch = ptr_node->pptr_chn[i];

    for (int j = 0; j < ptr_ch->cell_size; j += 8)
    {

      aux_idx_x = (ptr_ch->ptr_cell_idx_x[j] >> 1);
      aux_idx_y = (ptr_ch->ptr_cell_idx_y[j] >> 1);
      aux_idx_z = (ptr_ch->ptr_cell_idx_z[j] >> 1);

      box_idx_x_node = aux_idx_x - ptr_node->box_ts_x;
      box_idx_y_node = aux_idx_y - ptr_node->box_ts_y;
      box_idx_z_node = aux_idx_z - ptr_node->box_ts_z;

      if (ptr_node->pbc_crosses_sim_box_bdry == true)
      {
        if (aux_idx_x > ptr_node->box_max_x)
        {
          box_idx_x_node -= (1 << lv);
        }

        if (aux_idx_y > ptr_node->box_max_y)
        {
          box_idx_y_node -= (1 << lv);
        }

        if (aux_idx_z > ptr_node->box_max_z)
        {
          box_idx_z_node -= (1 << lv);
        }
      }

      box_idx_node = box_idx_x_node + box_idx_y_node * box_real_dim_X_node + box_idx_z_node * box_real_dim_X_times_Y_node;
      ptr_node->ptr_box[box_idx_node] = -3;
    }
  }
}

/**
 * @brief Cells in the parent node are labeled as refinement cells in the \ref
 * node.ptr_box "ptr_box"
 *
 * **SHORT DESCRIPTION**: Cells in the parent node are labeled as refinement
 * cells in the \ref node.ptr_box "ptr_box" (-1) and added to the array \ref
 * node.ptr_cell_ref "ptr_cell_ref" as positional indices of the cell arrays.
 *
 * **PREREQUISITES**: Always used.
 *
 * @param[in,out] ptr_node Pointer to node structure
 *
 * **RETURN**: The error status.
 *
 * **LONG DESCRIPTION**:
 *
 * Cells in the parent node are labeled as refinement cells in the \ref
 * node.ptr_box "ptr_box", i.e. with status of -1 in the box, and added to the
 * array \ref node.ptr_cell_ref "ptr_cell_ref" as positional indices of the cell
 * arrays (for example \ref node.ptr_cell_idx_x "ptr_cell_idx_x"). There are two
 * criteria to decide if a cell requires refinement. A cell will be refined if
 * *(a)* the cell in the parent node contains a grandchild cell, or if *(b)*
 * the cell satisfies the refinement criteria given by the user. In any
 * situation, the neighboring cells of the refined cell should be also refined
 * according to the n_exp parameter.
 *
 * The flux of this function can be seen in the figure (work in progress), and
 * it is explained below:
 *
 * - [0]  <b> THE fill_cell_ref() FUNCTION STARTS....</b>
 *
 * - [1]  Defining some internal useful parameters.
 *
 * - [2]  Firstly the function asks about the first criteria *(a)* of
 * refinement if the cell contains a grandchild cell. To do that, a "for" loop
 * over the child nodes is performed
 *
 * - [3]  Run a "for" loop over the grandchild nodes.
 *
 * - [4]  Run a "for" loop over the cells in every grandchild node. The loop 
 * runs each 8 cells which correspond to one cell in the child node, and it
 * belongs to one cell in the corresponding parent node.
 *
 * - [5]  The parent box index which contains this package of 8 grandchild
 * cells is computed.
 *
 * - [6]  If the cell status hasn't been modified yet, it changes from "Exist"
 * (-3) to "Requires-Refinement" (-1) in the box index.
 *
 * - [7]  The neighboring cells at the distance of n_exp are also modified in
 * their box indices if they haven't been modified yet.
 *
 * - [8]  At this point the three "for" loops end, the "for" loop of the child
 * nodes, grandchild nodes, and the grandchild cells. Now, the function
 * continues with the second criterion, it will ask *(b)* if the cells satisfy
 * the refinement criteria given by the user.
 *
 * - [9]  To perform this later, a "for" loop over the parent node cell is
 * performed.
 *
 * - [10] If the cell satisfies the refinement criteria, and its status hasn't
 * been modified yet, it changes from "Exist" (-3) to "Requires-Refinement" (-1)
 * in the box index
 *
 * - [11] The neighboring cells at the distance of n_exp are also modified in
 * their box indices if they haven't been modified yet.
 *
 * - [12] At this point the "for" loop ends, and every cell in the node which
 * will require refinement has been labeled to do that.
 *
 * - [13] The next step is to add those cells to the auxiliary array of
 * positional cell indices \ref node.ptr_cell_ref "ptr_cell_ref". To do that, a
 * "for" loop running over the parent node cells is performed.
 *
 * - [14] If the box index associated to the cell has the status of
 * "Requires-Refinement" (-1), the positional index of the cell is added to the
 * array \ref node.ptr_cell_ref "ptr_cell_ref".
 *
 * - [15] The "for" loop over the parent node cells ends
 *
 * - [16] <b> THE fill_cell_ref() FUNCTION ENDS....</b>
 *
 * **ILLUSTRATIVE EXAMPLES**:
 * - [a]  Trivial.
 *
 * **RATIONALES**:
 * - [a]  The reason why to decide to add the elements to the array \ref
 * node.ptr_cell_ref "ptr_cell_ref", is because we are avoiding to ask if there
 * is enough Capacity (see Key Concepts \ref Key_Concepts_Capacity "Capacity") 
 * in the array.
 *
 * **NOTES**:
 * - [a]
 */

static int fill_cell_ref(struct node *ptr_node)
{
  //* >> Adding cells which satisfy the refinement criteria or contains grandchild cells to the array ptr_cell_ref and chaning the box_aux the status of refinement -1 *//

  struct node *ptr_ch;      // Child node
  struct node *ptr_grandch; // Grandchild node

  int cell_ref_idx = 0; // Index of the position in the cell refined array, Size of the refinement cells array. Initial size of the cell refined array must be 0

  int box_idxNbr_node; // Box index in the neigborhood

  int cntr; // Counter used to add cell index to the ptr_cell_ref

  int box_idx_x_node;
  int box_idx_y_node;
  int box_idx_z_node;
  int box_idx_node;

  int lv = ptr_node->lv;

  int aux_idx_x;
  int aux_idx_y;
  int aux_idx_z;

  int box_real_dim_X_node = ptr_node->box_real_dim_x;
  int box_real_dim_X_times_Y_node = ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;

  //* >> Adding to refine (-1) to cells in the node wich contains grandchildren refinement cells
  for (int ch = 0; ch < ptr_node->chn_size; ch++) // Loop over children
  {
    ptr_ch = ptr_node->pptr_chn[ch];

    for (int grandch = 0; grandch < ptr_ch->chn_size; grandch++)
    {
      ptr_grandch = ptr_ch->pptr_chn[grandch];

      for (int cell_grandch = 0; cell_grandch < ptr_grandch->cell_size; cell_grandch += 8)
      {

        aux_idx_x = ptr_grandch->ptr_cell_idx_x[cell_grandch] >> 2;
        aux_idx_y = ptr_grandch->ptr_cell_idx_y[cell_grandch] >> 2;
        aux_idx_z = ptr_grandch->ptr_cell_idx_z[cell_grandch] >> 2;

        box_idx_x_node = aux_idx_x - ptr_node->box_ts_x;
        box_idx_y_node = aux_idx_y - ptr_node->box_ts_y;
        box_idx_z_node = aux_idx_z - ptr_node->box_ts_z;

        if (ptr_node->pbc_crosses_sim_box_bdry == true)
        {
          if (aux_idx_x > ptr_node->box_max_x)
          {
            box_idx_x_node -= (1 << lv);
          }

          if (aux_idx_y > ptr_node->box_max_y)
          {
            box_idx_y_node -= (1 << lv);
          }

          if (aux_idx_z > ptr_node->box_max_z)
          {
            box_idx_z_node -= (1 << lv);
          }
        }

        box_idx_node = box_idx_x_node + box_idx_y_node * box_real_dim_X_node + box_idx_z_node * box_real_dim_X_times_Y_node;

        if (ptr_node->ptr_box[box_idx_node] == -3) // Cell has not been added yet
        {
          cell_ref_idx++; // +1 to the total number of cells to be refined
          //* >> Chaning the cell box status from EXIST (-3) to REFINEMENT REQUIRED (-1) *//
          ptr_node->ptr_box[box_idx_node] = -1;
        }

        for (int kk = -n_exp; kk < n_exp + 1; kk++)
        {
          for (int jj = -n_exp; jj < n_exp + 1; jj++)
          {
            for (int ii = -n_exp; ii < n_exp + 1; ii++)
            {
              box_idxNbr_node = box_idx_node + ii + jj * box_real_dim_X_node + kk * box_real_dim_X_times_Y_node;

              if (bdry_cond_type == 0 && ptr_node->ptr_box[box_idxNbr_node] == -6)
              {
                if (ptr_node->ptr_box[box_idx_node + ii] == -6)
                {
                  box_idxNbr_node += ii < 0 ? (1 << lv) : -(1 << lv);
                }

                if (ptr_node->ptr_box[box_idx_node + jj * box_real_dim_X_node] == -6)
                {
                  box_idxNbr_node += jj < 0 ? (1 << lv) * box_real_dim_X_node : -(1 << lv) * box_real_dim_X_node;
                }

                if (ptr_node->ptr_box[box_idx_node + kk * box_real_dim_X_times_Y_node] == -6)
                {
                  box_idxNbr_node += kk < 0 ? (1 << lv) * box_real_dim_X_times_Y_node : -(1 << lv) * box_real_dim_X_times_Y_node;
                }
              }

              //* >> Asking if the neighboring cell has not been changed yet *//
              // The border (-2) of the simulation can be added
              if (ptr_node->ptr_box[box_idxNbr_node] == -3) // Cell has not been added yet
              {
                cell_ref_idx++; // +1 to the total number of cells to be refined
                //* >> Chaning the cell box status from EXIST (-3) to REFINEMENT REQUIRED (-1) *//
                ptr_node->ptr_box[box_idxNbr_node] = -1;
              }
            }
          }
        }
      }
    }
  }

  //* >> Changing the box status from EXIST (-3) to REFINEMENT REQUIRED (-1) of cells satisfying the refinement criterion *//
  for (int i = 0; i < ptr_node->cell_size; i++)
  {
    box_idx_node = ptr_node->ptr_box_idx[i];

    // Refinement criterion in the box_mass in no border box points
    if (ptr_node->ptr_cell_struct[box_idx_node].cell_mass >= ref_criterion_mass || ptr_node->ptr_cell_struct[box_idx_node].ptcl_size >= ref_criterion_ptcl) // Including border (-2)
    {
      // if (ptr_node->ptr_box[box_idx_node] == -3) // Cell has not been added yet
      if (ptr_node->ptr_box[box_idx_node] == -3) // Cell has not been added yet
      {
        cell_ref_idx++;
        //* >> Chaning the cell box status from EXIST (-3) to REFINEMENT REQUIRED (-1) *//
        ptr_node->ptr_box[box_idx_node] = -1;
      }

      //* >> Changing the neighboring cell status *//
      for (int kk = -n_exp; kk < n_exp + 1; kk++)
      {
        for (int jj = -n_exp; jj < n_exp + 1; jj++)
        {
          for (int ii = -n_exp; ii < n_exp + 1; ii++)
          {
            box_idxNbr_node = box_idx_node + ii + jj * box_real_dim_X_node + kk * box_real_dim_X_times_Y_node;

            if (bdry_cond_type == 0 && ptr_node->ptr_box[box_idxNbr_node] == -6)
            {
              if (ptr_node->ptr_box[box_idx_node + ii] == -6)
              {
                box_idxNbr_node += ii < 0 ? (1 << lv) : -(1 << lv);
              }

              if (ptr_node->ptr_box[box_idx_node + jj * box_real_dim_X_node] == -6)
              {
                box_idxNbr_node += jj < 0 ? (1 << lv) * box_real_dim_X_node : -(1 << lv) * box_real_dim_X_node;
              }

              if (ptr_node->ptr_box[box_idx_node + kk * box_real_dim_X_times_Y_node] == -6)
              {
                box_idxNbr_node += kk < 0 ? (1 << lv) * box_real_dim_X_times_Y_node : -(1 << lv) * box_real_dim_X_times_Y_node;
              }
            }

            //* >> Asking if the neighboring cell has not been changed yet *//
            if (ptr_node->ptr_box[box_idxNbr_node] == -3) // Cell has not been added yet
            {
              cell_ref_idx++;
              //* >> Chaning the cell box status from EXIST (-3) to REFINEMENT REQUIRED (-1) *//
              ptr_node->ptr_box[box_idxNbr_node] = -1;
            }
          }
        }
      }
    }
  }

  //* >> Space checking of the capacity of the refined cells *//
  if (space_check(&(ptr_node->cell_ref_cap), cell_ref_idx, 1.0f, "p1i1", &(ptr_node->ptr_cell_ref)) == _FAILURE_)
  {
    printf("Error, in space_check function\n");
    return _FAILURE_;
  }

  //* >> Adding the infomation about size of the ptr_cell_ref array *//
  ptr_node->cell_ref_size = cell_ref_idx;

  //* >> Adding cells refined to the array ptr_cell_ref *//
  cntr = 0; // Counter for the position in the cell refined array
  for (int i = 0; i < ptr_node->cell_size; i++)
  {
    box_idx_node = ptr_node->ptr_box_idx[i];

    if (ptr_node->ptr_box[box_idx_node] == -1) // Cell require refinement
    {
      ptr_node->ptr_cell_ref[cntr] = i;
      cntr++;
    }
  }

  return _SUCCESS_;
} // end function fill_cell_ref

/**
 * @brief Cells in the parent node labeled as refinement cells are separated
 * into different refinement zones.
 *
 * **SHORT DESCRIPTION**: Cells in the parent node labeled as cells that require
 * refinement are separated in the new refinement zones of the parent node.
 *
 * **PREREQUISITES**: Always used.
 *
 * @param[in,out] ptr_node Pointer to node structure
 *
 * **RETURN**: The error status.
 *
 * **LONG DESCRIPTION**:
 *
 * Cells in the parent node labeled in the box array \ref node.ptr_box "ptr_box"
 * with the status of "Requires-Refinement" (-1) are organized in different
 * groups which will be the new refinement zones of the parent node. This new
 * zones are created and stored in the pointer to the arrays \ref
 * node.pptr_zones "pptr_zones".
 *
 * Every refinement zone has the same properties as the previous one, i.e. it is
 * built by joining only the 6 closest neighboring cells localized in left, 
 * right, over, under, forward, and backward directions. This representation can 
 * be seen in the figure (work in progress). So, two different zones can not
 * have any of their face cells joined, but they can share edges or corners.
 *
 * Roughly speaking, the way to perform the creation of every zone consists of
 * two steps, *(a)* pick one cell that requires refinement but that is not yet 
 * in a refinement zone, and *(b)* create the zone using this initial cell. The
 * second step is performed by adding every cell to the zone array which is in
 * "face-contact" with another cell of this zone array.
 *
 * The flux of this  function can be seen in the figure (work in progress), and
 * it is explained below:
 *
 * - [0]  <b> THE fill_zones_ref() FUNCTION STARTS....</b>
 *
 * - [1]  Reset parameters if periodic boundary conditions are activated (pbc).
 *
 * - [2]  A "while" loop runs over all cells which require refinement, i.e.
 * cells belonging to the array \ref node.ptr_cell_ref "ptr_cell_ref", that
 * have not been added to any refinement zone yet. A counter is used to
 * perform this task.
 *
 * - [3]  In every step of the "while" loop over the cells, we ask if the cell
 * belongs to a subzone of refinement (box value \f$ \geq 0 \f$). If it belongs,
 * we continue asking to the next cell, if it does not belong (box value \f$ = 
 * -1 \f$), we change its box status to "zone_idx" and pass to the next step.
 *
 * - [4]  Having found a cell with no subzone, we are going to create an 
 * entirely new zone starting with this cell as the foundation stone. To perform 
 * this, a new "while" loop is executed running until there are no more cells 
 * without analysis in the zone, i.e. the final state of the block found is a 
 * solid piece completely isolated from the other zones, labeled with the 
 * "zone_idx" value in the parent node box.
 *
 * - [5]  At this point both "while" loops of the steps [2] and [3] end. Now,
 * using the box with this new information the zones array \ref
 * node.pptr_zones "pptr_zones" is filled with new zones of refinement. To do
 * this, a "for" loop over the refined cells array \ref node.ptr_cell_ref
 * "ptr_cell_ref" is performed.
 *
 * - [6]  Finally, the box status of the parent node is returned to its initial
 * state, putting the value of \e zone_idx in the corresponding cell of the new
 * box.
 *
 * - [7] <b> THE fill_zones_ref() FUNCTION ENDS....</b>
 *
 * **ILLUSTRATIVE EXAMPLES**:
 * - [a]  Trivial.
 *
 * **RATIONALES**:
 * - [a]  This design allows the contact of edges and corks between the 
 * different refinement zones, it also allows particles inside of a child node 
 * to jump to another child node of the same parent node in one time-step of the 
 * simulation. This leads to a computational spend (see particle_updating.c 
 * module) can be avoided  if we decided to change the design of the refinement 
 * zone considering now to ask to more cells than the closest 6 neighboring. 
 * However, this modification also implies the computational spend to ask to 
 * these new neighboring and also can increase the size of the refinement zone 
 * which in turn makes it difficult to load balance in multiple cores.
 *
 * **NOTES**:
 * - [a]
 */

static int fill_zones_ref(struct node *ptr_node)
{

  //* >> Filling the auxiliary diferent refinement zones *//
  int cntr_cell_add_all_zones = 0; // Counter Number of cells added to any refinement zone
  int cntr_cell_add;               // Counter number of cells added per inspection
  int cntr_insp;                   // Numer of inspected cells in the zone

  int zone_size;                          // Number of cells in the zone
  int zone_idx_max = ptr_node->zones_cap; // Maximum id of the zone. It is equal to to the capacity in the number of zones
  int zone_idx = 0;                       // Index of the zone. Initialized at zone 0

  int box_idx; // Box index

  int box_idxNbr_x_plus;  // Box index in the neigborhood on the right
  int box_idxNbr_x_minus; // Box index in the neigborhood on the left
  int box_idxNbr_y_plus;  // Box index in the neigborhood behind
  int box_idxNbr_y_minus; // Box index in the neigborhood in front
  int box_idxNbr_z_plus;  // Box index in the neigborhood up
  int box_idxNbr_z_minus; // Box index in the neigborhood down

  int cell_idx;         // The cell index is simply i of the for loop
  int cell_ref_idx = 0; // The index in the cell refined array ptr_cell_ref

  int lv = ptr_node->lv;

  int box_real_dim_X_node = ptr_node->box_real_dim_x;
  int box_real_dim_X_times_Y_node = ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;

  // The auxiliary array ptr_aux_idx will be used to store the box indices of the zone

  //* >> Space checking of auxiliary array ptr_aux_idx *//
  // The size will always be bigger or equal than the number of refined cells

  if (space_check(&(ptr_node->aux_idx_cap), ptr_node->cell_ref_cap, 1.0f, "p1i1", &(ptr_node->ptr_aux_idx)) == _FAILURE_)
  {
    printf("Error, in space_check function\n");
    return _FAILURE_;
  }

  if (bdry_cond_type == 0 && ptr_node->pbc_crosses_whole_sim_box == true)
  {

    // Reseting status of anomalies
    for (int i = 0; i < ptr_node->pbc_bool_bdry_anomalies_cap; i++)
    {
      ptr_node->ptr_pbc_bool_bdry_anomalies_x[i] = false;
      ptr_node->ptr_pbc_bool_bdry_anomalies_y[i] = false;
      ptr_node->ptr_pbc_bool_bdry_anomalies_z[i] = false;
    }

    if (space_check(&(ptr_node->pbc_bool_bdry_anomalies_cap), ptr_node->cell_ref_cap, 1.0f, "p3b1b1b1", &(ptr_node->ptr_pbc_bool_bdry_anomalies_x), &(ptr_node->ptr_pbc_bool_bdry_anomalies_y), &(ptr_node->ptr_pbc_bool_bdry_anomalies_z)) == _FAILURE_)
    {
      printf("Error, in space_check function\n");
      return _FAILURE_;
    }
  }

  //* >>  Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0) *//
  while (ptr_node->cell_ref_size > cntr_cell_add_all_zones) // The loop while works as long as the number of cell addeed is less than the total refined cells
  {
    // Notes that the initiality we inspect the elements in the refined cell array until an element has been found that is not found in any of the current refinement zones
    cell_idx = ptr_node->ptr_cell_ref[cell_ref_idx]; // Index of the cells array in the node
    box_idx = ptr_node->ptr_box_idx[cell_idx];

    if (ptr_node->ptr_box[box_idx] == -1) // A cell without zone has been founded
    {
      zone_size = 0; // Initial number of element in the zone

      //* >> Including the first element of the box to the auxiliary array ptr_aux_idx *//
      ptr_node->ptr_aux_idx[0] = box_idx;

      //* >>  Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0) *//
      ptr_node->ptr_box[box_idx] = zone_idx;

      zone_size++;               // +1 to the number of cells in the zone
      cntr_cell_add_all_zones++; // +1 to the number of cells added in total

      //* >> Building of the zone of ID = zone_idx *//
      cntr_insp = 0;                                                                     // Counter for the number of elements inspectioned in the current zone array
      while (zone_size > cntr_insp && ptr_node->cell_ref_size > cntr_cell_add_all_zones) // The loop ends when all the elements of the zone have been inspected or when the number of cells refined is equal to the total number of cells added.
      {
        // Note that the number of elements in the zone increases if the neighbors of the inspected cell must be added to the zone

        cntr_cell_add = 0; // Counter number of cells added per cell inspected
        box_idx = ptr_node->ptr_aux_idx[cntr_insp];

        box_idxNbr_x_plus = box_idx + 1;
        box_idxNbr_x_minus = box_idx - 1;
        box_idxNbr_y_plus = box_idx + box_real_dim_X_node;
        box_idxNbr_y_minus = box_idx - box_real_dim_X_node;
        box_idxNbr_z_plus = box_idx + box_real_dim_X_times_Y_node;
        box_idxNbr_z_minus = box_idx - box_real_dim_X_times_Y_node;

        if (bdry_cond_type == 0 && ptr_node->pbc_crosses_whole_sim_box == true)
        {
          if (ptr_node->pbc_crosses_whole_sim_box_x == true)
          {
            if (ptr_node->ptr_box[box_idxNbr_x_plus] == -6)
            {
              box_idxNbr_x_plus -= (1 << lv);
              if (ptr_node->ptr_box[box_idxNbr_x_plus] == -1 || ptr_node->ptr_box[box_idxNbr_x_plus] == zone_idx)
              {
                ptr_node->ptr_pbc_bool_bdry_anomalies_x[zone_idx] = true;
              }
            }
            else if (ptr_node->ptr_box[box_idxNbr_x_minus] == -6)
            {
              box_idxNbr_x_minus += (1 << lv);
              if (ptr_node->ptr_box[box_idxNbr_x_minus] == -1 || ptr_node->ptr_box[box_idxNbr_x_minus] == zone_idx)
              {
                ptr_node->ptr_pbc_bool_bdry_anomalies_x[zone_idx] = true;
              }
            }
          }
          if (ptr_node->pbc_crosses_whole_sim_box_y == true)
          {
            if (ptr_node->ptr_box[box_idxNbr_y_plus] == -6)
            {
              box_idxNbr_y_plus -= (1 << lv) * box_real_dim_X_node;
              if (ptr_node->ptr_box[box_idxNbr_y_plus] == -1 || ptr_node->ptr_box[box_idxNbr_y_plus] == zone_idx)
              {
                ptr_node->ptr_pbc_bool_bdry_anomalies_y[zone_idx] = true;
              }
            }
            else if (ptr_node->ptr_box[box_idxNbr_y_minus] == -6)
            {
              box_idxNbr_y_minus += (1 << lv) * box_real_dim_X_node;
              if (ptr_node->ptr_box[box_idxNbr_y_minus] == -1 || ptr_node->ptr_box[box_idxNbr_y_minus] == zone_idx)
              {
                ptr_node->ptr_pbc_bool_bdry_anomalies_y[zone_idx] = true;
              }
            }
          }

          if (ptr_node->pbc_crosses_whole_sim_box_z == true)
          {
            if (ptr_node->ptr_box[box_idxNbr_z_plus] == -6)
            {
              box_idxNbr_z_plus -= (1 << lv) * box_real_dim_X_times_Y_node;
              // printf("\nbox_idxNbr_z_plus == -6\n");
              // printf("ptr_node->ptr_box[box_idxNbr_z_plus] = %d\n", ptr_node->ptr_box[box_idxNbr_z_plus]);
              if (ptr_node->ptr_box[box_idxNbr_z_plus] == -1 || ptr_node->ptr_box[box_idxNbr_z_plus] == zone_idx)
              {
                ptr_node->ptr_pbc_bool_bdry_anomalies_z[zone_idx] = true;
              }
            }
            else if (ptr_node->ptr_box[box_idxNbr_z_minus] == -6)
            {
              box_idxNbr_z_minus += (1 << lv) * box_real_dim_X_times_Y_node;
              // printf("\nbox_idxNbr_z_minus == -6\n");
              // printf("ptr_node->ptr_box[box_idxNbr_z_minus] = %d\n", ptr_node->ptr_box[box_idxNbr_z_minus]);
              if (ptr_node->ptr_box[box_idxNbr_z_minus] == -1 || ptr_node->ptr_box[box_idxNbr_z_minus] == zone_idx)
              {
                ptr_node->ptr_pbc_bool_bdry_anomalies_z[zone_idx] = true;
              }
            }
          }
        }

        //* Checking the nearest 6 neighbors of face
        // First neighbor
        if (ptr_node->ptr_box[box_idxNbr_x_plus] == -1)
        {
          ptr_node->ptr_aux_idx[zone_size + cntr_cell_add] = box_idxNbr_x_plus; // Including the neighboring element of the box to the auxiliary array
          ptr_node->ptr_box[box_idxNbr_x_plus] = zone_idx;                      // Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0)
          cntr_cell_add++;                                                      // +1 to the number of cells added in the current inspection
        }
        // Second neighbor
        if (ptr_node->ptr_box[box_idxNbr_x_minus] == -1)
        {
          ptr_node->ptr_aux_idx[zone_size + cntr_cell_add] = box_idxNbr_x_minus; // Including the neighboring element of the box to the auxiliary array
          ptr_node->ptr_box[box_idxNbr_x_minus] = zone_idx;                      // Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0)
          cntr_cell_add++;                                                       // +1 to the number of cells added in the current inspection
        }
        // Third neighbor
        if (ptr_node->ptr_box[box_idxNbr_y_plus] == -1)
        {
          ptr_node->ptr_aux_idx[zone_size + cntr_cell_add] = box_idxNbr_y_plus; // Including the neighboring element of the box to the auxiliary array
          ptr_node->ptr_box[box_idxNbr_y_plus] = zone_idx;                      // Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0)
          cntr_cell_add++;                                                      // +1 to the number of cells added in the current inspection
        }
        // Fourth neighbor
        if (ptr_node->ptr_box[box_idxNbr_y_minus] == -1)
        {
          ptr_node->ptr_aux_idx[zone_size + cntr_cell_add] = box_idxNbr_y_minus; // Including the neighboring element of the box to the auxiliary array
          ptr_node->ptr_box[box_idxNbr_y_minus] = zone_idx;                      // Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0)
          cntr_cell_add++;                                                       // +1 to the number of cells added in the current inspection
        }
        // Fifth neighbor
        if (ptr_node->ptr_box[box_idxNbr_z_plus] == -1)
        {
          ptr_node->ptr_aux_idx[zone_size + cntr_cell_add] = box_idxNbr_z_plus; // Including the neighboring element of the box to the auxiliary array
          ptr_node->ptr_box[box_idxNbr_z_plus] = zone_idx;                      // Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0)
          cntr_cell_add++;                                                      // +1 to the number of cells added in the current inspection
        }
        // Sixth neighbor
        if (ptr_node->ptr_box[box_idxNbr_z_minus] == -1)
        {
          ptr_node->ptr_aux_idx[zone_size + cntr_cell_add] = box_idxNbr_z_minus; // Including the neighboring element of the box to the auxiliary array
          ptr_node->ptr_box[box_idxNbr_z_minus] = zone_idx;                      // Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0)
          cntr_cell_add++;                                                       // +1 to the number of cells added in the current inspection
        }
        zone_size += cntr_cell_add;               // Increasing the number of cells in the zone
        cntr_cell_add_all_zones += cntr_cell_add; // Increasing the number of cells added to the zones
        cntr_insp++;                              // Increasing the number of inspections
      }                                           // End "while" loop, now the box contains the the information about all cells of the zone "zone_idx"

      //* >> Space checking of refinement zones arrays, refinement capacity array and refinement size array *//
      if (zone_idx_max < zone_idx + 1)
      {
        if (space_check(&(zone_idx_max), zone_idx + 1, 2.0f, "p3i2i1i1", &(ptr_node->pptr_zones), &(ptr_node->ptr_zone_cap), &(ptr_node->ptr_zone_size)) == _FAILURE_)
        {
          printf("Error, in space_check function\n");
          return _FAILURE_;
        }
      }
      ptr_node->ptr_zone_size[zone_idx] = zone_size;
      zone_idx++; // Increasing the zone number
    }             // Zone defined in the box
    cell_ref_idx++;
  } // At this point the box_aux contains the information of all new refinement zones

  ptr_node->zones_cap = zone_idx_max; // Maximum amount of zones
  ptr_node->zones_size = zone_idx;    // Total amount of zones

  //* >> Space checking in each zone *//
  for (int i = 0; i < ptr_node->zones_size; i++)
  {
    if (space_check(&(ptr_node->ptr_zone_cap[i]), ptr_node->ptr_zone_size[i], 2.0f, "p1i1", &(ptr_node->pptr_zones[i])) == _FAILURE_)
    {
      printf("Error, in space_check function\n");
      return _FAILURE_;
    }
  }
  //* >> Initializing the ptr_aux_idx array to be used as a counter of elemnents in each zone*//
  // Notes that the number of refined cells is always bigger or equal than the number of refined zones, so that the capacity of ptr_aux_idx is always enough to counter the number of refined zones
  for (int i = 0; i < ptr_node->zones_size; i++)
  {
    ptr_node->ptr_aux_idx[i] = 0;
  }
  //* >> Adding the cells to the zone array pptr_zones *//
  for (int i = 0; i < ptr_node->cell_ref_size; i++)
  {
    cell_idx = ptr_node->ptr_cell_ref[i];
    box_idx = ptr_node->ptr_box_idx[cell_idx];
    zone_idx = ptr_node->ptr_box[box_idx];
    cntr_cell_add = ptr_node->ptr_aux_idx[zone_idx];          // Counter the element in the zone "zone_idx"
    ptr_node->pptr_zones[zone_idx][cntr_cell_add] = cell_idx; // Adding the index of the cell array in the block to the zone
    ptr_node->ptr_aux_idx[zone_idx] += 1;                     // Counter the number of elements added in the zone "zone_idx"
  }
  return _SUCCESS_;
} // End function fill_zones_ref

/**
 * @brief Create the relation between old zones and new zones of refinement in
 * the node.
 *
 * **SHORT DESCRIPTION**: Create the relation between old zones and new zones of
 * refinement in the node. There are 4 integer arrays of size equal to the number
 * of new zones of refinement. Every array contains the identifications of the
 * refinement zones either old or new ones, and they are sorted in different 
 * ways.
 *
 * **PREREQUISITES**: If there is at least a new refinement zone.
 *
 * @param[in,out] ptr_node Pointer to node structure
 *
 * **RETURN**: The error status.
 *
 * **LONG DESCRIPTION**:
 *
 * Create the relation between old zones and new zones of refinement in the 
 * node. To do this, 4 integer arrays are created, all of them with a size equal 
 * to the number of refinement zones. Each of these arrays is ordered in a 
 * particular way such that its nth element is connected to the other nth 
 * elements of the other arrays.
 *
 * The first array \ref node.ptr_links_old_ord_old "ptr_links_old_ord_old", is
 * organized in increasing order, and it contains the old refinement zones which
 * will be reused and if necessary the new nodes that are going to be created 
 * because the number of new zones is bigger than the number of old zones of
 * refinement.
 *
 * The second array \ref node.ptr_links_new_ord_new "ptr_links_new_ord_new", is
 * organized in increasing order, and it contains the new refinement zones. So
 * it always contains the consecutive integers from 0 to new_zone_max_id.
 *
 * The next 2 arrays are in charge of connecting the old zones with the new 
 * zones and vice versa.
 *
 * The third array \ref node.ptr_links_new_ord_old "ptr_links_new_ord_old", is
 * organized using the sort of the old sort array \ref 
 * node.ptr_links_old_ord_old "ptr_links_old_ord_old", and it contains the new 
 * refinement zones. This means that for example the nth element of the new zone 
 * array \ref node.ptr_links_new_ord_old "ptr_links_new_ord_old" is connected 
 * with the nth element of the old zone array \ref node.ptr_links_old_ord_old 
 * "ptr_links_old_ord_old".
 *
 * The fourth array \ref node.ptr_links_old_ord_new "ptr_links_old_ord_new" is 
 * the inverse of the third array, it is organized using the sort of the new 
 * sort array \ref node.ptr_links_new_ord_new "ptr_links_new_ord_new", and it 
 * contains the old refinement zones. This means that for example the nth 
 * element of the old zone array  \ref node.ptr_links_old_ord_new 
 * "ptr_links_old_ord_new" is connected with the nth element of the new zone 
 * array  \ref node.ptr_links_new_ord_new "ptr_links_new_ord_new".
 *
 * How the function associates an old refinement zone with a new one is
 * considering the maximum number of matching cells between the old zones and
 * the new zones.
 *
 * The flux of this function can be seen in the figure (work in progress), and
 * it is explained below:
 *
 * - [0]  <b> THE create_links() FUNCTION STARTS....</b>
 *
 * - [1]  Defining some internal useful parameters.
 *
 * - [2]  The first part of the function is in charge of associating the old 
 * node zones with the new ones using considering the maximum number of matching 
 * cells between the old zones and the new zones as a criterion of link.
 *
 * - [3]  Run a "for" loop over the child nodes, which corresponds to the 
 * number of old zones of refinement.
 *
 * - [4]  Run a "for" loop over the cells in every child node. The loop runs
 * each 8 cells which correspond to one cell in the parent node.
 *
 * - [5]  The parent box index associated with this package of 8 child cells is
 * computed. Moreover, for each child node, the number of cells matching between 
 * the old zone and the new zones is computed.
 *
 * - [6]  Both "for" loops over the cells and child nodes end.
 *
 * - [7]  For each child, the largest number of cells that match a single new 
 * zone is computed. The match with the new zones can be repeated, so for 
 * example, the old zone 0, and 1, could have their largest number equal to 50
 * and 80 respectively, but both old zones could match with the single new zone 
 * 3.
 *
 * - [8]  As the relation between the old zones and the new ones needs to be
 * injective, it is necessary to choose what old zones should link with a single
 * new zone without repetitions. To perform that, the old zone with the largest
 * number of matches between all the old zones is linked first with the 
 * corresponding new zone. If that this new zone corresponds to the largest 
 * number of matches of other old zones, their largest number is recomputed but 
 * discarding the new zone as an option. So, at the end of this process, every 
 * old zone will have at most a unique new zone of refinement.
 *
 * - [9]  Using array with the highest number of repetitions, and the 
 * corresponding new zones of matching, the link arrays with old order \ref 
 * node.ptr_links_old_ord_old "ptr_links_old_ord_old", and \ref 
 * node.ptr_links_new_ord_old "ptr_links_new_ord_old" are filled.
 *
 * - [10] If the number of new zones of refinement is bigger than the old number
 * of zones, the arrays with old order \ref node.ptr_links_old_ord_old 
 * "ptr_links_old_ord_old", and \ref node.ptr_links_new_ord_old 
 * "ptr_links_new_ord_old" are filled with the remaining new zones and with the 
 * new nodes identification to be created.
 *
 * - [11] Now, the old order arrays are completely filled but they still require 
 * to be sorted by the old order scheme. Thus, they are reorganized in such a 
 * way.
 *
 * - [12] Finally, using the old order arrays, the new order arrays are 
 * completely filled.
 *
 * - [13]  <b> THE create_links() FUNCTION ENDS....</b>
 *
 * **ILLUSTRATIVE EXAMPLES**:
 * - [a]  Trivial.
 *
 * **RATIONALES**:
 * - [a]  The criteria used to decide what old zone should be linked with a new
 * one, was initially simplified to only ask for the first cell which matched 
 * with a new zone. However, this criterion, although with a faster linking 
 * process, led to that on several occasions old zones were linked with new 
 * zones only by a few amounts of cells existing other old zones with a much 
 * higher number of matches for that new zones. So, in the end, the time 
 * consuming in to adapt these old zones to fit with this wrong new zone 
 * requires a lot of computation time.
 *
 * **NOTES**:
 * - [a]  a
 */

static int create_links(struct node *ptr_node)
{
  struct node *ptr_ch; // child node

  int cntr_links_plus;    // Counte links
  int cntr_link_elements; // Counter linked elements
  bool check_link;        // Check if the element is linked

  int box_value_new; // Value in the new box

  int box_idx_x_node; // Box index at X direction
  int box_idx_y_node; // Box index at Y direction
  int box_idx_z_node; // Box index at Z direction
  int box_idx_node;   // Box index

  int aux_idx_x;
  int aux_idx_y;
  int aux_idx_z;

  int aux_min;
  int element_idx;
  int aux_int;

  int lv = ptr_node->lv;

  int box_real_dim_X_node = ptr_node->box_real_dim_x;
  int box_real_dim_X_times_Y_node = ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;

  //* >> Creating the links between old and new refinement zones IDs *//

  int cntr_links = 0; // Counter the number of links between the new and old zones of refinement

  //* >> Space checking of the capacity of links order arrays *//
  if (space_check(&(ptr_node->links_cap), ptr_node->zones_size, 2.0f, "p4i1i1i1i1", &(ptr_node->ptr_links_old_ord_old), &(ptr_node->ptr_links_new_ord_old), &(ptr_node->ptr_links_old_ord_new), &(ptr_node->ptr_links_new_ord_new)) == _FAILURE_)
  {
    printf("Error, in space_check function\n");
    return _FAILURE_;
  }

  if (ptr_node->chn_size > 0)
  {
    int **pptr_cntr_zones;
    pptr_cntr_zones = (int **)malloc(ptr_node->chn_size * sizeof(int *));
    for (int i = 0; i < ptr_node->chn_size; i++)
    {
      pptr_cntr_zones[i] = (int *)calloc(ptr_node->zones_size, sizeof(int));
    }

    // Filling the repetitions zone number of each child
    for (int ch = 0; ch < ptr_node->chn_size; ch++)
    {
      ptr_ch = ptr_node->pptr_chn[ch];
      for (int cell_idx = 0; cell_idx < ptr_ch->cell_size; cell_idx += 8)
      {

        aux_idx_x = ptr_ch->ptr_cell_idx_x[cell_idx] >> 1;
        aux_idx_y = ptr_ch->ptr_cell_idx_y[cell_idx] >> 1;
        aux_idx_z = ptr_ch->ptr_cell_idx_z[cell_idx] >> 1;

        box_idx_x_node = aux_idx_x - ptr_node->box_ts_x;
        box_idx_y_node = aux_idx_y - ptr_node->box_ts_y;
        box_idx_z_node = aux_idx_z - ptr_node->box_ts_z;

        if (ptr_node->pbc_crosses_sim_box_bdry == true)
        {
          if (aux_idx_x > ptr_node->box_max_x)
          {
            box_idx_x_node -= (1 << lv);
          }

          if (aux_idx_y > ptr_node->box_max_y)
          {
            box_idx_y_node -= (1 << lv);
          }

          if (aux_idx_z > ptr_node->box_max_z)
          {
            box_idx_z_node -= (1 << lv);
          }
        }

        box_idx_node = box_idx_x_node + box_idx_y_node * box_real_dim_X_node + box_idx_z_node * box_real_dim_X_times_Y_node;
        box_value_new = ptr_node->ptr_box[box_idx_node];
        if (box_value_new >= 0)
        {
          pptr_cntr_zones[ch][box_value_new] += 1;
        }
      }
    }

    int *ptr_higher_repetitions, *ptr_zone_higher_repetitions;
    ptr_higher_repetitions = (int *)calloc(ptr_node->chn_size, sizeof(int));
    ptr_zone_higher_repetitions = (int *)calloc(ptr_node->chn_size, sizeof(int));

    // Finding the largest number of cells that match a single new zone

    for (int ch = 0; ch < ptr_node->chn_size; ch++)
    {
      ptr_higher_repetitions[ch] = pptr_cntr_zones[ch][0];

      for (int zone_idx = 1; zone_idx < ptr_node->zones_size; zone_idx++)
      {
        if (ptr_higher_repetitions[ch] < pptr_cntr_zones[ch][zone_idx])
        {
          ptr_higher_repetitions[ch] = pptr_cntr_zones[ch][zone_idx];
          ptr_zone_higher_repetitions[ch] = zone_idx;
        }
      }
    }

    // Reorganization of pptr_cntr_zones according to the biggest number of zone repetitions
    int aux_pos;
    int aux_int;
    int *ptr_ch_ID;
    ptr_ch_ID = (int *)malloc(ptr_node->chn_size * sizeof(int));

    for (int ch = 0; ch < ptr_node->chn_size; ch++)
    {
      ptr_ch_ID[ch] = ch;
    }

    for (int ch = 0; ch < ptr_node->chn_size; ch++)
    {
      aux_pos = ch;
      for (int i = ch + 1; i < ptr_node->chn_size; i++)
      {
        if (ptr_higher_repetitions[ptr_ch_ID[aux_pos]] < ptr_higher_repetitions[ptr_ch_ID[i]])
        {
          aux_pos = i;
        }
      }
      aux_int = ptr_ch_ID[ch];
      ptr_ch_ID[ch] = ptr_ch_ID[aux_pos];
      ptr_ch_ID[aux_pos] = aux_int;

      if (ptr_node->zones_size > ch)
      {
        for (int i = ch + 1; i < ptr_node->chn_size; i++)
        {
          pptr_cntr_zones[ptr_ch_ID[i]][ptr_zone_higher_repetitions[ptr_ch_ID[ch]]] = -1;
          if (ptr_zone_higher_repetitions[ptr_ch_ID[i]] == ptr_zone_higher_repetitions[ptr_ch_ID[ch]])
          {
            ptr_higher_repetitions[ptr_ch_ID[i]] = pptr_cntr_zones[ptr_ch_ID[i]][0];
            ptr_zone_higher_repetitions[ptr_ch_ID[i]] = 0;
            for (int zone_idx = 1; zone_idx < ptr_node->zones_size; zone_idx++)
            {
              if (ptr_higher_repetitions[ptr_ch_ID[i]] < pptr_cntr_zones[ptr_ch_ID[i]][zone_idx])
              {
                ptr_higher_repetitions[ptr_ch_ID[i]] = pptr_cntr_zones[ptr_ch_ID[i]][zone_idx];
                ptr_zone_higher_repetitions[ptr_ch_ID[i]] = zone_idx;
              }
            }
          }
        }
      }
      else
      {
        ch = ptr_node->chn_size;
      }
    }

    // cntr_links = 0 initially
    while (cntr_links < ptr_node->zones_size && cntr_links < ptr_node->chn_size)
    {
      ptr_node->ptr_links_old_ord_old[cntr_links] = ptr_ch_ID[cntr_links];
      ptr_node->ptr_links_new_ord_old[cntr_links] = ptr_zone_higher_repetitions[ptr_ch_ID[cntr_links]];
      cntr_links++;
    }

    for (int i = 0; i < ptr_node->chn_size; i++)
    {
      free(pptr_cntr_zones[i]);
    }
    free(pptr_cntr_zones);
    free(ptr_higher_repetitions);
    free(ptr_zone_higher_repetitions);
    free(ptr_ch_ID);
  }

  if (ptr_node->chn_size < ptr_node->zones_size)
  {
    cntr_link_elements = -1;
    for (int ch = ptr_node->chn_size; ch < ptr_node->zones_size; ch++)
    {
      ptr_node->ptr_links_old_ord_old[ch] = ch;

      check_link = true;
      while (check_link == true)
      {
        check_link = false;
        cntr_link_elements++;
        for (int j = 0; j < cntr_links; j++)
        {
          if (ptr_node->ptr_links_new_ord_old[j] == cntr_link_elements)
          {
            check_link = true;
            j = cntr_links;
          }
        }
      }
      ptr_node->ptr_links_new_ord_old[ch] = cntr_link_elements;
    }
  }

  //* >> Ordering links arrays *//
  cntr_links_plus = 0;
  for (int i = 0; i < ptr_node->zones_size - 1; i++)
  {
    aux_min = ptr_node->chn_size + ptr_node->zones_size;
    element_idx = -1;
    cntr_link_elements = i;
    while (aux_min > cntr_links_plus && cntr_link_elements < ptr_node->zones_size)
    {
      if (aux_min > ptr_node->ptr_links_old_ord_old[cntr_link_elements])
      {
        aux_min = ptr_node->ptr_links_old_ord_old[cntr_link_elements];
        element_idx = cntr_link_elements;
      }
      cntr_link_elements++;
    }

    if (ptr_node->ptr_links_old_ord_old[i] != aux_min)
    {
      //* >> Old *//
      aux_int = ptr_node->ptr_links_old_ord_old[i];
      ptr_node->ptr_links_old_ord_old[i] = aux_min;
      ptr_node->ptr_links_old_ord_old[element_idx] = aux_int;
      //* >> New *//
      aux_int = ptr_node->ptr_links_new_ord_old[i];
      ptr_node->ptr_links_new_ord_old[i] = ptr_node->ptr_links_new_ord_old[element_idx];
      ptr_node->ptr_links_new_ord_old[element_idx] = aux_int;
    }
    cntr_links_plus = aux_min + 1;
  }

  //* >> New order *//
  for (int i = 0; i < ptr_node->zones_size; i++)
  {
    ptr_node->ptr_links_new_ord_new[i] = i;
  }
  for (int i = 0; i < ptr_node->zones_size; i++)
  {
    for (int j = 0; j < ptr_node->zones_size; j++)
    {
      if (ptr_node->ptr_links_new_ord_old[j] == i)
      {
        ptr_node->ptr_links_old_ord_new[i] = ptr_node->ptr_links_old_ord_old[j];
        j = ptr_node->zones_size;
      }
    }
  }
  return _SUCCESS_;
}

/**
 * @brief Removing cells belonging to child nodes that no longer require
 * refinement in the parent node
 *
 * **SHORT DESCRIPTION**: Removing cells belonging to child nodes that no longer
 * require refinement in the parent node, namely, cells that were previously in
 * some old refinement zone but are no longer in any of the new refinement
 * zones. To be more precise, only cells belonging to child nodes that are
 * linked to the new refinement zones are removed.
 *
 * **PREREQUISITES**: If there is at least one new and one old refinement zone.
 *
 * @param[in,out] ptr_node Pointer to node structure
 *
 * **RETURN**: No parameter is returned.
 *
 * **LONG DESCRIPTION**:
 *
 * The following description is more precise than the *Short Description* above:
 *
 * This function is in charge of removing cells from linked child nodes. The 
 * cells removed are cells whose parent cells do not belong to any new 
 * refinement zone. Here, the meaning of "removing" a cell of the chile node 
 * consists in remove the cell and particles from the respective arrays, 
 * modifying its box status, and other local parameters of the child node.
 *
 * The flux of this function can be seen in the figure (work in progress), and
 * it is explained below:
 *
 * - [0]  <b> THE remov_cells_nolonger_require_refinement() FUNCTION STARTS....</b>
 *
 * - [1]  Defining some internal useful parameters.
 *
 * - [2]  Run a "for" loop over the new zones of refinement, but only consider
 * the new zones linked with existing child nodes.
 *
 * - [3]  Run a "for" loop over the cells of the child nodes. The loop runs
 * each 8 cells which correspond to one cell in the parent node.
 *
 * - [4]  The parent box index associated with this package of 8 child cells is
 * computed
 *
 * - [5]  If the box index of the parent node does not belong to any new 
 * refinement zone (value < 0), then the mass and number of particles, 
 * properties of the cell structure, and box status of the child nodes are 
 * updated, and cell are removed from the cell arrays of the child nodes.
 *
 * - [6]  Finnally, the cell *size* (see Key Concepts \ref Key_Concepts_Size "Size") of the
 * child nodes is also updated.
 *
 * - [7] <b> THE remov_cells_nolonger_require_refinement() FUNCTION ENDS....</b>
 *
 * **ILLUSTRATIVE EXAMPLES**:
 * - [a]  Trivial.
 *
 * **RATIONALES**:
 * - [a]  At the end the "for" loop over the new zones of refinement only asks
 * about the old linked zones of refinement. So, it is possible for some old 
 * child nodes to still exist with no null information about their cells but
 * these nodes will always be sent to the memory pool by the function
 * transferring_unused_child_node_to_memory_pool().
 *
 * **NOTES**:
 * - [a]
 */

static void remov_cells_nolonger_require_refinement(struct node *ptr_node)
{

  struct node *ptr_ch;

  int box_idx_ch; // Box index of the child cell

  int aux_int;

  int box_idx_x_node; // Box index in X direcction of the node cell
  int box_idx_y_node; // Box index in Y direcction of the node cell
  int box_idx_z_node; // Box index in Z direcction of the node cell
  int box_idx_node;   // Box index of the node cell

  int aux_idx_x;
  int aux_idx_y;
  int aux_idx_z;

  int lv = ptr_node->lv;

  int no_cells; // Total number of cells in the node

  int box_real_dim_X_node = ptr_node->box_real_dim_x;
  int box_real_dim_X_times_Y_node = ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;

  // Loop over new refinement zones
  for (int zone_idx = 0; zone_idx < ptr_node->zones_size; zone_idx++)
  {
    if (ptr_node->ptr_links_old_ord_old[zone_idx] < ptr_node->chn_size)
    {
      ptr_ch = ptr_node->pptr_chn[ptr_node->ptr_links_old_ord_old[zone_idx]];

      no_cells = ptr_ch->cell_size;
      for (int cell_idx = 0; cell_idx < no_cells; cell_idx += 8)
      {
        aux_idx_x = ptr_ch->ptr_cell_idx_x[cell_idx] >> 1;
        aux_idx_y = ptr_ch->ptr_cell_idx_y[cell_idx] >> 1;
        aux_idx_z = ptr_ch->ptr_cell_idx_z[cell_idx] >> 1;

        box_idx_x_node = aux_idx_x - ptr_node->box_ts_x;
        box_idx_y_node = aux_idx_y - ptr_node->box_ts_y;
        box_idx_z_node = aux_idx_z - ptr_node->box_ts_z;

        if (ptr_node->pbc_crosses_sim_box_bdry == true)
        {
          if (aux_idx_x > ptr_node->box_max_x)
          {
            box_idx_x_node -= (1 << lv);
          }

          if (aux_idx_y > ptr_node->box_max_y)
          {
            box_idx_y_node -= (1 << lv);
          }

          if (aux_idx_z > ptr_node->box_max_z)
          {
            box_idx_z_node -= (1 << lv);
          }
        }

        box_idx_node = box_idx_x_node + box_idx_y_node * box_real_dim_X_node + box_idx_z_node * box_real_dim_X_times_Y_node;
        //* >> The child cell no longer requires refinement *//
        if (ptr_node->ptr_box[box_idx_node] < 0)
        {
          //* >> Removing particles in the cell from the child node *//
          // Updating local mass, cell struct and box

          for (int j = cell_idx; j < cell_idx + 8; j++)
          {
            box_idx_ch = ptr_ch->ptr_box_idx[j];
            // box_idx_ch = (box_idx_x_ch + ii) + (box_idx_y_ch + jj) * ptr_ch->box_real_dim_x + (box_idx_z_ch + kk) * ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;
            ptr_ch->node_mass -= ptr_ch->ptr_cell_struct[box_idx_ch].cell_mass;
            ptr_ch->no_ptcl_full_node -= ptr_ch->ptr_cell_struct[box_idx_ch].ptcl_size;
            ptr_ch->ptr_cell_struct[box_idx_ch].cell_mass = 0;
            ptr_ch->ptr_cell_struct[box_idx_ch].ptcl_size = 0;
            ptr_ch->ptr_box[box_idx_ch] = -4;

            //* >> Removing the cells *//
            aux_int = no_cells - 8 + j - cell_idx;
            ptr_ch->ptr_cell_idx_x[j] = ptr_ch->ptr_cell_idx_x[aux_int];
            ptr_ch->ptr_cell_idx_y[j] = ptr_ch->ptr_cell_idx_y[aux_int];
            ptr_ch->ptr_cell_idx_z[j] = ptr_ch->ptr_cell_idx_z[aux_int];
            ptr_ch->ptr_box_idx[j] = ptr_ch->ptr_box_idx[aux_int];
          }

          no_cells -= 8;
          cell_idx -= 8;
        }
      }
      ptr_ch->cell_size = no_cells;
    }
  }
}

/**
 * @brief Old child nodes are adapted to fit with their new refinement zones
 * associated.
 *
 * **SHORT DESCRIPTION**: The already child nodes corresponding to the old
 * refinement zones linked with new refinement zones are adapted in their
 * parameter to fit with the new zone of refinement.
 *
 * **PREREQUISITES**: If there is at least one new and one old refinement zone.
 *
 * @param[in,out] ptr_node Pointer to node structure
 *
 * **RETURN**: The error status.
 *
 * **LONG DESCRIPTION**:
 *
 * The already child nodes corresponding to the old
 * refinement zones linked with new refinement zones are adapted in their
 * parameter to fit with the new zone of refinement. Some of those parameters
 * are boundary flags, box dimensions and spacial location, translation 
 * constants between cell and box arrays, the density array, and others.
 *
 * There is a very important flag in the node called \ref node.box_check_fit
 * "box_check_fit". Essentially, this flag tells if the node has been modified
 * in its box dimensions. Thus, when an old refinement zone can not fit the new
 * refinement zone linked with it, its box requires to be adapted in its
 * dimensions, and thus the connection between cells and the box of the child 
 * node needs to be updated. To know, if the child node box "fits" the
 * new refinement zone linked, is necessary to compute the dimensions of the new
 * refinement zone and compare them with the current box dimensions.
 *
 * The flux of this function can be seen in the figure (work in progress), and
 * it is explained below:
 *
 * - [0]  <b> THE adapt_child_nodes() FUNCTION STARTS....</b>
 *
 * - [1]  Defining some internal useful parameters.
 *
 * - [2]  Run a "for" loop over the new zones of refinement, but only consider
 * the new zones linked with existing child nodes.
 *
 * - [3]  If periodic boundary conditions (pbc) is activated, and if any boundary
 * flag anomalies \ref node.ptr_pbc_bool_bdry_anomalies_x 
 * "ptr_pbc_bool_bdry_anomalies_x" (\ref node.ptr_pbc_bool_bdry_anomalies_y "y", 
 * \ref node.ptr_pbc_bool_bdry_anomalies_z "z")  is activated, then the minimum 
 * and maximum of the subzones is computed through the function 
 * find_min_max_subzones_ref_PERIODIC_BOUNDARY().
 *
 * - [4]  For every dirección X, Y, and Z, the minimum and maximum, and the
 * boundary flags are computed and stored in local parameters.
 *
 * - [5]  The minimum and maximum of the previous point were computed in the 
 * parent node, so then, they are transformed to the next level of refinement to 
 * be accord with the respective child node.
 *
 * - [6]  These new local parameters obtained above are compared with the old 
 * minimum, maximum, and box dimensions of the child node to then decide if the 
 * old box will fit the new one. For the periodic boundary conditions (pbc), are 
 * necessary to use other local and node parameters to check this point. The 
 * node parameter which will store this fit between the old and the new box is a 
 * boolean one called \ref node.box_check_fit "box_check_fit".
 *
 * - [7]  If \ref node.box_check_fit "box_check_fit" is "true", and the periodic
 * boundary conditions (pbc) are activated, it is necessary to check special 
 * situations with require changes in the traslation constants between the cell 
 * indices and the box indices (see \ref node.box_ts_x "box_ts_x" (\ref 
 * node.box_ts_y "y", \ref node.box_ts_z "z") to know the meaning of the 
 * translation constants.
 *
 * - [8]  In special situations when \ref node.box_check_fit "box_check_fit" is
 * "true" it is necessary to reset the boundary cells box status from "out of 
 * the simulation box" (-5 or -6) to the status of simply "No-Exist" (-4).
 *
 * - [9]  Finally, if the  \ref node.box_check_fit "box_check_fit" is false, it 
 * is necessary to update the new box dimensions parameter of the child node and
 * store a copy of the cell structure which also will require to be updated 
 * through the future functions of the tree_adaptaion.c module.
 *
 * - [10] Finally, the some last parameters of the child node are updated 
 * whatever the status of \ref node.box_check_fit "box_check_fit". These 
 * parameters are size of the grid points \ref node.grid_intr_size 
 * "grid_intr_size", \ref node.grid_bdry_size "grid_bdry_size", and \ref 
 * node.grid_sim_bdry_size "grid_sim_bdry_size", the grid density \ref 
 * node.ptr_d "ptr_d", and the boundary flags. 
 *
 * - [11] <b> THE adapt_child_nodes() FUNCTION ENDS....</b>
 *
 * **ILLUSTRATIVE EXAMPLES**:
 * - [a]  Trivial.
 *
 * **RATIONALES**:
 * - [a]  xxxx
 *
 * **NOTES**:
 * - [a]
 */

static int adapt_child_nodes(struct node *ptr_node)
{
  struct node *ptr_ch;

  // struct cell_struct cell_struct_aux;

  struct cell_struct *ptr_cell_struct_aux;

  // int box_idx_ch_old;   // Box index of the child cell

  int new_box_min_x;
  int new_box_min_y;
  int new_box_min_z;
  int new_box_max_x;
  int new_box_max_y;
  int new_box_max_z;

  bool new_sim_bdry_contact;
  bool new_sim_bdry_contact_x;
  bool new_sim_bdry_contact_y;
  bool new_sim_bdry_contact_z;

  // The following parameters are only used for the Periodic boundary conditions (bdry_cond_type = 0), pbc = periodic boundary conditions
  bool new_pbc_crosses_sim_box_bdry;
  bool new_pbc_crosses_sim_box_bdry_x;
  bool new_pbc_crosses_sim_box_bdry_y;
  bool new_pbc_crosses_sim_box_bdry_z;

  bool new_pbc_crosses_whole_sim_box;
  bool new_pbc_crosses_whole_sim_box_x;
  bool new_pbc_crosses_whole_sim_box_y;
  bool new_pbc_crosses_whole_sim_box_z;

  int cell_idx; // The cell index is simply i of the for loop

  int new_zone_idx; // ID of the new zone of refinement

  int size;

  // int aux_int;

  int aux_fix_min_x;
  int aux_fix_min_y;
  int aux_fix_min_z;
  int aux_fix_max_x;
  int aux_fix_max_y;
  int aux_fix_max_z;

  int box_idx_ch;
  int box_real_dim_X_ch;
  int box_real_dim_X_times_Y_ch;

  int pos_x; // Distance between the real box and the Smallest Box when the last is localized in the middle of the real one
  int pos_y;
  int pos_z;

  bool check; // check subzone analized

  int aux_subzones_analized;

  int lv = ptr_node->lv;

  int bder_box = 1 > n_exp ? 1 : n_exp;

  bool movement_box_x;
  bool movement_box_y;
  bool movement_box_z;

  int aux_box_ts;

  // Loop over new refinement zones
  for (int zone_idx = 0; zone_idx < ptr_node->zones_size; zone_idx++)
  {
    if (ptr_node->ptr_links_old_ord_old[zone_idx] < ptr_node->chn_size)
    {
      new_zone_idx = ptr_node->ptr_links_new_ord_old[zone_idx];
      ptr_ch = ptr_node->pptr_chn[ptr_node->ptr_links_old_ord_old[zone_idx]];

      // Updating the ID of the child nodes
      ptr_ch->ID = new_zone_idx;

      //* CELLS *//
      size = 8 * ptr_node->ptr_zone_size[new_zone_idx] + ptr_ch->cell_size;

      //* >> Space checking of cells indices of the child node
      if (space_check(&(ptr_ch->cell_cap), size, 1.0f, "p4i1i1i1i1", &(ptr_ch->ptr_cell_idx_x), &(ptr_ch->ptr_cell_idx_y), &(ptr_ch->ptr_cell_idx_z), &(ptr_ch->ptr_box_idx)) == _FAILURE_)
      {
        printf("Error, in space_check function\n");
        return _FAILURE_;
      }

      //* BOXES *//
      ptr_ch->box_check_fit = true;

      new_box_min_x = INT_MAX;
      new_box_min_y = INT_MAX;
      new_box_min_z = INT_MAX;
      new_box_max_x = INT_MIN;
      new_box_max_y = INT_MIN;
      new_box_max_z = INT_MIN;

      // Reset  boundary simulation flags
      //  if (ptr_node->sim_bdry_contact == true)
      //  {
      //      //* >> Boundary of the simulation box *//
      //      ptr_ch->sim_bdry_contact = false;
      //      ptr_ch->sim_bdry_contact_x = false;
      //      ptr_ch->sim_bdry_contact_y = false;
      //      ptr_ch->sim_bdry_contact_z = false;

      //     if (ptr_node->pbc_crosses_sim_box_bdry == true)
      //     {
      //         ptr_ch->pbc_crosses_sim_box_bdry = false;
      //         ptr_ch->pbc_crosses_sim_box_bdry_x = false; // when one node croses the simulation box at X direcction
      //         ptr_ch->pbc_crosses_sim_box_bdry_y = false;
      //         ptr_ch->pbc_crosses_sim_box_bdry_z = false;

      //         if (ptr_node->pbc_crosses_whole_sim_box == true)
      //         {
      //             ptr_ch->pbc_crosses_whole_sim_box = false;
      //             ptr_ch->pbc_crosses_whole_sim_box_x = false;
      //             ptr_ch->pbc_crosses_whole_sim_box_y = false;
      //             ptr_ch->pbc_crosses_whole_sim_box_z = false;
      //         }
      //     }
      // }

      new_sim_bdry_contact = false;
      new_sim_bdry_contact_x = false;
      new_sim_bdry_contact_y = false;
      new_sim_bdry_contact_z = false;

      // The following parameters are only used for the Periodic boundary conditions (bdry_cond_type = 0), pbc = periodic boundary conditions
      new_pbc_crosses_sim_box_bdry = false;
      new_pbc_crosses_sim_box_bdry_x = false;
      new_pbc_crosses_sim_box_bdry_y = false;
      new_pbc_crosses_sim_box_bdry_z = false;

      new_pbc_crosses_whole_sim_box = false;
      new_pbc_crosses_whole_sim_box_x = false;
      new_pbc_crosses_whole_sim_box_y = false;
      new_pbc_crosses_whole_sim_box_z = false;

      ptr_ch->pbc_correction_due_pbc_flag = false;

      movement_box_x = false;
      movement_box_y = false;
      movement_box_z = false;

      // printf("\nAdapt child, lv = %d, ID = %d\n", ptr_ch->lv,ptr_ch->ID);

      //* >> The MIN and MAX of the set containig only the new cells to be refined
      // MIN and MAX cell indices values of the node.
      if (bdry_cond_type == 0 && ptr_node->pbc_crosses_whole_sim_box == true &&
          (ptr_node->ptr_pbc_bool_bdry_anomalies_x[new_zone_idx] == true ||
           ptr_node->ptr_pbc_bool_bdry_anomalies_y[new_zone_idx] == true ||
           ptr_node->ptr_pbc_bool_bdry_anomalies_z[new_zone_idx] == true))
      {

        // printf("there are boundary anomalies\n");
        // printf("ptr_node->ptr_pbc_bool_bdry_anomalies_x[new_zone_idx] = %d\n", ptr_node->ptr_pbc_bool_bdry_anomalies_x[new_zone_idx]);
        // printf("ptr_node->ptr_pbc_bool_bdry_anomalies_y[new_zone_idx] = %d\n", ptr_node->ptr_pbc_bool_bdry_anomalies_y[new_zone_idx]);
        // printf("ptr_node->ptr_pbc_bool_bdry_anomalies_z[new_zone_idx] = %d\n", ptr_node->ptr_pbc_bool_bdry_anomalies_z[new_zone_idx]);

        // printf("Calling min and max subzones functionn\n");

        if (find_min_max_subzones_ref_PERIODIC_BOUNDARY(ptr_node, new_zone_idx) == _FAILURE_)
        {
          printf("Error at function fill_subzones_ref_PERIODIC_BOUNDARY()\n");
          return _FAILURE_;
        }

        // printf("ptr_node->pbc_subzones_size = %d\n", ptr_node->pbc_subzones_size);
      }

      //* >> X axis
      if (ptr_node->pbc_crosses_whole_sim_box_x == true)
      {
        if (ptr_node->ptr_pbc_bool_bdry_anomalies_x[new_zone_idx] == false)
        {
          for (int j = 0; j < ptr_node->ptr_zone_size[new_zone_idx]; j++)
          {
            cell_idx = ptr_node->pptr_zones[new_zone_idx][j]; // Cell index in the parent zone
            if (new_box_min_x > ptr_node->ptr_cell_idx_x[cell_idx])
            {
              new_box_min_x = ptr_node->ptr_cell_idx_x[cell_idx];
            }
            if (new_box_max_x < ptr_node->ptr_cell_idx_x[cell_idx])
            {
              new_box_max_x = ptr_node->ptr_cell_idx_x[cell_idx];
            }
          }
        }
        else if (ptr_node->pbc_subzones_size == 1)
        {
          new_box_min_x = 0;
          new_box_max_x = (1 << lv) - 1;
        }
        else
        {
          // Matching the subzones
          aux_subzones_analized = ptr_node->pbc_subzones_size;
          check = true;
          while (aux_subzones_analized > 1 && check == true)
          {
            check = false;
            for (int i = 0; i < 2; i++)
            {
              for (int j = i + 1; j < aux_subzones_analized; j++)
              {
                if ((ptr_node->ptr_pbc_max_subzones_x[j] + 1 >= ptr_node->ptr_pbc_min_subzones_x[i] && ptr_node->ptr_pbc_min_subzones_x[j] <= 1 + ptr_node->ptr_pbc_max_subzones_x[i]))
                {
                  ptr_node->ptr_pbc_min_subzones_x[i] = ptr_node->ptr_pbc_min_subzones_x[i] < ptr_node->ptr_pbc_min_subzones_x[j] ? ptr_node->ptr_pbc_min_subzones_x[i] : ptr_node->ptr_pbc_min_subzones_x[j];
                  ptr_node->ptr_pbc_max_subzones_x[i] = ptr_node->ptr_pbc_max_subzones_x[i] > ptr_node->ptr_pbc_max_subzones_x[j] ? ptr_node->ptr_pbc_max_subzones_x[i] : ptr_node->ptr_pbc_max_subzones_x[j];

                  // Reducing the total amount of analized subzones
                  ptr_node->ptr_pbc_min_subzones_x[j] = ptr_node->ptr_pbc_min_subzones_x[aux_subzones_analized - 1];
                  ptr_node->ptr_pbc_max_subzones_x[j] = ptr_node->ptr_pbc_max_subzones_x[aux_subzones_analized - 1];
                  aux_subzones_analized--;
                  j--;
                  check = true;
                }
              }
            }
          }

          if (aux_subzones_analized > 2 || aux_subzones_analized < 1)
          {
            printf("error, aux_subzones analized is bigger than 2 or lower than 1, it is equal to = %d\n", aux_subzones_analized);
          }
          // printf("aux_subzones analized = %d\n", aux_subzones_analized);

          // min and max between all subzones
          if (aux_subzones_analized == 1)
          {
            // new_box_min_x = ptr_node->ptr_pbc_min_subzones[0];
            // new_box_max_x = ptr_node->ptr_pbc_max_subzones[0];
            new_box_min_x = 0;
            new_box_max_x = (1 << lv) - 1;
          }
          else
          {
            new_box_min_x = ptr_node->ptr_pbc_min_subzones_x[0] == 0 ? ptr_node->ptr_pbc_min_subzones_x[1] - (1 << lv) : ptr_node->ptr_pbc_min_subzones_x[0] - (1 << lv);
            new_box_max_x = ptr_node->ptr_pbc_min_subzones_x[0] == 0 ? ptr_node->ptr_pbc_max_subzones_x[0] : ptr_node->ptr_pbc_max_subzones_x[1];
          }
        }

        // Analysis of min and max of the refinement zone
        // Case refinement zone crosses the whole simulation box
        if (new_box_max_x - new_box_min_x == (1 << lv) - 1)
        {
          new_sim_bdry_contact = true;
          new_sim_bdry_contact_x = true;

          new_pbc_crosses_sim_box_bdry = true;
          new_pbc_crosses_sim_box_bdry_x = true;

          new_pbc_crosses_whole_sim_box = true;
          new_pbc_crosses_whole_sim_box_x = true;

          // ptr_ch->sim_bdry_contact = true;
          // ptr_ch->sim_bdry_contact_x = true;

          // ptr_ch->pbc_crosses_sim_box_bdry = true;
          // ptr_ch->pbc_crosses_sim_box_bdry_x = true;

          // ptr_ch->pbc_crosses_whole_sim_box = true;
          // ptr_ch->pbc_crosses_whole_sim_box_x = true;
        }
        // Case refinement zone crosses the simulation boundary
        else if (new_box_min_x < 0)
        {
          new_sim_bdry_contact = true;
          new_sim_bdry_contact_x = true;

          new_pbc_crosses_sim_box_bdry = true;
          new_pbc_crosses_sim_box_bdry_x = true;

          // ptr_ch->sim_bdry_contact = true;
          // ptr_ch->sim_bdry_contact_x = true;

          // ptr_ch->pbc_crosses_sim_box_bdry = true;
          // ptr_ch->pbc_crosses_sim_box_bdry_x = true;
        }
        // Case refinement zone touches the simulation boundary
        else if (new_box_min_x == 0 || new_box_max_x == (1 << lv) - 1)
        {
          new_sim_bdry_contact = true;
          new_sim_bdry_contact_x = true;

          // ptr_ch->sim_bdry_contact = true;
          // ptr_ch->sim_bdry_contact_x = true;
        }

      } // End if(ptr_node->pbc_crosses_whole_sim_box_x == true)
      // Case Parent x axis crosses the simulation box
      else if (ptr_node->pbc_crosses_sim_box_bdry_x == true)
      {
        for (int j = 0; j < ptr_node->ptr_zone_size[new_zone_idx]; j++)
        {
          cell_idx = ptr_node->pptr_zones[new_zone_idx][j]; // Cell index in the parent zone

          if (ptr_node->ptr_cell_idx_x[cell_idx] > ptr_node->box_max_x)
          {
            new_box_min_x = new_box_min_x < ptr_node->ptr_cell_idx_x[cell_idx] - (1 << lv) ? new_box_min_x : ptr_node->ptr_cell_idx_x[cell_idx] - (1 << lv);
            new_box_max_x = new_box_max_x > ptr_node->ptr_cell_idx_x[cell_idx] - (1 << lv) ? new_box_max_x : ptr_node->ptr_cell_idx_x[cell_idx] - (1 << lv);
          }
          else
          {
            new_box_min_x = new_box_min_x < ptr_node->ptr_cell_idx_x[cell_idx] ? new_box_min_x : ptr_node->ptr_cell_idx_x[cell_idx];
            new_box_max_x = new_box_max_x > ptr_node->ptr_cell_idx_x[cell_idx] ? new_box_max_x : ptr_node->ptr_cell_idx_x[cell_idx];
          }
        }

        // Analysis of min and max of the refinement zone
        if (new_box_min_x < 0 && new_box_max_x >= 0)
        {
          new_sim_bdry_contact = true;
          new_sim_bdry_contact_x = true;

          new_pbc_crosses_sim_box_bdry = true;
          new_pbc_crosses_sim_box_bdry_x = true;

          // ptr_ch->sim_bdry_contact = true;
          // ptr_ch->sim_bdry_contact_x = true;

          // ptr_ch->pbc_crosses_sim_box_bdry = true;
          // ptr_ch->pbc_crosses_sim_box_bdry_x = true;
        }
        else if (new_box_min_x == 0)
        {
          new_sim_bdry_contact = true;
          new_sim_bdry_contact_x = true;

          // ptr_ch->sim_bdry_contact = true;
          // ptr_ch->sim_bdry_contact_x = true;
        }
        else if (new_box_max_x <= -1)
        {
          if (new_box_max_x == -1)
          {
            new_sim_bdry_contact = true;
            new_sim_bdry_contact_x = true;

            // ptr_ch->sim_bdry_contact = true;
            // ptr_ch->sim_bdry_contact_x = true;
          }

          new_box_min_x += (1 << lv);
          new_box_max_x += (1 << lv);
        }
      } // End else if(ptr_node->pbc_crosses_sim_box_bdry = true)
      // Case Parent x axis do not cross the simulation box
      else
      {
        for (int j = 0; j < ptr_node->ptr_zone_size[new_zone_idx]; j++)
        {
          cell_idx = ptr_node->pptr_zones[new_zone_idx][j]; // Cell index in the parent zone
          if (new_box_min_x > ptr_node->ptr_cell_idx_x[cell_idx])
          {
            new_box_min_x = ptr_node->ptr_cell_idx_x[cell_idx];
          }
          if (new_box_max_x < ptr_node->ptr_cell_idx_x[cell_idx])
          {
            new_box_max_x = ptr_node->ptr_cell_idx_x[cell_idx];
          }
        }

        // Analysis of min and max of the refinement zone
        if (new_box_min_x == 0 || new_box_max_x == (1 << lv) - 1)
        {
          new_sim_bdry_contact = true;
          new_sim_bdry_contact_x = true;

          // ptr_ch->sim_bdry_contact = true;
          // ptr_ch->sim_bdry_contact_x = true;
        }
      }

      //* >> Y axis
      if (ptr_node->pbc_crosses_whole_sim_box_y == true)
      {
        if (ptr_node->ptr_pbc_bool_bdry_anomalies_y[new_zone_idx] == false)
        {
          for (int j = 0; j < ptr_node->ptr_zone_size[new_zone_idx]; j++)
          {
            cell_idx = ptr_node->pptr_zones[new_zone_idx][j]; // Cell index in the parent zone
            if (new_box_min_y > ptr_node->ptr_cell_idx_y[cell_idx])
            {
              new_box_min_y = ptr_node->ptr_cell_idx_y[cell_idx];
            }
            if (new_box_max_y < ptr_node->ptr_cell_idx_y[cell_idx])
            {
              new_box_max_y = ptr_node->ptr_cell_idx_y[cell_idx];
            }
          }
        }
        else if (ptr_node->pbc_subzones_size == 1)
        {
          new_box_min_y = 0;
          new_box_max_y = (1 << lv) - 1;
        }
        else
        {
          // Matching the subzones
          aux_subzones_analized = ptr_node->pbc_subzones_size;
          check = true;
          while (aux_subzones_analized > 1 && check == true)
          {
            check = false;
            for (int i = 0; i < 2; i++)
            {
              for (int j = i + 1; j < aux_subzones_analized; j++)
              {
                if ((ptr_node->ptr_pbc_max_subzones_y[j] + 1 >= ptr_node->ptr_pbc_min_subzones_y[i] && ptr_node->ptr_pbc_min_subzones_y[j] <= 1 + ptr_node->ptr_pbc_max_subzones_y[i]))
                {
                  ptr_node->ptr_pbc_min_subzones_y[i] = ptr_node->ptr_pbc_min_subzones_y[i] < ptr_node->ptr_pbc_min_subzones_y[j] ? ptr_node->ptr_pbc_min_subzones_y[i] : ptr_node->ptr_pbc_min_subzones_y[j];
                  ptr_node->ptr_pbc_max_subzones_y[i] = ptr_node->ptr_pbc_max_subzones_y[i] > ptr_node->ptr_pbc_max_subzones_y[j] ? ptr_node->ptr_pbc_max_subzones_y[i] : ptr_node->ptr_pbc_max_subzones_y[j];

                  // Reducing the total amount of analized subzones
                  ptr_node->ptr_pbc_min_subzones_y[j] = ptr_node->ptr_pbc_min_subzones_y[aux_subzones_analized - 1];
                  ptr_node->ptr_pbc_max_subzones_y[j] = ptr_node->ptr_pbc_max_subzones_y[aux_subzones_analized - 1];
                  aux_subzones_analized--;
                  j--;
                  check = true;
                }
              }
            }
          }

          // min and max between all subzones
          if (aux_subzones_analized == 1)
          {
            new_box_min_y = 0;
            new_box_max_y = (1 << lv) - 1;
          }
          else
          {
            new_box_min_y = ptr_node->ptr_pbc_min_subzones_y[0] == 0 ? ptr_node->ptr_pbc_min_subzones_y[1] - (1 << lv) : ptr_node->ptr_pbc_min_subzones_y[0] - (1 << lv);
            new_box_max_y = ptr_node->ptr_pbc_min_subzones_y[0] == 0 ? ptr_node->ptr_pbc_max_subzones_y[0] : ptr_node->ptr_pbc_max_subzones_y[1];
          }
        }

        // Analysis of min and max of the refinement zone
        // Case refinement zone crosses the whole simulation box
        if (new_box_max_y - new_box_min_y == (1 << lv) - 1)
        {
          new_sim_bdry_contact = true;
          new_sim_bdry_contact_y = true;

          new_pbc_crosses_sim_box_bdry = true;
          new_pbc_crosses_sim_box_bdry_y = true;

          new_pbc_crosses_whole_sim_box = true;
          new_pbc_crosses_whole_sim_box_y = true;

          // ptr_ch->sim_bdry_contact = true;
          // ptr_ch->sim_bdry_contact_y = true;

          // ptr_ch->pbc_crosses_sim_box_bdry = true;
          // ptr_ch->pbc_crosses_sim_box_bdry_y = true;

          // ptr_ch->pbc_crosses_whole_sim_box = true;
          // ptr_ch->pbc_crosses_whole_sim_box_y = true;
        }
        // Case refinement zone crosses the simulation boundary
        else if (new_box_min_y < 0)
        {
          new_sim_bdry_contact = true;
          new_sim_bdry_contact_y = true;

          new_pbc_crosses_sim_box_bdry = true;
          new_pbc_crosses_sim_box_bdry_y = true;

          // ptr_ch->sim_bdry_contact = true;
          // ptr_ch->sim_bdry_contact_y = true;

          // ptr_ch->pbc_crosses_sim_box_bdry = true;
          // ptr_ch->pbc_crosses_sim_box_bdry_y = true;
        }
        // Case refinement zone touches the simulation boundary
        else if (new_box_min_y == 0 || new_box_max_y == (1 << lv) - 1)
        {
          new_sim_bdry_contact = true;
          new_sim_bdry_contact_y = true;

          // ptr_ch->sim_bdry_contact = true;
          // ptr_ch->sim_bdry_contact_y = true;
        }

      } // End if(ptr_node->pbc_crosses_whole_sim_box_y == true)
      // Case Parent y axis crosses the simulation box
      else if (ptr_node->pbc_crosses_sim_box_bdry_y == true)
      {
        for (int j = 0; j < ptr_node->ptr_zone_size[new_zone_idx]; j++)
        {
          cell_idx = ptr_node->pptr_zones[new_zone_idx][j]; // Cell index in the parent zone

          if (ptr_node->ptr_cell_idx_y[cell_idx] > ptr_node->box_max_y)
          {
            new_box_min_y = new_box_min_y < ptr_node->ptr_cell_idx_y[cell_idx] - (1 << lv) ? new_box_min_y : ptr_node->ptr_cell_idx_y[cell_idx] - (1 << lv);
            new_box_max_y = new_box_max_y > ptr_node->ptr_cell_idx_y[cell_idx] - (1 << lv) ? new_box_max_y : ptr_node->ptr_cell_idx_y[cell_idx] - (1 << lv);
          }
          else
          {
            new_box_min_y = new_box_min_y < ptr_node->ptr_cell_idx_y[cell_idx] ? new_box_min_y : ptr_node->ptr_cell_idx_y[cell_idx];
            new_box_max_y = new_box_max_y > ptr_node->ptr_cell_idx_y[cell_idx] ? new_box_max_y : ptr_node->ptr_cell_idx_y[cell_idx];
          }
        }

        // Analysis of min and max of the refinement zone
        if (new_box_min_y < 0 && new_box_max_y >= 0)
        {
          new_sim_bdry_contact = true;
          new_sim_bdry_contact_y = true;

          new_pbc_crosses_sim_box_bdry = true;
          new_pbc_crosses_sim_box_bdry_y = true;

          // ptr_ch->sim_bdry_contact = true;
          // ptr_ch->sim_bdry_contact_y = true;

          // ptr_ch->pbc_crosses_sim_box_bdry = true;
          // ptr_ch->pbc_crosses_sim_box_bdry_y = true;
        }
        else if (new_box_min_y == 0)
        {
          new_sim_bdry_contact = true;
          new_sim_bdry_contact_y = true;

          // ptr_ch->sim_bdry_contact = true;
          // ptr_ch->sim_bdry_contact_y = true;
        }
        else if (new_box_max_y <= -1)
        {
          if (new_box_max_y == -1)
          {
            new_sim_bdry_contact = true;
            new_sim_bdry_contact_y = true;

            // ptr_ch->sim_bdry_contact = true;
            // ptr_ch->sim_bdry_contact_y = true;
          }

          new_box_min_y += (1 << lv);
          new_box_max_y += (1 << lv);
        }
      } // End else if(ptr_node->pbc_crosses_sim_box_bdry = true)
      // Case Parent x axis do not cross the simulation box
      else
      {
        for (int j = 0; j < ptr_node->ptr_zone_size[new_zone_idx]; j++)
        {
          cell_idx = ptr_node->pptr_zones[new_zone_idx][j]; // Cell index in the parent zone
          if (new_box_min_y > ptr_node->ptr_cell_idx_y[cell_idx])
          {
            new_box_min_y = ptr_node->ptr_cell_idx_y[cell_idx];
          }
          if (new_box_max_y < ptr_node->ptr_cell_idx_y[cell_idx])
          {
            new_box_max_y = ptr_node->ptr_cell_idx_y[cell_idx];
          }
        }

        // Analysis of min and max of the refinement zone
        if (new_box_min_y == 0 || new_box_max_y == (1 << lv) - 1)
        {
          new_sim_bdry_contact = true;
          new_sim_bdry_contact_y = true;

          // ptr_ch->sim_bdry_contact = true;
          // ptr_ch->sim_bdry_contact_y = true;
        }
      }

      //* >> Z axis
      if (ptr_node->pbc_crosses_whole_sim_box_z == true)
      {
        if (ptr_node->ptr_pbc_bool_bdry_anomalies_z[new_zone_idx] == false)
        {
          for (int j = 0; j < ptr_node->ptr_zone_size[new_zone_idx]; j++)
          {
            cell_idx = ptr_node->pptr_zones[new_zone_idx][j]; // Cell index in the parent zone
            if (new_box_min_z > ptr_node->ptr_cell_idx_z[cell_idx])
            {
              new_box_min_z = ptr_node->ptr_cell_idx_z[cell_idx];
            }
            if (new_box_max_z < ptr_node->ptr_cell_idx_z[cell_idx])
            {
              new_box_max_z = ptr_node->ptr_cell_idx_z[cell_idx];
            }
          }
        }
        else if (ptr_node->pbc_subzones_size == 1)
        {
          new_box_min_z = 0;
          new_box_max_z = (1 << lv) - 1;
        }
        else
        {
          // Matching the subzones
          aux_subzones_analized = ptr_node->pbc_subzones_size;
          check = true;
          while (aux_subzones_analized > 1 && check == true)
          {
            check = false;
            for (int i = 0; i < 2; i++)
            {
              for (int j = i + 1; j < aux_subzones_analized; j++)
              {
                if ((ptr_node->ptr_pbc_max_subzones_z[j] + 1 >= ptr_node->ptr_pbc_min_subzones_z[i] && ptr_node->ptr_pbc_min_subzones_z[j] <= 1 + ptr_node->ptr_pbc_max_subzones_z[i]))
                {
                  ptr_node->ptr_pbc_min_subzones_z[i] = ptr_node->ptr_pbc_min_subzones_z[i] < ptr_node->ptr_pbc_min_subzones_z[j] ? ptr_node->ptr_pbc_min_subzones_z[i] : ptr_node->ptr_pbc_min_subzones_z[j];
                  ptr_node->ptr_pbc_max_subzones_z[i] = ptr_node->ptr_pbc_max_subzones_z[i] > ptr_node->ptr_pbc_max_subzones_z[j] ? ptr_node->ptr_pbc_max_subzones_z[i] : ptr_node->ptr_pbc_max_subzones_z[j];

                  // Reducing the total amount of analized subzones
                  ptr_node->ptr_pbc_min_subzones_z[j] = ptr_node->ptr_pbc_min_subzones_z[aux_subzones_analized - 1];
                  ptr_node->ptr_pbc_max_subzones_z[j] = ptr_node->ptr_pbc_max_subzones_z[aux_subzones_analized - 1];
                  aux_subzones_analized--;
                  j--;
                  check = true;
                }
              }
            }
          }

          // printf("aux_subzones_analized = %d\n", aux_subzones_analized);

          // min and max between all subzones
          if (aux_subzones_analized == 1)
          {
            new_box_min_z = 0;
            new_box_max_z = (1 << lv) - 1;
          }
          else
          {
            new_box_min_z = ptr_node->ptr_pbc_min_subzones_z[0] == 0 ? ptr_node->ptr_pbc_min_subzones_z[1] - (1 << lv) : ptr_node->ptr_pbc_min_subzones_z[0] - (1 << lv);
            new_box_max_z = ptr_node->ptr_pbc_min_subzones_z[0] == 0 ? ptr_node->ptr_pbc_max_subzones_z[0] : ptr_node->ptr_pbc_max_subzones_z[1];
          }
        }

        // Analysis of min and max of the refinement zone
        // Case refinement zone crosses the whole simulation box
        if (new_box_max_z - new_box_min_z == (1 << lv) - 1)
        {
          new_sim_bdry_contact = true;
          new_sim_bdry_contact_z = true;

          new_pbc_crosses_sim_box_bdry = true;
          new_pbc_crosses_sim_box_bdry_z = true;

          new_pbc_crosses_whole_sim_box = true;
          new_pbc_crosses_whole_sim_box_z = true;

          // ptr_ch->sim_bdry_contact = true;
          // ptr_ch->sim_bdry_contact_z = true;

          // ptr_ch->pbc_crosses_sim_box_bdry = true;
          // ptr_ch->pbc_crosses_sim_box_bdry_z = true;

          // ptr_ch->pbc_crosses_whole_sim_box = true;
          // ptr_ch->pbc_crosses_whole_sim_box_z = true;
        }
        // Case refinement zone crosses the simulation boundary
        else if (new_box_min_z < 0)
        {
          new_sim_bdry_contact = true;
          new_sim_bdry_contact_z = true;

          new_pbc_crosses_sim_box_bdry = true;
          new_pbc_crosses_sim_box_bdry_z = true;

          // ptr_ch->sim_bdry_contact = true;
          // ptr_ch->sim_bdry_contact_z = true;

          // ptr_ch->pbc_crosses_sim_box_bdry = true;
          // ptr_ch->pbc_crosses_sim_box_bdry_z = true;
        }
        // Case refinement zone touches the simulation boundary
        else if (new_box_min_z == 0 || new_box_max_z == (1 << lv) - 1)
        {
          new_sim_bdry_contact = true;
          new_sim_bdry_contact_z = true;

          // ptr_ch->sim_bdry_contact = true;
          // ptr_ch->sim_bdry_contact_z = true;
        }

      } // End if(ptr_node->pbc_crosses_whole_sim_box_z == true)
      // Case Parent x axis crosses the simulation box
      else if (ptr_node->pbc_crosses_sim_box_bdry_z == true)
      {
        for (int j = 0; j < ptr_node->ptr_zone_size[new_zone_idx]; j++)
        {
          cell_idx = ptr_node->pptr_zones[new_zone_idx][j]; // Cell index in the parent zone

          if (ptr_node->ptr_cell_idx_z[cell_idx] > ptr_node->box_max_z)
          {
            new_box_min_z = new_box_min_z < ptr_node->ptr_cell_idx_z[cell_idx] - (1 << lv) ? new_box_min_z : ptr_node->ptr_cell_idx_z[cell_idx] - (1 << lv);
            new_box_max_z = new_box_max_z > ptr_node->ptr_cell_idx_z[cell_idx] - (1 << lv) ? new_box_max_z : ptr_node->ptr_cell_idx_z[cell_idx] - (1 << lv);
          }
          else
          {
            new_box_min_z = new_box_min_z < ptr_node->ptr_cell_idx_z[cell_idx] ? new_box_min_z : ptr_node->ptr_cell_idx_z[cell_idx];
            new_box_max_z = new_box_max_z > ptr_node->ptr_cell_idx_z[cell_idx] ? new_box_max_z : ptr_node->ptr_cell_idx_z[cell_idx];
          }
        }

        // Analysis of min and max of the refinement zone
        if (new_box_min_z < 0 && new_box_max_z >= 0)
        {
          new_sim_bdry_contact = true;
          new_sim_bdry_contact_z = true;

          new_pbc_crosses_sim_box_bdry = true;
          new_pbc_crosses_sim_box_bdry_z = true;

          // ptr_ch->sim_bdry_contact = true;
          // ptr_ch->sim_bdry_contact_z = true;

          // ptr_ch->pbc_crosses_sim_box_bdry = true;
          // ptr_ch->pbc_crosses_sim_box_bdry_z = true;
        }
        else if (new_box_min_z == 0)
        {
          new_sim_bdry_contact = true;
          new_sim_bdry_contact_z = true;

          // ptr_ch->sim_bdry_contact = true;
          // ptr_ch->sim_bdry_contact_z = true;
        }
        else if (new_box_max_z <= -1)
        {
          if (new_box_max_z == -1)
          {
            new_sim_bdry_contact = true;
            new_sim_bdry_contact_z = true;

            // ptr_ch->sim_bdry_contact = true;
            // ptr_ch->sim_bdry_contact_z = true;
          }

          new_box_min_z += (1 << lv);
          new_box_max_z += (1 << lv);
        }
      } // End else if(ptr_node->pbc_crosses_sim_box_bdry = true)
      // Case Parent x axis do not cross the simulation box
      else
      {
        for (int j = 0; j < ptr_node->ptr_zone_size[new_zone_idx]; j++)
        {
          cell_idx = ptr_node->pptr_zones[new_zone_idx][j]; // Cell index in the parent zone
          if (new_box_min_z > ptr_node->ptr_cell_idx_z[cell_idx])
          {
            new_box_min_z = ptr_node->ptr_cell_idx_z[cell_idx];
          }
          if (new_box_max_z < ptr_node->ptr_cell_idx_z[cell_idx])
          {
            new_box_max_z = ptr_node->ptr_cell_idx_z[cell_idx];
          }
        }

        // Analysis of min and max of the refinement zone
        if (new_box_min_z == 0 || new_box_max_z == (1 << lv) - 1)
        {
          new_sim_bdry_contact = true;
          new_sim_bdry_contact_z = true;

          // ptr_ch->sim_bdry_contact = true;
          // ptr_ch->sim_bdry_contact_z = true;
        }
      }

      /////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////

      // Changing the min and max of the "Smallest Box" from parent units to child units
      new_box_min_x = 2 * new_box_min_x;
      new_box_min_y = 2 * new_box_min_y;
      new_box_min_z = 2 * new_box_min_z;
      new_box_max_x = 2 * new_box_max_x + 1;
      new_box_max_y = 2 * new_box_max_y + 1;
      new_box_max_z = 2 * new_box_max_z + 1;

      aux_fix_min_x = new_box_min_x - ptr_ch->box_ts_x;
      aux_fix_min_y = new_box_min_y - ptr_ch->box_ts_y;
      aux_fix_min_z = new_box_min_z - ptr_ch->box_ts_z;

      aux_fix_max_x = new_box_max_x - ptr_ch->box_ts_x;
      aux_fix_max_y = new_box_max_y - ptr_ch->box_ts_y;
      aux_fix_max_z = new_box_max_z - ptr_ch->box_ts_z;

      if (bdry_cond_type == 0 && ptr_node->pbc_crosses_sim_box_bdry == true)
      {
        //  if (((ptr_ch->pbc_crosses_whole_sim_box_x == true && new_pbc_crosses_sim_box_bdry_x == true && new_pbc_crosses_whole_sim_box_x == false) ||
        //      (new_pbc_crosses_whole_sim_box_x == true && ptr_ch->pbc_crosses_sim_box_bdry_x == true && ptr_ch->pbc_crosses_whole_sim_box_x == false)) == false)
        // if (ptr_ch->box_real_dim_x <= (1 << (lv + 1)))
        //{
        if (ptr_ch->box_min_x < 0 && new_box_min_x >= 0)
        {
          if (new_box_min_x - ptr_ch->box_max_x > ptr_ch->box_min_x - new_box_max_x + (1 << (lv + 1)))
          {
            aux_fix_min_x -= (1 << (lv + 1));
            aux_fix_max_x -= (1 << (lv + 1));
            ptr_ch->box_ts_x += (1 << (lv + 1));
            movement_box_x = true;
          }
        }
        else if (ptr_ch->box_min_x >= 0 && new_box_min_x < 0)
        {
          if (ptr_ch->box_min_x - new_box_max_x > new_box_min_x - ptr_ch->box_max_x + (1 << (lv + 1)))
          {
            aux_fix_min_x += (1 << (lv + 1));
            aux_fix_max_x += (1 << (lv + 1));
            ptr_ch->box_ts_x -= (1 << (lv + 1));
            movement_box_x = true;
          }
        }
        //}

        // if (((ptr_ch->pbc_crosses_whole_sim_box_y == true && new_pbc_crosses_sim_box_bdry_y == true && new_pbc_crosses_whole_sim_box_y == false) ||
        //      (new_pbc_crosses_whole_sim_box_y == true && ptr_ch->pbc_crosses_sim_box_bdry_y == true && ptr_ch->pbc_crosses_whole_sim_box_y == false)) == false)
        // {
        if (ptr_ch->box_min_y < 0 && new_box_min_y >= 0)
        {
          if (new_box_min_y - ptr_ch->box_max_y > ptr_ch->box_min_y - new_box_max_y + (1 << (lv + 1)))
          {
            aux_fix_min_y -= (1 << (lv + 1));
            aux_fix_max_y -= (1 << (lv + 1));
            ptr_ch->box_ts_y += (1 << (lv + 1));
            movement_box_y = true;
          }
        }
        else if (ptr_ch->box_min_y >= 0 && new_box_min_y < 0)
        {
          if (ptr_ch->box_min_y - new_box_max_y > new_box_min_y - ptr_ch->box_max_y + (1 << (lv + 1)))
          {
            aux_fix_min_y += (1 << (lv + 1));
            aux_fix_max_y += (1 << (lv + 1));
            ptr_ch->box_ts_y -= (1 << (lv + 1));
            movement_box_y = true;
          }
        }
        // }

        // if (((ptr_ch->pbc_crosses_whole_sim_box_z == true && new_pbc_crosses_sim_box_bdry_z == true && new_pbc_crosses_whole_sim_box_z == false) ||
        //      (new_pbc_crosses_whole_sim_box_z == true && ptr_ch->pbc_crosses_sim_box_bdry_z == true && ptr_ch->pbc_crosses_whole_sim_box_z == false)) == false)
        // {
        if (ptr_ch->box_min_z < 0 && new_box_min_z >= 0)
        {
          if (new_box_min_z - ptr_ch->box_max_z > ptr_ch->box_min_z - new_box_max_z + (1 << (lv + 1)))
          {
            aux_fix_min_z -= (1 << (lv + 1));
            aux_fix_max_z -= (1 << (lv + 1));
            ptr_ch->box_ts_z += (1 << (lv + 1));
            movement_box_z = true;
          }
        }
        else if (ptr_ch->box_min_z >= 0 && new_box_min_z < 0)
        {
          if (ptr_ch->box_min_z - new_box_max_z > new_box_min_z - ptr_ch->box_max_z + (1 << (lv + 1)))
          {
            aux_fix_min_z += (1 << (lv + 1));
            aux_fix_max_z += (1 << (lv + 1));
            ptr_ch->box_ts_z -= (1 << (lv + 1));
            movement_box_z = true;
          }
        }
        //}
      }

      if (aux_fix_min_x < bder_box ||
          aux_fix_min_y < bder_box ||
          aux_fix_min_z < bder_box ||
          aux_fix_max_x > ptr_ch->box_real_dim_x - bder_box ||
          aux_fix_max_y > ptr_ch->box_real_dim_y - bder_box ||
          aux_fix_max_z > ptr_ch->box_real_dim_z - bder_box)
      {
        ptr_ch->box_check_fit = false;
      }

      // Correction flag in case the new cells require recompute his box indices when they stay in the box
      if (bdry_cond_type == 0 && ptr_ch->box_check_fit == true)
      {
        // if (((ptr_ch->pbc_crosses_whole_sim_box == true || ptr_ch->pbc_crosses_sim_box_bdry == false) && new_pbc_crosses_sim_box_bdry == true && new_pbc_crosses_whole_sim_box == false) ||
        //     ((new_pbc_crosses_whole_sim_box == true || new_pbc_crosses_sim_box_bdry == false) && ptr_ch->pbc_crosses_sim_box_bdry == true && ptr_ch->pbc_crosses_whole_sim_box == false))
        // {
        //     ptr_ch->pbc_correction_due_pbc_flag = true;
        // }

        // X axis
        if (new_box_min_x >= 0 && ptr_ch->box_min_x < 0)
        {
          if (new_box_min_x <= ptr_ch->box_max_x && new_box_max_x >= ptr_ch->box_min_x + (1 << (lv + 1)))
          {
            ptr_ch->pbc_correction_due_pbc_flag = true;
          }
        }
        else if (new_box_min_x < 0 && ptr_ch->box_min_x >= 0)
        {
          if (ptr_ch->box_min_x <= new_box_max_x && ptr_ch->box_max_x >= new_box_min_x + (1 << (lv + 1)))
          {
            ptr_ch->pbc_correction_due_pbc_flag = true;
          }
        }
        else if (new_box_min_x < 0 && ptr_ch->box_min_x < 0)
        {
          if (new_box_min_x + (1 << (lv + 1)) <= ptr_ch->box_max_x || ptr_ch->box_min_x + (1 << (lv + 1)) <= new_box_max_x)
          {
            ptr_ch->pbc_correction_due_pbc_flag = true;
          }
        }

        // Y axis
        if (new_box_min_y >= 0 && ptr_ch->box_min_y < 0)
        {
          if (new_box_min_y <= ptr_ch->box_max_y && new_box_max_y >= ptr_ch->box_min_y + (1 << (lv + 1)))
          {
            ptr_ch->pbc_correction_due_pbc_flag = true;
          }
        }
        else if (new_box_min_y < 0 && ptr_ch->box_min_y >= 0)
        {
          if (ptr_ch->box_min_y <= new_box_max_y && ptr_ch->box_max_y >= new_box_min_y + (1 << (lv + 1)))
          {
            ptr_ch->pbc_correction_due_pbc_flag = true;
          }
        }
        else if (new_box_min_y < 0 && ptr_ch->box_min_y < 0)
        {
          if (new_box_min_y + (1 << (lv + 1)) <= ptr_ch->box_max_y || ptr_ch->box_min_y + (1 << (lv + 1)) <= new_box_max_y)
          {
            ptr_ch->pbc_correction_due_pbc_flag = true;
          }
        }

        // Z axis
        if (new_box_min_z >= 0 && ptr_ch->box_min_z < 0)
        {
          if (new_box_min_z <= ptr_ch->box_max_z && new_box_max_z >= ptr_ch->box_min_z + (1 << (lv + 1)))
          {
            ptr_ch->pbc_correction_due_pbc_flag = true;
          }
        }
        else if (new_box_min_z < 0 && ptr_ch->box_min_z >= 0)
        {
          if (ptr_ch->box_min_z <= new_box_max_z && ptr_ch->box_max_z >= new_box_min_z + (1 << (lv + 1)))
          {
            ptr_ch->pbc_correction_due_pbc_flag = true;
          }
        }
        else if (new_box_min_z < 0 && ptr_ch->box_min_z < 0)
        {
          if (new_box_min_z + (1 << (lv + 1)) <= ptr_ch->box_max_z || ptr_ch->box_min_z + (1 << (lv + 1)) <= new_box_max_z)
          {
            ptr_ch->pbc_correction_due_pbc_flag = true;
          }
        }
      }

      // Removing -6 and -5 to the box boundary when it corresponds
      if (ptr_ch->box_check_fit == true && ptr_ch->sim_bdry_contact == true)
      {

        box_real_dim_X_ch = ptr_ch->box_real_dim_x;
        box_real_dim_X_times_Y_ch = ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;

        // Removing status of -6
        if (bdry_cond_type == 0 && ptr_ch->pbc_crosses_whole_sim_box == true)
        {
          // X axis
          if (ptr_ch->pbc_crosses_whole_sim_box_x == true && new_pbc_crosses_whole_sim_box_x == false)
          {

            aux_box_ts = ptr_ch->box_ts_x;

            if (movement_box_x == true)
            {
              aux_box_ts += new_box_min_x < 0 ? (1 << (lv + 1)) : -(1 << (lv + 1));
            }

            for (int k = 0; k < ptr_ch->box_real_dim_z; k++)
            {
              for (int j = 0; j < ptr_ch->box_real_dim_y; j++)
              {
                // for (int i = 0; i < (ptr_ch->box_real_dim_x - ptr_ch->box_dim_x) / 2; i++)
                for (int i = 0; i < ptr_ch->box_min_x - aux_box_ts; i++)
                {
                  box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
                  ptr_ch->ptr_box[box_idx_ch] = -4;
                }

                // for (int i = (ptr_ch->box_real_dim_x + ptr_ch->box_dim_x) / 2; i < ptr_ch->box_real_dim_x; i++)
                for (int i = ptr_ch->box_max_x + 1 - aux_box_ts; i < ptr_ch->box_real_dim_x; i++)
                {
                  box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
                  ptr_ch->ptr_box[box_idx_ch] = -4;
                }
              }
            }
          }

          // Y axis
          if (ptr_ch->pbc_crosses_whole_sim_box_y == true && new_pbc_crosses_whole_sim_box_y == false)
          {

            aux_box_ts = ptr_ch->box_ts_y;

            if (movement_box_y == true)
            {
              aux_box_ts += new_box_min_y < 0 ? (1 << (lv + 1)) : -(1 << (lv + 1));
            }

            for (int k = 0; k < ptr_ch->box_real_dim_z; k++)
            {
              for (int i = 0; i < ptr_ch->box_real_dim_x; i++)
              {
                // for (int j = 0; j < (ptr_ch->box_real_dim_y - ptr_ch->box_dim_y) / 2; j++)
                for (int j = 0; j < ptr_ch->box_min_y - aux_box_ts; j++)
                {
                  box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
                  ptr_ch->ptr_box[box_idx_ch] = -4;
                }

                // for (int j = (ptr_ch->box_real_dim_y + ptr_ch->box_dim_y) / 2; j < ptr_ch->box_real_dim_y; j++)
                for (int j = ptr_ch->box_max_y + 1 - aux_box_ts; j < ptr_ch->box_real_dim_y; j++)
                {
                  box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
                  ptr_ch->ptr_box[box_idx_ch] = -4;
                }
              }
            }
          }

          // Z axis
          if (ptr_ch->pbc_crosses_whole_sim_box_z == true && new_pbc_crosses_whole_sim_box_z == false)
          {

            aux_box_ts = ptr_ch->box_ts_z;

            if (movement_box_z == true)
            {
              aux_box_ts += new_box_min_z < 0 ? (1 << (lv + 1)) : -(1 << (lv + 1));
            }

            for (int j = 0; j < ptr_ch->box_real_dim_y; j++)
            {
              for (int i = 0; i < ptr_ch->box_real_dim_x; i++)
              {
                // for (int k = 0; k < (ptr_ch->box_real_dim_z - ptr_ch->box_dim_z) / 2; k++)
                for (int k = 0; k < ptr_ch->box_min_z - aux_box_ts; k++)
                {
                  box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
                  ptr_ch->ptr_box[box_idx_ch] = -4;
                }

                // for (int k = (ptr_ch->box_real_dim_z + ptr_ch->box_dim_z) / 2; k < ptr_ch->box_real_dim_z; k++)
                for (int k = ptr_ch->box_max_z + 1 - aux_box_ts; k < ptr_ch->box_real_dim_z; k++)
                {
                  box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
                  ptr_ch->ptr_box[box_idx_ch] = -4;
                }
              }
            }
          }
        }
        // Removing status of -5
        else if (bdry_cond_type != 0)
        {
          // X axis
          if (ptr_ch->sim_bdry_contact_x == true && new_sim_bdry_contact_x == false)
          {

            if (ptr_ch->box_min_x == 0)
            {
              for (int k = 0; k < ptr_ch->box_real_dim_z; k++)
              {
                for (int j = 0; j < ptr_ch->box_real_dim_y; j++)
                {
                  // for (int i = 0; i < (ptr_ch->box_real_dim_x - ptr_ch->box_dim_x) / 2; i++)
                  for (int i = 0; i < ptr_ch->box_min_x - ptr_ch->box_ts_x; i++)
                  {
                    box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
                    ptr_ch->ptr_box[box_idx_ch] = -4;
                  }
                }
              }
            }

            if (ptr_ch->box_max_x == (1 << (lv + 1)) - 1)
            {
              for (int k = 0; k < ptr_ch->box_real_dim_z; k++)
              {
                for (int j = 0; j < ptr_ch->box_real_dim_y; j++)
                {
                  // for (int i = (ptr_ch->box_real_dim_x + ptr_ch->box_dim_x) / 2; i < ptr_ch->box_real_dim_x; i++)
                  for (int i = ptr_ch->box_max_x + 1 - ptr_ch->box_ts_x; i < ptr_ch->box_real_dim_x; i++)
                  {
                    box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
                    ptr_ch->ptr_box[box_idx_ch] = -4;
                  }
                }
              }
            }
          }

          // Y axis
          if (ptr_ch->sim_bdry_contact_y == true && new_sim_bdry_contact_y == false)
          {

            if (ptr_ch->box_min_y == 0)
            {
              for (int k = 0; k < ptr_ch->box_real_dim_z; k++)
              {
                for (int i = 0; i < ptr_ch->box_real_dim_x; i++)
                {
                  // for (int j = 0; j < (ptr_ch->box_real_dim_y - ptr_ch->box_dim_y) / 2; j++)
                  for (int j = 0; j < ptr_ch->box_min_y - ptr_ch->box_ts_y; j++)
                  {
                    box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
                    ptr_ch->ptr_box[box_idx_ch] = -4;
                  }
                }
              }
            }

            if (ptr_ch->box_max_y == (1 << (lv + 1)) - 1)
            {
              for (int k = 0; k < ptr_ch->box_real_dim_z; k++)
              {
                for (int i = 0; i < ptr_ch->box_real_dim_x; i++)
                {
                  // for (int j = (ptr_ch->box_real_dim_y + ptr_ch->box_dim_y) / 2; j < ptr_ch->box_real_dim_y; j++)
                  for (int j = ptr_ch->box_max_y + 1 - ptr_ch->box_ts_y; j < ptr_ch->box_real_dim_y; j++)
                  {
                    box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
                    ptr_ch->ptr_box[box_idx_ch] = -4;
                  }
                }
              }
            }
          }

          // Z axis
          if (ptr_ch->sim_bdry_contact_z == true && new_sim_bdry_contact_z == false)
          {
            if (ptr_ch->box_min_z == 0)
            {
              for (int j = 0; j < ptr_ch->box_real_dim_y; j++)
              {
                for (int i = 0; i < ptr_ch->box_real_dim_x; i++)
                {
                  // for (int k = 0; k < (ptr_ch->box_real_dim_z - ptr_ch->box_dim_z) / 2; k++)
                  for (int k = 0; k < ptr_ch->box_min_z - ptr_ch->box_ts_z; k++)
                  {
                    box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
                    ptr_ch->ptr_box[box_idx_ch] = -4;
                  }
                }
              }
            }

            if (ptr_ch->box_max_z == (1 << (lv + 1)) - 1)
            {
              for (int j = 0; j < ptr_ch->box_real_dim_y; j++)
              {
                for (int i = 0; i < ptr_ch->box_real_dim_x; i++)
                {
                  // for (int k = (ptr_ch->box_real_dim_z + ptr_ch->box_dim_z) / 2; k < ptr_ch->box_real_dim_z; k++)
                  for (int k = ptr_ch->box_max_z + 1 - ptr_ch->box_ts_z; k < ptr_ch->box_real_dim_z; k++)
                  {
                    box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
                    ptr_ch->ptr_box[box_idx_ch] = -4;
                  }
                }
              }
            }
          }
        }
      }

      // if (new_sim_bdry_contact == true)
      // {
      // printf("check fit = %d\n", ptr_ch->box_check_fit);
      // printf("old box_min_x = %d, box_max_x = %d\n", ptr_ch->box_min_x, ptr_ch->box_max_x);
      // printf("old box_min_y = %d, box_max_y = %d\n", ptr_ch->box_min_y, ptr_ch->box_max_y);
      // printf("old box_min_z = %d, box_max_z = %d\n\n", ptr_ch->box_min_z, ptr_ch->box_max_z);
      // printf("new box_min_x = %d, box_max_x = %d\n", new_box_min_x, new_box_max_x);
      // printf("new box_min_y = %d, box_max_y = %d\n", new_box_min_y, new_box_max_y);
      // printf("new box_min_z = %d, box_max_z = %d\n\n", new_box_min_z, new_box_max_z);
      //}

      //* >> Updating the The MIN and MAX of the set containig only the new cells to be refined
      ptr_ch->box_min_x = new_box_min_x;
      ptr_ch->box_min_y = new_box_min_y;
      ptr_ch->box_min_z = new_box_min_z;
      ptr_ch->box_max_x = new_box_max_x;
      ptr_ch->box_max_y = new_box_max_y;
      ptr_ch->box_max_z = new_box_max_z;

      // Size of the new "Smallest Box"
      ptr_ch->box_dim_x = new_box_max_x - new_box_min_x + 1;
      ptr_ch->box_dim_y = new_box_max_y - new_box_min_y + 1;
      ptr_ch->box_dim_z = new_box_max_z - new_box_min_z + 1;

      //* >> The new box does not fit in the old box *//
      if (ptr_ch->box_check_fit == false)
      {
        // Real dimensions of the old box
        // ptr_ch->box_real_dim_x_old = ptr_ch->box_real_dim_x;
        // ptr_ch->box_real_dim_y_old = ptr_ch->box_real_dim_y;
        // ptr_ch->box_real_dim_z_old = ptr_ch->box_real_dim_z;

        // Translations between cell array and old box (Used to access only to the mass box)
        // ptr_ch->box_ts_x_old = ptr_ch->box_ts_x; // Every cell in the level l in the box must be subtracted this value to obtain the box index
        // ptr_ch->box_ts_y_old = ptr_ch->box_ts_y;
        // ptr_ch->box_ts_z_old = ptr_ch->box_ts_z;

        // Real dimensions of the new box
        ptr_ch->box_real_dim_x = 5 > (n_exp - 1) ? (ptr_ch->box_dim_x + 10) : (ptr_ch->box_dim_x + 2 * n_exp - 2);
        ptr_ch->box_real_dim_y = 5 > (n_exp - 1) ? (ptr_ch->box_dim_y + 10) : (ptr_ch->box_dim_y + 2 * n_exp - 2);
        ptr_ch->box_real_dim_z = 5 > (n_exp - 1) ? (ptr_ch->box_dim_z + 10) : (ptr_ch->box_dim_z + 2 * n_exp - 2);

        // Translations between cell array and new box
        pos_x = (ptr_ch->box_real_dim_x - ptr_ch->box_dim_x) / 2; // Half of the distance of the box side less the "Smallest Box" side
        ptr_ch->box_ts_x = new_box_min_x - pos_x;                 // Every cell in the level l in the box must be subtracted this value to obtain the box index
        pos_y = (ptr_ch->box_real_dim_y - ptr_ch->box_dim_y) / 2;
        ptr_ch->box_ts_y = new_box_min_y - pos_y;
        pos_z = (ptr_ch->box_real_dim_z - ptr_ch->box_dim_z) / 2;
        ptr_ch->box_ts_z = new_box_min_z - pos_z;

        // //* >> Transfer from cell structure to auxiliary cell structure *//
        // ptr_cell_struct_aux = ptr_ch->ptr_cell_struct_old;
        // ptr_ch->ptr_cell_struct_old = ptr_ch->ptr_cell_struct;
        // ptr_ch->ptr_cell_struct = ptr_cell_struct_aux;

        // size = ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y * ptr_ch->box_real_dim_z;

        // //* >> Space checking of boxes and cell structures
        // if (space_check(&(ptr_ch->box_cap), size, 1.0f, "p4i1i1c1c1", &(ptr_ch->ptr_box), &(ptr_ch->ptr_box_old), &(ptr_ch->ptr_cell_struct), &(ptr_ch->ptr_cell_struct_old)) == _FAILURE_)
        // {
        //     printf("Error, in space_check function\n");
        //     return _FAILURE_;
        // }

        //* >> Transfer from cell structure to auxiliary cell structure *//
        //* >> Space checking of cell structure old
        if (space_check(&(ptr_ch->cell_struct_old_cap), ptr_ch->box_cap, 1.0f, "p1c1", &(ptr_ch->ptr_cell_struct_old)) == _FAILURE_)
        {
          printf("Error, in space_check function\n");
          return _FAILURE_;
        }

        ptr_cell_struct_aux = ptr_ch->ptr_cell_struct_old;
        ptr_ch->ptr_cell_struct_old = ptr_ch->ptr_cell_struct;
        ptr_ch->ptr_cell_struct = ptr_cell_struct_aux;

        size = ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y * ptr_ch->box_real_dim_z;

        //* >> Space checking of boxes and cell structures
        if (space_check(&(ptr_ch->box_cap), size, 1.0f, "p3i1i1c1", &(ptr_ch->ptr_box), &(ptr_ch->ptr_box_old), &(ptr_ch->ptr_cell_struct)) == _FAILURE_)
        {
          printf("Error, in space_check function\n");
          return _FAILURE_;
        }

        // Initializing box at NO-EXIST (-4)
        for (int j = 0; j < size; j++)
        {
          ptr_ch->ptr_box[j] = -4;
        }

        //* Potential, Acceleration and density capacity of the grid *//
        size = (ptr_ch->box_real_dim_x + 1) * (ptr_ch->box_real_dim_y + 1) * (ptr_ch->box_real_dim_z + 1);
        //* >> Space checking
        if (space_check(&(ptr_ch->grid_properties_cap), size, 1.0f, "p6v1v1v1v1v1v1", &(ptr_ch->ptr_pot), &(ptr_ch->ptr_pot_old), &(ptr_ch->ptr_ax), &(ptr_ch->ptr_ay), &(ptr_ch->ptr_az), &(ptr_ch->ptr_d)) == _FAILURE_)
        {
          printf("Error, in space_check function\n");
          return _FAILURE_;
        }
      }

      // //* >> Adding -6 and -5 to the box boundary when it corresponds
      // if (new_sim_bdry_contact == true)
      // {

      //     box_real_dim_X_ch = ptr_ch->box_real_dim_x;
      //     box_real_dim_X_times_Y_ch = ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;

      //     // Adding status of -6
      //     if (bdry_cond_type == 0 && new_pbc_crosses_whole_sim_box == true)
      //     {
      //         // X axis
      //         if (ptr_ch->pbc_crosses_whole_sim_box_x == false && new_pbc_crosses_whole_sim_box_x == true)
      //         {
      //             for (int k = 0; k < ptr_ch->box_real_dim_z; k++)
      //             {
      //                 for (int j = 0; j < ptr_ch->box_real_dim_y; j++)
      //                 {
      //                     for (int i = 0; i < (ptr_ch->box_real_dim_x - ptr_ch->box_dim_x) / 2; i++)
      //                     {
      //                         box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
      //                         ptr_ch->ptr_box[box_idx_ch] = -6;
      //                     }

      //                     for (int i = (ptr_ch->box_real_dim_x + ptr_ch->box_dim_x) / 2; i < ptr_ch->box_real_dim_x; i++)
      //                     {
      //                         box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
      //                         ptr_ch->ptr_box[box_idx_ch] = -6;
      //                     }
      //                 }
      //             }
      //         }

      //         // Y axis
      //         if (ptr_ch->pbc_crosses_whole_sim_box_y == false && new_pbc_crosses_whole_sim_box_y == true)
      //         {
      //             for (int k = 0; k < ptr_ch->box_real_dim_z; k++)
      //             {
      //                 for (int i = 0; i < ptr_ch->box_real_dim_x; i++)
      //                 {
      //                     for (int j = 0; j < (ptr_ch->box_real_dim_y - ptr_ch->box_dim_y) / 2; j++)
      //                     {
      //                         box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
      //                         ptr_ch->ptr_box[box_idx_ch] = -6;
      //                     }

      //                     for (int j = (ptr_ch->box_real_dim_y + ptr_ch->box_dim_y) / 2; j < ptr_ch->box_real_dim_y; j++)
      //                     {
      //                         box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
      //                         ptr_ch->ptr_box[box_idx_ch] = -6;
      //                     }
      //                 }
      //             }
      //         }

      //         // Z axis
      //         if (ptr_ch->pbc_crosses_whole_sim_box_z == false && new_pbc_crosses_whole_sim_box_z == true)
      //         {
      //             for (int j = 0; j < ptr_ch->box_real_dim_y; j++)
      //             {
      //                 for (int i = 0; i < ptr_ch->box_real_dim_x; i++)
      //                 {
      //                     for (int k = 0; k < (ptr_ch->box_real_dim_z - ptr_ch->box_dim_z) / 2; k++)
      //                     {
      //                         box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
      //                         ptr_ch->ptr_box[box_idx_ch] = -6;
      //                     }

      //                     for (int k = (ptr_ch->box_real_dim_z + ptr_ch->box_dim_z) / 2; k < ptr_ch->box_real_dim_z; k++)
      //                     {
      //                         box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
      //                         ptr_ch->ptr_box[box_idx_ch] = -6;
      //                     }
      //                 }
      //             }
      //         }
      //     }
      //     // Removing status of -5
      //     else if (bdry_cond_type != 0)
      //     {
      //         // X axis
      //         if (ptr_ch->sim_bdry_contact_x == false && new_sim_bdry_contact_x == true)
      //         {
      //             for (int k = 0; k < ptr_ch->box_real_dim_z; k++)
      //             {
      //                 for (int j = 0; j < ptr_ch->box_real_dim_y; j++)
      //                 {
      //                     for (int i = 0; i < (ptr_ch->box_real_dim_x - ptr_ch->box_dim_x) / 2; i++)
      //                     {
      //                         box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
      //                         ptr_ch->ptr_box[box_idx_ch] = -5;
      //                     }

      //                     for (int i = (ptr_ch->box_real_dim_x + ptr_ch->box_dim_x) / 2; i < ptr_ch->box_real_dim_x; i++)
      //                     {
      //                         box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
      //                         ptr_ch->ptr_box[box_idx_ch] = -5;
      //                     }
      //                 }
      //             }
      //         }

      //         // Y axis
      //         if (ptr_ch->sim_bdry_contact_y == false && new_sim_bdry_contact_y == true)
      //         {
      //             for (int k = 0; k < ptr_ch->box_real_dim_z; k++)
      //             {
      //                 for (int i = 0; i < ptr_ch->box_real_dim_x; i++)
      //                 {
      //                     for (int j = 0; j < (ptr_ch->box_real_dim_y - ptr_ch->box_dim_y) / 2; j++)
      //                     {
      //                         box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
      //                         ptr_ch->ptr_box[box_idx_ch] = -5;
      //                     }

      //                     for (int j = (ptr_ch->box_real_dim_y + ptr_ch->box_dim_y) / 2; j < ptr_ch->box_real_dim_y; j++)
      //                     {
      //                         box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
      //                         ptr_ch->ptr_box[box_idx_ch] = -5;
      //                     }
      //                 }
      //             }
      //         }

      //         // Z axis
      //         if (ptr_ch->sim_bdry_contact_z == false && new_sim_bdry_contact_z == true)
      //         {
      //             for (int j = 0; j < ptr_ch->box_real_dim_y; j++)
      //             {
      //                 for (int i = 0; i < ptr_ch->box_real_dim_x; i++)
      //                 {
      //                     for (int k = 0; k < (ptr_ch->box_real_dim_z - ptr_ch->box_dim_z) / 2; k++)
      //                     {
      //                         box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
      //                         ptr_ch->ptr_box[box_idx_ch] = -5;
      //                     }

      //                     for (int k = (ptr_ch->box_real_dim_z + ptr_ch->box_dim_z) / 2; k < ptr_ch->box_real_dim_z; k++)
      //                     {
      //                         box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
      //                         ptr_ch->ptr_box[box_idx_ch] = -5;
      //                     }
      //                 }
      //             }
      //         }
      //     }
      // }

      //* >> Grid points *//
      ptr_ch->grid_intr_size = 0;
      ptr_ch->grid_bdry_size = 0;
      ptr_ch->grid_sim_bdry_size = 0;

      //* >> Reset grid properties: Density  *//
      size = (ptr_ch->box_real_dim_x + 1) * (ptr_ch->box_real_dim_y + 1) * (ptr_ch->box_real_dim_z + 1);
      for (int j = 0; j < size; j++)
      {
        // ptr_ch->ptr_pot[j] = 0;
        // ptr_ch->ptr_ax[j] = 0;
        // ptr_ch->ptr_ay[j] = 0;
        // ptr_ch->ptr_az[j] = 0;
        ptr_ch->ptr_d[j] = 0;
      }

      //* >> Boundary flags of the simulation box *//
      ptr_ch->sim_bdry_contact = new_sim_bdry_contact;
      ptr_ch->sim_bdry_contact_x = new_sim_bdry_contact_x;
      ptr_ch->sim_bdry_contact_y = new_sim_bdry_contact_y;
      ptr_ch->sim_bdry_contact_z = new_sim_bdry_contact_z;

      ptr_ch->pbc_crosses_sim_box_bdry = new_pbc_crosses_sim_box_bdry;
      ptr_ch->pbc_crosses_sim_box_bdry_x = new_pbc_crosses_sim_box_bdry_x; // when one node croses the simulation box at X direcction
      ptr_ch->pbc_crosses_sim_box_bdry_y = new_pbc_crosses_sim_box_bdry_y;
      ptr_ch->pbc_crosses_sim_box_bdry_z = new_pbc_crosses_sim_box_bdry_z;

      ptr_ch->pbc_crosses_whole_sim_box = new_pbc_crosses_whole_sim_box;
      ptr_ch->pbc_crosses_whole_sim_box_x = new_pbc_crosses_whole_sim_box_x;
      ptr_ch->pbc_crosses_whole_sim_box_y = new_pbc_crosses_whole_sim_box_y;
      ptr_ch->pbc_crosses_whole_sim_box_z = new_pbc_crosses_whole_sim_box_z;

      // printf("box_real_dim_x = %d\n",ptr_ch->box_real_dim_x);
      // printf("box_real_dim_y = %d\n",ptr_ch->box_real_dim_y);
      // printf("box_real_dim_z = %d\n",ptr_ch->box_real_dim_z);
      // printf("box_ts_x = %d\n", ptr_ch->box_ts_x);
      // printf("box_ts_y = %d\n", ptr_ch->box_ts_y);
      // printf("box_ts_z = %d\n\n", ptr_ch->box_ts_z);

      // printf("movement_box_x = %d\n", movement_box_x);
      // printf("movement_box_y = %d\n", movement_box_y);
      // printf("movement_box_z = %d\n", movement_box_z);
      // printf("corrections = %d\n", ptr_ch->pbc_correction_due_pbc_flag);
      // printf("ptr_ch->sim_bdry_contact = %d\n", ptr_ch->sim_bdry_contact);
      // printf("ptr_ch->sim_bdry_contact_x = %d\n", ptr_ch->sim_bdry_contact_x);
      // printf("ptr_ch->sim_bdry_contact_y = %d\n", ptr_ch->sim_bdry_contact_y);
      // printf("ptr_ch->sim_bdry_contact_z = %d\n", ptr_ch->sim_bdry_contact_z);
      // printf("ptr_ch->pbc_crosses_sim_box_bdry = %d\n", ptr_ch->pbc_crosses_sim_box_bdry);
      // printf("ptr_ch->pbc_crosses_sim_box_bdry_x = %d\n", ptr_ch->pbc_crosses_sim_box_bdry_x);
      // printf("ptr_ch->pbc_crosses_sim_box_bdry_y = %d\n", ptr_ch->pbc_crosses_sim_box_bdry_y);
      // printf("ptr_ch->pbc_crosses_sim_box_bdry_z = %d\n", ptr_ch->pbc_crosses_sim_box_bdry_z);
      // printf("ptr_ch->pbc_crosses_whole_sim_box = %d\n", ptr_ch->pbc_crosses_whole_sim_box);
      // printf("ptr_ch->pbc_crosses_whole_sim_box_x = %d\n", ptr_ch->pbc_crosses_whole_sim_box_x);
      // printf("ptr_ch->pbc_crosses_whole_sim_box_y = %d\n", ptr_ch->pbc_crosses_whole_sim_box_y);
      // printf("ptr_ch->pbc_crosses_whole_sim_box_z = %d\n", ptr_ch->pbc_crosses_whole_sim_box_z);
    }
  }

  return _SUCCESS_;
}

/**
 * @brief New child nodes are created to fit with their new refinement zones.
 *
 * **SHORT DESCRIPTION**: New child nodes are created to fit with refinement 
 * zones that were not linked to the old child nodes
 *
 * **PREREQUISITES**: If the number of new refinement zones is bigger than the
 * number of old child nodes.
 *
 * @param[in,out] ptr_node Pointer to node structure
 *
 * **RETURN**: The error status.
 *
 * **LONG DESCRIPTION**:
 *
 * In the case there are more new refinement zones than old child nodes, new
 * child nodes need to be created and prepared to fit the cells of the 
 * refinement zone. The parameters to be fitted are the same as the function 
 * tree_adaptation(), so at the end, we will obtain new child nodes exactly in 
 * the same position as the old adapted nodes.
 *
 * Unlike the function function tree_adaptation(), the fit box flag parameter
 * \ref node.box_check_fit "box_check_fit" will always be equal to "true".
 *
 * The flux of this function can be seen in the figure (work in progress), and
 * it is explained below:
 *
 * - [0]  <b> THE create_new_child_nodes() FUNCTION STARTS....</b>
 *
 * - [1]  Defining some internal useful parameters.
 *
 * - [2]  Run a "for" loop over the difference between new refinement zones and
 * old child nodes.
 *
 * - [3]  The initial child node is obtained through the external function
 * new_node(), which reuses already allocated nodes that have been reset.
 *
 * - [4]  If periodic boundary conditions (pbc) is activated, and if any 
 * boundary flag anomalies \ref node.ptr_pbc_bool_bdry_anomalies_x 
 * "ptr_pbc_bool_bdry_anomalies_x" (\ref node.ptr_pbc_bool_bdry_anomalies_y "y",
 * \ref node.ptr_pbc_bool_bdry_anomalies_z "z")  is activated, then the minimum
 * and maximum of the subzones is computed through the function
 * find_min_max_subzones_ref_PERIODIC_BOUNDARY().
 *
 * - [5]  For every dirección X, Y, and Z, the minimum and maximum, and the
 * boundary flags are computed and stored in local parameters.
 *
 * - [6]  The minimum and maximum of the previous point were computed in the
 * parent node, so then, they are transformed to the next level of refinement to
 * be accord with the respective child node.
 *
 * - [7]  The box dimensions parameters and traslation constants (see \ref 
 * node.box_ts_x "box_ts_x" (\ref node.box_ts_y "y", \ref node.box_ts_z "z") are 
 * computed.
 *
 * - [8]  The box array is initialized at "No-Exist" status (-4).
 *
 * - [9]  Finally, the child node is linked to the parent node through the
 * struct node pointers parameters \ref node.pptr_chn "pptr_chn" and \ref 
 * node.ptr_pt "ptr_pt".
 *
 * - [10] <b> THE create_new_child_nodes() FUNCTION ENDS....</b>
 *
 * **ILLUSTRATIVE EXAMPLES**:
 * - [a]  Trivial.
 *
 * **RATIONALES**:
 * - [a]  xxxx
 *
 * **NOTES**:
 * - [a]
 */

static int create_new_child_nodes(struct node *ptr_node)
{

  struct node *ptr_ch; // child node
  int cell_idx;        // The cell index is simply i of the for loop

  int pos_x; // Distance between the real box and the Smallest Box when the last is localized in the middle of the real one
  int pos_y;
  int pos_z;

  // int box_idx_ch;
  // int box_real_dim_X_ch;
  // int box_real_dim_X_times_Y_ch;

  int new_zone_idx; // ID of the new zone of refinement

  int lv = ptr_node->lv;

  bool check; // check subzone analized

  int aux_subzones_analized;

  int size; // Size or number of elements in some array in child nodes

  //* >> Space checking in the number of child nodes of ptr_node
  if (space_check(&(ptr_node->chn_cap), ptr_node->zones_size, 2.0f, "p1n2", &(ptr_node->pptr_chn)) == _FAILURE_)
  {
    printf("Error, in space_check function\n");
    return _FAILURE_;
  }

  for (int zone_idx = ptr_node->chn_size; zone_idx < ptr_node->zones_size; zone_idx++)
  {
    ptr_ch = new_node();

    new_zone_idx = ptr_node->ptr_links_new_ord_old[zone_idx];

    //* >> Global node properties *//
    // ptr_ch->ID = zone_idx;
    ptr_ch->ID = new_zone_idx;
    ptr_ch->lv = lv + 1;
    //* >> Cells in the node *//
    size = 8 * ptr_node->ptr_zone_size[new_zone_idx];

    //* >> Space checking of cells indices of the child node
    if (space_check(&(ptr_ch->cell_cap), size, 1.0f, "p4i1i1i1i1", &(ptr_ch->ptr_cell_idx_x), &(ptr_ch->ptr_cell_idx_y), &(ptr_ch->ptr_cell_idx_z), &(ptr_ch->ptr_box_idx)) == _FAILURE_)
    {
      printf("Error, in space_check function\n");
      return _FAILURE_;
    }

    // printf("\n Creating child, lv = %d, ID = %d\n", ptr_ch->lv,ptr_ch->ID);

    //* >> Boxes *//
    // MIN and MAX cell indices values of the node.
    if (bdry_cond_type == 0 && ptr_node->pbc_crosses_whole_sim_box == true &&
        (ptr_node->ptr_pbc_bool_bdry_anomalies_x[new_zone_idx] == true ||
         ptr_node->ptr_pbc_bool_bdry_anomalies_y[new_zone_idx] == true ||
         ptr_node->ptr_pbc_bool_bdry_anomalies_z[new_zone_idx] == true))
    {

      // printf("there are boundary anomalies\n");
      // printf("ptr_node->ptr_pbc_bool_bdry_anomalies_x[new_zone_idx] = %d\n", ptr_node->ptr_pbc_bool_bdry_anomalies_x[new_zone_idx]);
      // printf("ptr_node->ptr_pbc_bool_bdry_anomalies_y[new_zone_idx] = %d\n", ptr_node->ptr_pbc_bool_bdry_anomalies_y[new_zone_idx]);
      // printf("ptr_node->ptr_pbc_bool_bdry_anomalies_z[new_zone_idx] = %d\n", ptr_node->ptr_pbc_bool_bdry_anomalies_z[new_zone_idx]);

      // printf("Calling min and max subzones functionn\n");

      if (find_min_max_subzones_ref_PERIODIC_BOUNDARY(ptr_node, new_zone_idx) == _FAILURE_)
      {
        printf("Error at function fill_subzones_ref_PERIODIC_BOUNDARY()\n");
        return _FAILURE_;
      }

      // printf("ptr_node->pbc_subzones_size = %d\n", ptr_node->pbc_subzones_size);
    }

    //* >> X axis
    if (ptr_node->pbc_crosses_whole_sim_box_x == true)
    {
      if (ptr_node->ptr_pbc_bool_bdry_anomalies_x[new_zone_idx] == false)
      {
        for (int j = 0; j < ptr_node->ptr_zone_size[new_zone_idx]; j++)
        {
          cell_idx = ptr_node->pptr_zones[new_zone_idx][j]; // Cell index in the parent zone
          if (ptr_ch->box_min_x > ptr_node->ptr_cell_idx_x[cell_idx])
          {
            ptr_ch->box_min_x = ptr_node->ptr_cell_idx_x[cell_idx];
          }
          if (ptr_ch->box_max_x < ptr_node->ptr_cell_idx_x[cell_idx])
          {
            ptr_ch->box_max_x = ptr_node->ptr_cell_idx_x[cell_idx];
          }
        }
      }
      else if (ptr_node->pbc_subzones_size == 1)
      {
        ptr_ch->box_min_x = 0;
        ptr_ch->box_max_x = (1 << lv) - 1;
      }
      else
      {
        // Matching the subzones
        aux_subzones_analized = ptr_node->pbc_subzones_size;
        check = true;
        while (aux_subzones_analized > 1 && check == true)
        {
          check = false;
          for (int i = 0; i < 2; i++)
          {
            for (int j = i + 1; j < aux_subzones_analized; j++)
            {
              // if ((ptr_node->ptr_pbc_min_subzones_x[j] <= ptr_node->ptr_pbc_min_subzones_x[i] && ptr_node->ptr_pbc_max_subzones_x[j] >= ptr_node->ptr_pbc_min_subzones_x[i]) ||
              //     (ptr_node->ptr_pbc_max_subzones_x[j] >= ptr_node->ptr_pbc_max_subzones_x[i] && ptr_node->ptr_pbc_min_subzones_x[j] <= ptr_node->ptr_pbc_max_subzones_x[i]))
              if ((ptr_node->ptr_pbc_max_subzones_x[j] + 1 >= ptr_node->ptr_pbc_min_subzones_x[i] && ptr_node->ptr_pbc_min_subzones_x[j] <= 1 + ptr_node->ptr_pbc_max_subzones_x[i]))
              {
                ptr_node->ptr_pbc_min_subzones_x[i] = ptr_node->ptr_pbc_min_subzones_x[i] < ptr_node->ptr_pbc_min_subzones_x[j] ? ptr_node->ptr_pbc_min_subzones_x[i] : ptr_node->ptr_pbc_min_subzones_x[j];
                ptr_node->ptr_pbc_max_subzones_x[i] = ptr_node->ptr_pbc_max_subzones_x[i] > ptr_node->ptr_pbc_max_subzones_x[j] ? ptr_node->ptr_pbc_max_subzones_x[i] : ptr_node->ptr_pbc_max_subzones_x[j];

                // Reducing the total amount of analized subzones
                ptr_node->ptr_pbc_min_subzones_x[j] = ptr_node->ptr_pbc_min_subzones_x[aux_subzones_analized - 1];
                ptr_node->ptr_pbc_max_subzones_x[j] = ptr_node->ptr_pbc_max_subzones_x[aux_subzones_analized - 1];
                aux_subzones_analized--;
                j--;
                check = true;
              }
            }
          }
        }

        if (aux_subzones_analized > 2 || aux_subzones_analized < 1)
        {
          printf("error, aux_subzones analized is bigger than 2 or lower than 1, it is equal to = %d\n", aux_subzones_analized);
        }

        // min and max between all subzones
        if (aux_subzones_analized == 1)
        {
          // ptr_ch->box_min_x = ptr_node->ptr_pbc_min_subzones[0];
          // ptr_ch->box_max_x = ptr_node->ptr_pbc_max_subzones[0];
          ptr_ch->box_min_x = 0;
          ptr_ch->box_max_x = (1 << lv) - 1;
        }
        else
        {
          ptr_ch->box_min_x = ptr_node->ptr_pbc_min_subzones_x[0] == 0 ? ptr_node->ptr_pbc_min_subzones_x[1] - (1 << lv) : ptr_node->ptr_pbc_min_subzones_x[0] - (1 << lv);
          ptr_ch->box_max_x = ptr_node->ptr_pbc_min_subzones_x[0] == 0 ? ptr_node->ptr_pbc_max_subzones_x[0] : ptr_node->ptr_pbc_max_subzones_x[1];
        }
      }

      // Analysis of min and max of the refinement zone
      // Case refinement zone crosses the whole simulation box
      if (ptr_ch->box_max_x - ptr_ch->box_min_x == (1 << lv) - 1)
      {
        ptr_ch->sim_bdry_contact = true;
        ptr_ch->sim_bdry_contact_x = true;

        ptr_ch->pbc_crosses_sim_box_bdry = true;
        ptr_ch->pbc_crosses_sim_box_bdry_x = true;

        ptr_ch->pbc_crosses_whole_sim_box = true;
        ptr_ch->pbc_crosses_whole_sim_box_x = true;
      }
      // Case refinement zone crosses the simulation boundary
      else if (ptr_ch->box_min_x < 0)
      {
        ptr_ch->sim_bdry_contact = true;
        ptr_ch->sim_bdry_contact_x = true;

        ptr_ch->pbc_crosses_sim_box_bdry = true;
        ptr_ch->pbc_crosses_sim_box_bdry_x = true;
      }
      // Case refinement zone touches the simulation boundary
      else if (ptr_ch->box_min_x == 0 || ptr_ch->box_max_x == (1 << lv) - 1)
      {
        ptr_ch->sim_bdry_contact = true;
        ptr_ch->sim_bdry_contact_x = true;
      }

    } // End if(ptr_node->pbc_crosses_whole_sim_box_x == true)
    // Case Parent x axis crosses the simulation box
    else if (ptr_node->pbc_crosses_sim_box_bdry_x == true)
    {
      for (int j = 0; j < ptr_node->ptr_zone_size[new_zone_idx]; j++)
      {
        cell_idx = ptr_node->pptr_zones[new_zone_idx][j]; // Cell index in the parent zone

        if (ptr_node->ptr_cell_idx_x[cell_idx] > ptr_node->box_max_x)
        {
          ptr_ch->box_min_x = ptr_ch->box_min_x < ptr_node->ptr_cell_idx_x[cell_idx] - (1 << lv) ? ptr_ch->box_min_x : ptr_node->ptr_cell_idx_x[cell_idx] - (1 << lv);
          ptr_ch->box_max_x = ptr_ch->box_max_x > ptr_node->ptr_cell_idx_x[cell_idx] - (1 << lv) ? ptr_ch->box_max_x : ptr_node->ptr_cell_idx_x[cell_idx] - (1 << lv);
        }
        else
        {
          ptr_ch->box_min_x = ptr_ch->box_min_x < ptr_node->ptr_cell_idx_x[cell_idx] ? ptr_ch->box_min_x : ptr_node->ptr_cell_idx_x[cell_idx];
          ptr_ch->box_max_x = ptr_ch->box_max_x > ptr_node->ptr_cell_idx_x[cell_idx] ? ptr_ch->box_max_x : ptr_node->ptr_cell_idx_x[cell_idx];
        }
      }

      // Analysis of min and max of the refinement zone
      if (ptr_ch->box_min_x < 0 && ptr_ch->box_max_x >= 0)
      {
        ptr_ch->sim_bdry_contact = true;
        ptr_ch->sim_bdry_contact_x = true;

        ptr_ch->pbc_crosses_sim_box_bdry = true;
        ptr_ch->pbc_crosses_sim_box_bdry_x = true;
      }
      else if (ptr_ch->box_min_x == 0)
      {
        ptr_ch->sim_bdry_contact = true;
        ptr_ch->sim_bdry_contact_x = true;
      }
      else if (ptr_ch->box_max_x <= -1)
      {
        if (ptr_ch->box_max_x == -1)
        {
          ptr_ch->sim_bdry_contact = true;
          ptr_ch->sim_bdry_contact_x = true;
        }

        ptr_ch->box_min_x += (1 << lv);
        ptr_ch->box_max_x += (1 << lv);
      }
    } // End else if(ptr_node->pbc_crosses_sim_box_bdry = true)
    // Case Parent x axis do not cross the simulation box
    else
    {
      for (int j = 0; j < ptr_node->ptr_zone_size[new_zone_idx]; j++)
      {
        cell_idx = ptr_node->pptr_zones[new_zone_idx][j]; // Cell index in the parent zone
        if (ptr_ch->box_min_x > ptr_node->ptr_cell_idx_x[cell_idx])
        {
          ptr_ch->box_min_x = ptr_node->ptr_cell_idx_x[cell_idx];
        }
        if (ptr_ch->box_max_x < ptr_node->ptr_cell_idx_x[cell_idx])
        {
          ptr_ch->box_max_x = ptr_node->ptr_cell_idx_x[cell_idx];
        }
      }

      // Analysis of min and max of the refinement zone
      if (ptr_ch->box_min_x == 0 || ptr_ch->box_max_x == (1 << lv) - 1)
      {
        ptr_ch->sim_bdry_contact = true;
        ptr_ch->sim_bdry_contact_x = true;
      }
    }

    //* >> Y axis
    if (ptr_node->pbc_crosses_whole_sim_box_y == true)
    {
      if (ptr_node->ptr_pbc_bool_bdry_anomalies_y[new_zone_idx] == false)
      {
        for (int j = 0; j < ptr_node->ptr_zone_size[new_zone_idx]; j++)
        {
          cell_idx = ptr_node->pptr_zones[new_zone_idx][j]; // Cell index in the parent zone
          if (ptr_ch->box_min_y > ptr_node->ptr_cell_idx_y[cell_idx])
          {
            ptr_ch->box_min_y = ptr_node->ptr_cell_idx_y[cell_idx];
          }
          if (ptr_ch->box_max_y < ptr_node->ptr_cell_idx_y[cell_idx])
          {
            ptr_ch->box_max_y = ptr_node->ptr_cell_idx_y[cell_idx];
          }
        }
      }
      else if (ptr_node->pbc_subzones_size == 1)
      {
        ptr_ch->box_min_y = 0;
        ptr_ch->box_max_y = (1 << lv) - 1;
      }
      else
      {
        // Matching the subzones
        aux_subzones_analized = ptr_node->pbc_subzones_size;
        check = true;
        while (aux_subzones_analized > 1 && check == true)
        {
          check = false;
          for (int i = 0; i < 2; i++)
          {
            for (int j = i + 1; j < aux_subzones_analized; j++)
            {
              // if ((ptr_node->ptr_pbc_min_subzones_y[j] <= ptr_node->ptr_pbc_min_subzones_y[i] && ptr_node->ptr_pbc_max_subzones_y[j] >= ptr_node->ptr_pbc_min_subzones_y[i]) ||
              //     (ptr_node->ptr_pbc_max_subzones_y[j] >= ptr_node->ptr_pbc_max_subzones_y[i] && ptr_node->ptr_pbc_min_subzones_y[j] <= ptr_node->ptr_pbc_max_subzones_y[i]))
              if ((ptr_node->ptr_pbc_max_subzones_y[j] + 1 >= ptr_node->ptr_pbc_min_subzones_y[i] && ptr_node->ptr_pbc_min_subzones_y[j] <= 1 + ptr_node->ptr_pbc_max_subzones_y[i]))
              {
                ptr_node->ptr_pbc_min_subzones_y[i] = ptr_node->ptr_pbc_min_subzones_y[i] < ptr_node->ptr_pbc_min_subzones_y[j] ? ptr_node->ptr_pbc_min_subzones_y[i] : ptr_node->ptr_pbc_min_subzones_y[j];
                ptr_node->ptr_pbc_max_subzones_y[i] = ptr_node->ptr_pbc_max_subzones_y[i] > ptr_node->ptr_pbc_max_subzones_y[j] ? ptr_node->ptr_pbc_max_subzones_y[i] : ptr_node->ptr_pbc_max_subzones_y[j];

                // Reducing the total amount of analized subzones
                ptr_node->ptr_pbc_min_subzones_y[j] = ptr_node->ptr_pbc_min_subzones_y[aux_subzones_analized - 1];
                ptr_node->ptr_pbc_max_subzones_y[j] = ptr_node->ptr_pbc_max_subzones_y[aux_subzones_analized - 1];
                aux_subzones_analized--;
                j--;
                check = true;
              }
            }
          }
        }

        // min and max between all subzones
        if (aux_subzones_analized == 1)
        {
          ptr_ch->box_min_y = 0;
          ptr_ch->box_max_y = (1 << lv) - 1;
        }
        else
        {
          ptr_ch->box_min_y = ptr_node->ptr_pbc_min_subzones_y[0] == 0 ? ptr_node->ptr_pbc_min_subzones_y[1] - (1 << lv) : ptr_node->ptr_pbc_min_subzones_y[0] - (1 << lv);
          ptr_ch->box_max_y = ptr_node->ptr_pbc_min_subzones_y[0] == 0 ? ptr_node->ptr_pbc_max_subzones_y[0] : ptr_node->ptr_pbc_max_subzones_y[1];
        }
      }

      // Analysis of min and max of the refinement zone
      // Case refinement zone crosses the whole simulation box
      if (ptr_ch->box_max_y - ptr_ch->box_min_y == (1 << lv) - 1)
      {
        ptr_ch->sim_bdry_contact = true;
        ptr_ch->sim_bdry_contact_y = true;

        ptr_ch->pbc_crosses_sim_box_bdry = true;
        ptr_ch->pbc_crosses_sim_box_bdry_y = true;

        ptr_ch->pbc_crosses_whole_sim_box = true;
        ptr_ch->pbc_crosses_whole_sim_box_y = true;
      }
      // Case refinement zone crosses the simulation boundary
      else if (ptr_ch->box_min_y < 0)
      {
        ptr_ch->sim_bdry_contact = true;
        ptr_ch->sim_bdry_contact_y = true;

        ptr_ch->pbc_crosses_sim_box_bdry = true;
        ptr_ch->pbc_crosses_sim_box_bdry_y = true;
      }
      // Case refinement zone touches the simulation boundary
      else if (ptr_ch->box_min_y == 0 || ptr_ch->box_max_y == (1 << lv) - 1)
      {
        ptr_ch->sim_bdry_contact = true;
        ptr_ch->sim_bdry_contact_y = true;
      }

    } // End if(ptr_node->pbc_crosses_whole_sim_box_y == true)
    // Case Parent y axis crosses the simulation box
    else if (ptr_node->pbc_crosses_sim_box_bdry_y == true)
    {
      for (int j = 0; j < ptr_node->ptr_zone_size[new_zone_idx]; j++)
      {
        cell_idx = ptr_node->pptr_zones[new_zone_idx][j]; // Cell index in the parent zone

        if (ptr_node->ptr_cell_idx_y[cell_idx] > ptr_node->box_max_y)
        {
          ptr_ch->box_min_y = ptr_ch->box_min_y < ptr_node->ptr_cell_idx_y[cell_idx] - (1 << lv) ? ptr_ch->box_min_y : ptr_node->ptr_cell_idx_y[cell_idx] - (1 << lv);
          ptr_ch->box_max_y = ptr_ch->box_max_y > ptr_node->ptr_cell_idx_y[cell_idx] - (1 << lv) ? ptr_ch->box_max_y : ptr_node->ptr_cell_idx_y[cell_idx] - (1 << lv);
        }
        else
        {
          ptr_ch->box_min_y = ptr_ch->box_min_y < ptr_node->ptr_cell_idx_y[cell_idx] ? ptr_ch->box_min_y : ptr_node->ptr_cell_idx_y[cell_idx];
          ptr_ch->box_max_y = ptr_ch->box_max_y > ptr_node->ptr_cell_idx_y[cell_idx] ? ptr_ch->box_max_y : ptr_node->ptr_cell_idx_y[cell_idx];
        }
      }

      // Analysis of min and max of the refinement zone
      if (ptr_ch->box_min_y < 0 && ptr_ch->box_max_y >= 0)
      {
        ptr_ch->sim_bdry_contact = true;
        ptr_ch->sim_bdry_contact_y = true;

        ptr_ch->pbc_crosses_sim_box_bdry = true;
        ptr_ch->pbc_crosses_sim_box_bdry_y = true;
      }
      else if (ptr_ch->box_min_y == 0)
      {
        ptr_ch->sim_bdry_contact = true;
        ptr_ch->sim_bdry_contact_y = true;
      }
      else if (ptr_ch->box_max_y <= -1)
      {
        if (ptr_ch->box_max_y == -1)
        {
          ptr_ch->sim_bdry_contact = true;
          ptr_ch->sim_bdry_contact_y = true;
        }

        ptr_ch->box_min_y += (1 << lv);
        ptr_ch->box_max_y += (1 << lv);
      }
    } // End else if(ptr_node->pbc_crosses_sim_box_bdry = true)
    // Case Parent x axis do not cross the simulation box
    else
    {
      for (int j = 0; j < ptr_node->ptr_zone_size[new_zone_idx]; j++)
      {
        cell_idx = ptr_node->pptr_zones[new_zone_idx][j]; // Cell index in the parent zone
        if (ptr_ch->box_min_y > ptr_node->ptr_cell_idx_y[cell_idx])
        {
          ptr_ch->box_min_y = ptr_node->ptr_cell_idx_y[cell_idx];
        }
        if (ptr_ch->box_max_y < ptr_node->ptr_cell_idx_y[cell_idx])
        {
          ptr_ch->box_max_y = ptr_node->ptr_cell_idx_y[cell_idx];
        }
      }

      // Analysis of min and max of the refinement zone
      if (ptr_ch->box_min_y == 0 || ptr_ch->box_max_y == (1 << lv) - 1)
      {
        ptr_ch->sim_bdry_contact = true;
        ptr_ch->sim_bdry_contact_y = true;
      }
    }

    //* >> Z axis
    if (ptr_node->pbc_crosses_whole_sim_box_z == true)
    {
      if (ptr_node->ptr_pbc_bool_bdry_anomalies_z[new_zone_idx] == false)
      {
        for (int j = 0; j < ptr_node->ptr_zone_size[new_zone_idx]; j++)
        {
          cell_idx = ptr_node->pptr_zones[new_zone_idx][j]; // Cell index in the parent zone
          if (ptr_ch->box_min_z > ptr_node->ptr_cell_idx_z[cell_idx])
          {
            ptr_ch->box_min_z = ptr_node->ptr_cell_idx_z[cell_idx];
          }
          if (ptr_ch->box_max_z < ptr_node->ptr_cell_idx_z[cell_idx])
          {
            ptr_ch->box_max_z = ptr_node->ptr_cell_idx_z[cell_idx];
          }
        }
      }
      else if (ptr_node->pbc_subzones_size == 1)
      {
        ptr_ch->box_min_z = 0;
        ptr_ch->box_max_z = (1 << lv) - 1;
      }
      else
      {
        // Matching the subzones
        aux_subzones_analized = ptr_node->pbc_subzones_size;
        check = true;
        while (aux_subzones_analized > 1 && check == true)
        {
          check = false;
          for (int i = 0; i < 2; i++)
          {
            for (int j = i + 1; j < aux_subzones_analized; j++)
            {
              // if ((ptr_node->ptr_pbc_min_subzones_z[j] <= ptr_node->ptr_pbc_min_subzones_z[i] && ptr_node->ptr_pbc_max_subzones_z[j] >= ptr_node->ptr_pbc_min_subzones_z[i]) ||
              //     (ptr_node->ptr_pbc_max_subzones_z[j] >= ptr_node->ptr_pbc_max_subzones_z[i] && ptr_node->ptr_pbc_min_subzones_z[j] <= ptr_node->ptr_pbc_max_subzones_z[i]))
              if ((ptr_node->ptr_pbc_max_subzones_z[j] + 1 >= ptr_node->ptr_pbc_min_subzones_z[i] && ptr_node->ptr_pbc_min_subzones_z[j] <= 1 + ptr_node->ptr_pbc_max_subzones_z[i]))
              {
                ptr_node->ptr_pbc_min_subzones_z[i] = ptr_node->ptr_pbc_min_subzones_z[i] < ptr_node->ptr_pbc_min_subzones_z[j] ? ptr_node->ptr_pbc_min_subzones_z[i] : ptr_node->ptr_pbc_min_subzones_z[j];
                ptr_node->ptr_pbc_max_subzones_z[i] = ptr_node->ptr_pbc_max_subzones_z[i] > ptr_node->ptr_pbc_max_subzones_z[j] ? ptr_node->ptr_pbc_max_subzones_z[i] : ptr_node->ptr_pbc_max_subzones_z[j];

                // Reducing the total amount of analized subzones
                ptr_node->ptr_pbc_min_subzones_z[j] = ptr_node->ptr_pbc_min_subzones_z[aux_subzones_analized - 1];
                ptr_node->ptr_pbc_max_subzones_z[j] = ptr_node->ptr_pbc_max_subzones_z[aux_subzones_analized - 1];
                aux_subzones_analized--;
                j--;
                check = true;
              }
            }
          }
        }

        // min and max between all subzones
        if (aux_subzones_analized == 1)
        {
          ptr_ch->box_min_z = 0;
          ptr_ch->box_max_z = (1 << lv) - 1;
        }
        else
        {
          ptr_ch->box_min_z = ptr_node->ptr_pbc_min_subzones_z[0] == 0 ? ptr_node->ptr_pbc_min_subzones_z[1] - (1 << lv) : ptr_node->ptr_pbc_min_subzones_z[0] - (1 << lv);
          ptr_ch->box_max_z = ptr_node->ptr_pbc_min_subzones_z[0] == 0 ? ptr_node->ptr_pbc_max_subzones_z[0] : ptr_node->ptr_pbc_max_subzones_z[1];
        }
      }

      // Analysis of min and max of the refinement zone
      // Case refinement zone crosses the whole simulation box
      if (ptr_ch->box_max_z - ptr_ch->box_min_z == (1 << lv) - 1)
      {
        ptr_ch->sim_bdry_contact = true;
        ptr_ch->sim_bdry_contact_z = true;

        ptr_ch->pbc_crosses_sim_box_bdry = true;
        ptr_ch->pbc_crosses_sim_box_bdry_z = true;

        ptr_ch->pbc_crosses_whole_sim_box = true;
        ptr_ch->pbc_crosses_whole_sim_box_z = true;
      }
      // Case refinement zone crosses the simulation boundary
      else if (ptr_ch->box_min_z < 0)
      {
        ptr_ch->sim_bdry_contact = true;
        ptr_ch->sim_bdry_contact_z = true;

        ptr_ch->pbc_crosses_sim_box_bdry = true;
        ptr_ch->pbc_crosses_sim_box_bdry_z = true;
      }
      // Case refinement zone touches the simulation boundary
      else if (ptr_ch->box_min_z == 0 || ptr_ch->box_max_z == (1 << lv) - 1)
      {
        ptr_ch->sim_bdry_contact = true;
        ptr_ch->sim_bdry_contact_z = true;
      }

    } // End if(ptr_node->pbc_crosses_whole_sim_box_z == true)
    // Case Parent x axis crosses the simulation box
    else if (ptr_node->pbc_crosses_sim_box_bdry_z == true)
    {
      for (int j = 0; j < ptr_node->ptr_zone_size[new_zone_idx]; j++)
      {
        cell_idx = ptr_node->pptr_zones[new_zone_idx][j]; // Cell index in the parent zone

        if (ptr_node->ptr_cell_idx_z[cell_idx] > ptr_node->box_max_z)
        {
          ptr_ch->box_min_z = ptr_ch->box_min_z < ptr_node->ptr_cell_idx_z[cell_idx] - (1 << lv) ? ptr_ch->box_min_z : ptr_node->ptr_cell_idx_z[cell_idx] - (1 << lv);
          ptr_ch->box_max_z = ptr_ch->box_max_z > ptr_node->ptr_cell_idx_z[cell_idx] - (1 << lv) ? ptr_ch->box_max_z : ptr_node->ptr_cell_idx_z[cell_idx] - (1 << lv);
        }
        else
        {
          ptr_ch->box_min_z = ptr_ch->box_min_z < ptr_node->ptr_cell_idx_z[cell_idx] ? ptr_ch->box_min_z : ptr_node->ptr_cell_idx_z[cell_idx];
          ptr_ch->box_max_z = ptr_ch->box_max_z > ptr_node->ptr_cell_idx_z[cell_idx] ? ptr_ch->box_max_z : ptr_node->ptr_cell_idx_z[cell_idx];
        }
      }

      // Analysis of min and max of the refinement zone
      if (ptr_ch->box_min_z < 0 && ptr_ch->box_max_z >= 0)
      {
        ptr_ch->sim_bdry_contact = true;
        ptr_ch->sim_bdry_contact_z = true;

        ptr_ch->pbc_crosses_sim_box_bdry = true;
        ptr_ch->pbc_crosses_sim_box_bdry_z = true;
      }
      else if (ptr_ch->box_min_z == 0)
      {
        ptr_ch->sim_bdry_contact = true;
        ptr_ch->sim_bdry_contact_z = true;
      }
      else if (ptr_ch->box_max_z <= -1)
      {
        if (ptr_ch->box_max_z == -1)
        {
          ptr_ch->sim_bdry_contact = true;
          ptr_ch->sim_bdry_contact_z = true;
        }

        ptr_ch->box_min_z += (1 << lv);
        ptr_ch->box_max_z += (1 << lv);
      }
    } // End else if(ptr_node->pbc_crosses_sim_box_bdry = true)
    // Case Parent x axis do not cross the simulation box
    else
    {
      for (int j = 0; j < ptr_node->ptr_zone_size[new_zone_idx]; j++)
      {
        cell_idx = ptr_node->pptr_zones[new_zone_idx][j]; // Cell index in the parent zone
        if (ptr_ch->box_min_z > ptr_node->ptr_cell_idx_z[cell_idx])
        {
          ptr_ch->box_min_z = ptr_node->ptr_cell_idx_z[cell_idx];
        }
        if (ptr_ch->box_max_z < ptr_node->ptr_cell_idx_z[cell_idx])
        {
          ptr_ch->box_max_z = ptr_node->ptr_cell_idx_z[cell_idx];
        }
      }

      // Analysis of min and max of the refinement zone
      if (ptr_ch->box_min_z == 0 || ptr_ch->box_max_z == (1 << lv) - 1)
      {
        ptr_ch->sim_bdry_contact = true;
        ptr_ch->sim_bdry_contact_z = true;
      }
    }

    // Changing the min and max of the "Smallest Box" from parent units to child units
    ptr_ch->box_min_x = 2 * ptr_ch->box_min_x;
    ptr_ch->box_min_y = 2 * ptr_ch->box_min_y;
    ptr_ch->box_min_z = 2 * ptr_ch->box_min_z;
    ptr_ch->box_max_x = 2 * ptr_ch->box_max_x + 1;
    ptr_ch->box_max_y = 2 * ptr_ch->box_max_y + 1;
    ptr_ch->box_max_z = 2 * ptr_ch->box_max_z + 1;

    // Size of the "Smallest Box"
    ptr_ch->box_dim_x = ptr_ch->box_max_x - ptr_ch->box_min_x + 1;
    ptr_ch->box_dim_y = ptr_ch->box_max_y - ptr_ch->box_min_y + 1;
    ptr_ch->box_dim_z = ptr_ch->box_max_z - ptr_ch->box_min_z + 1;

    // Real dimensions of the boxcap = ptr_ch->box_cap;
    ptr_ch->box_real_dim_x = 5 > (n_exp - 1) ? (ptr_ch->box_dim_x + 10) : (ptr_ch->box_dim_x + 2 * n_exp - 2);
    ptr_ch->box_real_dim_y = 5 > (n_exp - 1) ? (ptr_ch->box_dim_y + 10) : (ptr_ch->box_dim_y + 2 * n_exp - 2);
    ptr_ch->box_real_dim_z = 5 > (n_exp - 1) ? (ptr_ch->box_dim_z + 10) : (ptr_ch->box_dim_z + 2 * n_exp - 2);

    // ptr_ch->box_real_dim_x_old = ptr_ch->box_real_dim_x;
    // ptr_ch->box_real_dim_y_old = ptr_ch->box_real_dim_y;
    // ptr_ch->box_real_dim_z_old = ptr_ch->box_real_dim_z;

    // Translations between cell array and box
    pos_x = (ptr_ch->box_real_dim_x - ptr_ch->box_dim_x) / 2; // Half of the distance of the box side less the "Smallest Box" side
    ptr_ch->box_ts_x = ptr_ch->box_min_x - pos_x;             // Every cell in the level l in the box must be subtracted this value to obtain the box index
    pos_y = (ptr_ch->box_real_dim_y - ptr_ch->box_dim_y) / 2;
    ptr_ch->box_ts_y = ptr_ch->box_min_y - pos_y;
    pos_z = (ptr_ch->box_real_dim_z - ptr_ch->box_dim_z) / 2;
    ptr_ch->box_ts_z = ptr_ch->box_min_z - pos_z;

    // Filling the boxes
    size = ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y * ptr_ch->box_real_dim_z;
    //* >> Space checking of boxes: box, box_aux and box_mass_aux
    if (space_check(&(ptr_ch->box_cap), size, 1.0f, "p4i1i1c1c1", &(ptr_ch->ptr_box), &(ptr_ch->ptr_box_old), &(ptr_ch->ptr_cell_struct), &(ptr_ch->ptr_cell_struct_old)) == _FAILURE_)
    {
      printf("Error, in space_check function\n");
      return _FAILURE_;
    }

    // Putting the value of NO-EXIST (-4) in every box index
    for (int j = 0; j < size; j++)
    {
      ptr_ch->ptr_box[j] = -4;
    }

    // // Adding -6 to the boundary when it corresponds
    // if (bdry_cond_type == 0 && ptr_ch->pbc_crosses_whole_sim_box == true)
    // {

    //     box_real_dim_X_ch = ptr_ch->box_real_dim_x;
    //     box_real_dim_X_times_Y_ch = ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;

    //     if (ptr_ch->pbc_crosses_whole_sim_box_x == true)
    //     {
    //         for (int k = 0; k < ptr_ch->box_real_dim_z; k++)
    //         {
    //             for (int j = 0; j < ptr_ch->box_real_dim_y; j++)
    //             {
    //                 for (int i = 0; i < (ptr_ch->box_real_dim_x - ptr_ch->box_dim_x) / 2; i++)
    //                 {
    //                     box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
    //                     ptr_ch->ptr_box[box_idx_ch] = -6;
    //                 }

    //                 for (int i = (ptr_ch->box_real_dim_x + ptr_ch->box_dim_x) / 2; i < ptr_ch->box_real_dim_x; i++)
    //                 {
    //                     box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
    //                     ptr_ch->ptr_box[box_idx_ch] = -6;
    //                 }
    //             }
    //         }
    //     }

    //     if (ptr_ch->pbc_crosses_whole_sim_box_y == true)
    //     {
    //         for (int k = 0; k < ptr_ch->box_real_dim_z; k++)
    //         {
    //             for (int i = 0; i < ptr_ch->box_real_dim_x; i++)
    //             {
    //                 for (int j = 0; j < (ptr_ch->box_real_dim_y - ptr_ch->box_dim_y) / 2; j++)
    //                 {
    //                     box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
    //                     ptr_ch->ptr_box[box_idx_ch] = -6;
    //                 }

    //                 for (int j = (ptr_ch->box_real_dim_y + ptr_ch->box_dim_y) / 2; j < ptr_ch->box_real_dim_y; j++)
    //                 {
    //                     box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
    //                     ptr_ch->ptr_box[box_idx_ch] = -6;
    //                 }
    //             }
    //         }
    //     }

    //     if (ptr_ch->pbc_crosses_whole_sim_box_z == true)
    //     {
    //         for (int j = 0; j < ptr_ch->box_real_dim_y; j++)
    //         {
    //             for (int i = 0; i < ptr_ch->box_real_dim_x; i++)
    //             {
    //                 for (int k = 0; k < (ptr_ch->box_real_dim_z - ptr_ch->box_dim_z) / 2; k++)
    //                 {
    //                     box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
    //                     ptr_ch->ptr_box[box_idx_ch] = -6;
    //                 }

    //                 for (int k = (ptr_ch->box_real_dim_z + ptr_ch->box_dim_z) / 2; k < ptr_ch->box_real_dim_z; k++)
    //                 {
    //                     box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
    //                     ptr_ch->ptr_box[box_idx_ch] = -6;
    //                 }
    //             }
    //         }
    //     }
    // }
    // // Adding -5 to the boundary when it corresponds
    // else if (bdry_cond_type != 0 && ptr_ch->sim_bdry_contact == true)
    // {

    //     box_real_dim_X_ch = ptr_ch->box_real_dim_x;
    //     box_real_dim_X_times_Y_ch = ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;

    //     if (ptr_ch->sim_bdry_contact_x == true)
    //     {
    //         for (int k = 0; k < ptr_ch->box_real_dim_z; k++)
    //         {
    //             for (int j = 0; j < ptr_ch->box_real_dim_y; j++)
    //             {
    //                 for (int i = 0; i < (ptr_ch->box_real_dim_x - ptr_ch->box_dim_x) / 2; i++)
    //                 {
    //                     box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
    //                     ptr_ch->ptr_box[box_idx_ch] = -5;
    //                 }

    //                 for (int i = (ptr_ch->box_real_dim_x + ptr_ch->box_dim_x) / 2; i < ptr_ch->box_real_dim_x; i++)
    //                 {
    //                     box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
    //                     ptr_ch->ptr_box[box_idx_ch] = -5;
    //                 }
    //             }
    //         }
    //     }

    //     if (ptr_ch->sim_bdry_contact_y == true)
    //     {
    //         for (int k = 0; k < ptr_ch->box_real_dim_z; k++)
    //         {
    //             for (int i = 0; i < ptr_ch->box_real_dim_x; i++)
    //             {
    //                 for (int j = 0; j < (ptr_ch->box_real_dim_y - ptr_ch->box_dim_y) / 2; j++)
    //                 {
    //                     box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
    //                     ptr_ch->ptr_box[box_idx_ch] = -5;
    //                 }

    //                 for (int j = (ptr_ch->box_real_dim_y + ptr_ch->box_dim_y) / 2; j < ptr_ch->box_real_dim_y; j++)
    //                 {
    //                     box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
    //                     ptr_ch->ptr_box[box_idx_ch] = -5;
    //                 }
    //             }
    //         }
    //     }

    //     if (ptr_ch->sim_bdry_contact_z == true)
    //     {
    //         for (int j = 0; j < ptr_ch->box_real_dim_y; j++)
    //         {
    //             for (int i = 0; i < ptr_ch->box_real_dim_x; i++)
    //             {
    //                 for (int k = 0; k < (ptr_ch->box_real_dim_z - ptr_ch->box_dim_z) / 2; k++)
    //                 {
    //                     box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
    //                     ptr_ch->ptr_box[box_idx_ch] = -5;
    //                 }

    //                 for (int k = (ptr_ch->box_real_dim_z + ptr_ch->box_dim_z) / 2; k < ptr_ch->box_real_dim_z; k++)
    //                 {
    //                     box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
    //                     ptr_ch->ptr_box[box_idx_ch] = -5;
    //                 }
    //             }
    //         }
    //     }
    // }

    //* >> Potential, acceleration and density of the grid *//
    size = (ptr_ch->box_real_dim_x + 1) * (ptr_ch->box_real_dim_y + 1) * (ptr_ch->box_real_dim_z + 1);
    //* >> Space checking
    if (space_check(&(ptr_ch->grid_properties_cap), size, 1.0f, "p6v1v1v1v1v1v1", &(ptr_ch->ptr_pot), &(ptr_ch->ptr_pot_old), &(ptr_ch->ptr_ax), &(ptr_ch->ptr_ay), &(ptr_ch->ptr_az), &(ptr_ch->ptr_d)) == _FAILURE_)
    {
      printf("Error, in space_check function\n");
      return _FAILURE_;
    }

    //* >> Tree structure *//
    ptr_ch->ptr_pt = ptr_node;
    ptr_node->pptr_chn[zone_idx] = ptr_ch;

    // printf("check fit = %d\n", ptr_ch->box_check_fit);
    // printf("box_min_x = %d, box_max_x = %d\n", ptr_ch->box_min_x, ptr_ch->box_max_x);
    // printf("box_min_y = %d, box_max_y = %d\n", ptr_ch->box_min_y, ptr_ch->box_max_y);
    // printf("box_min_z = %d, box_max_z = %d\n\n", ptr_ch->box_min_z, ptr_ch->box_max_z);
    // printf("box_real_dim_x = %d\n", ptr_ch->box_real_dim_x);
    // printf("box_real_dim_y = %d\n", ptr_ch->box_real_dim_y);
    // printf("box_real_dim_z = %d\n", ptr_ch->box_real_dim_z);
    // printf("box_ts_x = %d\n", ptr_ch->box_ts_x);
    // printf("box_ts_y = %d\n", ptr_ch->box_ts_y);
    // printf("box_ts_z = %d\n", ptr_ch->box_ts_z);
    // printf("corrections = %d\n", ptr_ch->pbc_correction_due_pbc_flag);
    // printf("ptr_ch->sim_bdry_contact = %d\n", ptr_ch->sim_bdry_contact);
    // printf("ptr_ch->sim_bdry_contact_x = %d\n", ptr_ch->sim_bdry_contact_x);
    // printf("ptr_ch->sim_bdry_contact_y = %d\n", ptr_ch->sim_bdry_contact_y);
    // printf("ptr_ch->sim_bdry_contact_z = %d\n", ptr_ch->sim_bdry_contact_z);
    // printf("ptr_ch->pbc_crosses_sim_box_bdry = %d\n", ptr_ch->pbc_crosses_sim_box_bdry);
    // printf("ptr_ch->pbc_crosses_sim_box_bdry_x = %d\n", ptr_ch->pbc_crosses_sim_box_bdry_x);
    // printf("ptr_ch->pbc_crosses_sim_box_bdry_y = %d\n", ptr_ch->pbc_crosses_sim_box_bdry_y);
    // printf("ptr_ch->pbc_crosses_sim_box_bdry_z = %d\n", ptr_ch->pbc_crosses_sim_box_bdry_z);
    // printf("ptr_ch->pbc_crosses_whole_sim_box = %d\n", ptr_ch->pbc_crosses_whole_sim_box);
    // printf("ptr_ch->pbc_crosses_whole_sim_box_x = %d\n", ptr_ch->pbc_crosses_whole_sim_box_x);
    // printf("ptr_ch->pbc_crosses_whole_sim_box_y = %d\n", ptr_ch->pbc_crosses_whole_sim_box_y);
    // printf("ptr_ch->pbc_crosses_whole_sim_box_z = %d\n", ptr_ch->pbc_crosses_whole_sim_box_z);
  }
  return _SUCCESS_;
}

/**
 * @brief Transferring refined parent cells from the old to the new child nodes.
 *
 * **SHORT DESCRIPTION**: Transferring cells that require be refined from the 
 * old child nodes to which they are to the new child nodes to which they will 
 * belong.
 *
 * **PREREQUISITES**: If there is at least one new and one old refinement zone.
 *
 * @param[in,out] ptr_node Pointer to node structure
 *
 * **RETURN**: The error status.
 *
 * \anchor tree_acaptaion__create_new_child_nodes__LONG_DESCRIPTION **LONG DESCRIPTION**:
 *
 * Transferring cells that require be refined from the old child nodes to which 
 * they are to the new child nodes to which they will belong. Furthermore, the 
 * new child box status is also updated with the new information of child nodes, 
 * thus also like the mass and number of particles of the child nodes and other 
 * child node parameters.
 *
 * Moreover, when a child node has its fit box flag parameter \ref
 * node.box_check_fit "box_check_fit" equal to "false", it is also necessary to
 * move cells from its old cell structure \ref node.ptr_cell_struct_old 
 * "ptr_cell_struct_old" to the new cell structure \ref node.ptr_cell_struct 
 * "ptr_cell_struct" if the cell still belongs to that child node.
 *
 * This function is separated into 2 big parts. The first one consists in to 
 * move cells between couples of child nodes which represent a new refinement 
 * zone. The second part consists in to move cells between couples of child 
 * nodes, but only the receiving cell represents a new refinement zone, the 
 * sender child node finished this function will be removed from the parent node 
 * and sent to the memory pool array to be used in the future by another parent 
 * node in any level of refinement.
 *
 * The flux of this function can be seen in the figure (work in progress), and
 * it is explained below:
 *
 * - [0]  <b> THE transferring_old_child_to_new_child() FUNCTION STARTS....</b>
 *
 * - [1]  Defining some internal useful parameters.
 *
 * - [2]  As explained above in the \ref 
 * tree_acaptaion__create_new_child_nodes__LONG_DESCRIPTION "LONG DESCRIPTION", 
 * the first part of the function consists in to move cells between future child 
 * nodes. To do that, we run a "while" loop over every old child node with is 
 * linked with a new refinement zone; a *zone_idx* counter is used to perform 
 * that loop.
 *
 * - [3]  If the sender child node (A) has its box fit paramter \ref 
 * node.box_check_fit "box_check_fit" equal to "true", then we proceed.
 *
 * - [4]  Run a "for" loop over the cells of the child node A. The loop runs
 * each 8 cells which correspond to one new refined cell in the parent node.
 * Remember that all cells that did not require refinement have already been 
 * removed through the function remov_cells_nolonger_require_refinement().
 *
 * - [5]  The parent box index associated with this package of 8 child cells is
 * computed.
 *
 * - [6]  In case the cell package will stay in the same child node no actions
 * will be taken except in the exceptional case that periodic boundary 
 * conditions (pbc) and the child node flag parameter \ref
 * node.pbc_correction_due_pbc_flag "pbc_correction_due_pbc_flag" are activated.
 * In this exceptional case could be necessary to update the cell box status 
 * and, the cell structure.
 *
 * - [7]  In case the other situation occurs in step [6], where the cell package
 * requires to be moved to another child node, then the cell box indices of the
 * child node A and the receiving child node (B) are computed, the cells are
 * transferred (from A to B), and the information of the child nodes is updated.
 *
 * - [8]  At this point, the "if statement" of step [3] is finalized. So, now
 * the other option is addressed when the sender child node (A) has its box fit
 * parameter \ref node.box_check_fit "box_check_fit" equal to "false".
 *
 * - [9]  In this case, we repeat almost the same idea as the "true" case, but 
 * it is necessary to use the auxiliary old cell structure of the child node A 
 * \ref node.ptr_cell_struct_old "ptr_cell_struct_old" instead of the new one 
 * \ref node.ptr_cell_struct "ptr_cell_struct". Moreover, in the case where the 
 * cell package will stay in the same child node always will require to be 
 * updated from the old cell structure to the new cell structure of the child
 * node A.
 *
 * - [10] At this point, the "while" loop of step [2] is finalized. So, now the
 * second part of the function is initialized. Here cells are moved between
 * couples of child nodes but the sender child node (A) will be removed from the
 * parent node to later be transferred to the memory pool to be reused in the 
 * future. The process is the same as the first step but easier. Here it is not
 * necessary to be worried about updating the information of the sender child
 * node because it will be removed from the parent node.
 *
 * - [11] <b> THE transferring_old_child_to_new_child() FUNCTION ENDS....</b>
 *
 * **ILLUSTRATIVE EXAMPLES**:
 * - [a]  Trivial.
 *
 * **RATIONALES**:
 * - [a]  The design of the updated is almost the same in the parts of this 
 * function with a few different actions depending on a few flags parameters. 
 * The reason why we decided to separate this process into several parts instead 
 * of joining all together and putting the "if" statements very inside of the 
 * loops is because we are speeding up the performance of the function in 
 * exchange for having a longer script.
 *
 * **NOTES**:
 * - [a]
 */

static int transferring_old_child_to_new_child(struct node *ptr_node)
{
  struct node *ptr_ch_A;
  struct node *ptr_ch_B;

  struct cell_struct cell_struct_aux;

  int no_cells_ch_A; // Total number of cells in the child node A
  int no_cells_ch_B; // Total number of cells in the child node A

  int cntr_ch; // Counter the number of children

  int new_zone_idx; // ID of the new zone of refinement

  int box_idxNbr_ch_A_old;

  int box_idx_x_node; // Box index in X direcction of the node cell
  int box_idx_y_node; // Box index in Y direcction of the node cell
  int box_idx_z_node; // Box index in Z direcction of the node cell
  int box_idx_node;   // Box index of the node cell

  int box_idx_x_ch_A; // Box index in X direcction of the child cell
  int box_idx_y_ch_A; // Box index in Y direcction of the child cell
  int box_idx_z_ch_A; // Box index in Z direcction of the child cell
  int box_idx_x_ch_B; // Box index in X direcction of the child cell
  int box_idx_y_ch_B; // Box index in Y direcction of the child cell
  int box_idx_z_ch_B; // Box index in Z direcction of the child cell
  int box_idx_ch_A;   // Box index of the child cell
  int box_idx_ch_B;   // Box index of the child cell

  int box_idxNbr_ch_A;
  int box_idxNbr_ch_B;

  int aux_int;
  int aux_int_A;
  int aux_int_B;
  int aux_int_X;

  int box_real_dim_X_node = ptr_node->box_real_dim_x;
  int box_real_dim_X_times_Y_node = ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;

  int box_real_dim_X_ch_A;
  int box_real_dim_X_times_Y_ch_A;

  int box_real_dim_X_ch_B;
  int box_real_dim_X_times_Y_ch_B;

  int aux_idx_x;
  int aux_idx_y;
  int aux_idx_z;

  int lv = ptr_node->lv;

  // Loop over new refinement zones,
  // Case of old child nodes to be reused

  int zone_idx = 0;
  // while (zone_idx < ptr_node->zones_size && ptr_node->ptr_links_old_ord_old[zone_idx] < ptr_node->chn_size)
  while (zone_idx < ptr_node->zones_size && zone_idx < ptr_node->chn_size)
  {
    new_zone_idx = ptr_node->ptr_links_new_ord_old[zone_idx];

    ptr_ch_A = ptr_node->pptr_chn[ptr_node->ptr_links_old_ord_old[zone_idx]];
    no_cells_ch_A = ptr_ch_A->cell_size;

    if (ptr_ch_A->box_check_fit == true)
    {
      for (int cell_idx = 0; cell_idx < no_cells_ch_A; cell_idx += 8)
      {

        aux_idx_x = ptr_ch_A->ptr_cell_idx_x[cell_idx] >> 1;
        aux_idx_y = ptr_ch_A->ptr_cell_idx_y[cell_idx] >> 1;
        aux_idx_z = ptr_ch_A->ptr_cell_idx_z[cell_idx] >> 1;

        box_idx_x_node = aux_idx_x - ptr_node->box_ts_x;
        box_idx_y_node = aux_idx_y - ptr_node->box_ts_y;
        box_idx_z_node = aux_idx_z - ptr_node->box_ts_z;

        if (ptr_node->pbc_crosses_sim_box_bdry == true)
        {
          if (aux_idx_x > ptr_node->box_max_x)
          {
            box_idx_x_node -= (1 << lv);
          }

          if (aux_idx_y > ptr_node->box_max_y)
          {
            box_idx_y_node -= (1 << lv);
          }

          if (aux_idx_z > ptr_node->box_max_z)
          {
            box_idx_z_node -= (1 << lv);
          }
        }

        box_idx_node = box_idx_x_node + box_idx_y_node * box_real_dim_X_node + box_idx_z_node * box_real_dim_X_times_Y_node;

        if (ptr_node->ptr_box[box_idx_node] != new_zone_idx)
        {
          box_idx_ch_A = ptr_ch_A->ptr_box_idx[cell_idx];

          //  ptr_ch_B: Is the child node where the information from the child A need to be sent
          ptr_ch_B = ptr_node->pptr_chn[ptr_node->ptr_links_old_ord_new[ptr_node->ptr_box[box_idx_node]]];

          //* >> CELLS *//
          no_cells_ch_B = ptr_ch_B->cell_size;

          box_real_dim_X_ch_B = ptr_ch_B->box_real_dim_x;
          box_real_dim_X_times_Y_ch_B = ptr_ch_B->box_real_dim_x * ptr_ch_B->box_real_dim_y;

          aux_idx_x = ptr_ch_A->ptr_cell_idx_x[cell_idx];
          aux_idx_y = ptr_ch_A->ptr_cell_idx_y[cell_idx];
          aux_idx_z = ptr_ch_A->ptr_cell_idx_z[cell_idx];

          box_idx_x_ch_B = aux_idx_x - ptr_ch_B->box_ts_x;
          box_idx_y_ch_B = aux_idx_y - ptr_ch_B->box_ts_y;
          box_idx_z_ch_B = aux_idx_z - ptr_ch_B->box_ts_z;

          if (ptr_ch_B->pbc_crosses_sim_box_bdry == true)
          {
            if (aux_idx_x > ptr_ch_B->box_max_x)
            {
              box_idx_x_ch_B -= (1 << (lv + 1));
            }

            if (aux_idx_y > ptr_ch_B->box_max_y)
            {
              box_idx_y_ch_B -= (1 << (lv + 1));
            }

            if (aux_idx_z > ptr_ch_B->box_max_z)
            {
              box_idx_z_ch_B -= (1 << (lv + 1));
            }
          }

          box_idx_ch_B = box_idx_x_ch_B + box_idx_y_ch_B * box_real_dim_X_ch_B + box_idx_z_ch_B * box_real_dim_X_times_Y_ch_B;

          // box_idx_x_ch_B = ptr_ch_A->ptr_cell_idx_x[cell_idx] - ptr_ch_B->box_ts_x;
          // box_idx_y_ch_B = ptr_ch_A->ptr_cell_idx_y[cell_idx] - ptr_ch_B->box_ts_y;
          // box_idx_z_ch_B = ptr_ch_A->ptr_cell_idx_z[cell_idx] - ptr_ch_B->box_ts_z;
          // box_idx_ch_B = box_idx_x_ch_B + box_idx_y_ch_B * box_real_dim_X_ch_B + box_idx_z_ch_B * box_real_dim_X_times_Y_ch_B;

          for (int kk = 0; kk < 2; kk++)
          {
            for (int jj = 0; jj < 2; jj++)
            {
              for (int ii = 0; ii < 2; ii++)
              {
                aux_int = ii + jj * 2 + kk * 4;
                aux_int_A = cell_idx + aux_int;
                aux_int_B = no_cells_ch_B + aux_int;
                aux_int_X = no_cells_ch_A - 8 + aux_int;
                // Transferring from child A to child B
                ptr_ch_B->ptr_cell_idx_x[aux_int_B] = ptr_ch_A->ptr_cell_idx_x[aux_int_A];
                ptr_ch_B->ptr_cell_idx_y[aux_int_B] = ptr_ch_A->ptr_cell_idx_y[aux_int_A];
                ptr_ch_B->ptr_cell_idx_z[aux_int_B] = ptr_ch_A->ptr_cell_idx_z[aux_int_A];
                ptr_ch_B->ptr_box_idx[aux_int_B] = box_idx_ch_B + ii + jj * box_real_dim_X_ch_B + kk * box_real_dim_X_times_Y_ch_B;

                // Updating the box status of child node A and B
                box_idxNbr_ch_A = ptr_ch_A->ptr_box_idx[aux_int_A];
                box_idxNbr_ch_B = ptr_ch_B->ptr_box_idx[aux_int_B];
                ptr_ch_A->ptr_box[box_idxNbr_ch_A] = -4; // Putting the status of NO EXIST (-4) in the child node cell
                ptr_ch_B->ptr_box[box_idxNbr_ch_B] = -3; // Putting the status of EXIST (-3) in the child node cell

                //* >> Cell structure and local mass *//
                cell_struct_aux = ptr_ch_B->ptr_cell_struct[box_idxNbr_ch_B];
                ptr_ch_B->ptr_cell_struct[box_idxNbr_ch_B] = ptr_ch_A->ptr_cell_struct[box_idxNbr_ch_A];
                ptr_ch_A->ptr_cell_struct[box_idxNbr_ch_A] = cell_struct_aux;

                // Local mass
                ptr_ch_B->node_mass += ptr_ch_B->ptr_cell_struct[box_idxNbr_ch_B].cell_mass;
                ptr_ch_B->no_ptcl_full_node += ptr_ch_B->ptr_cell_struct[box_idxNbr_ch_B].ptcl_size;

                ptr_ch_A->node_mass -= ptr_ch_B->ptr_cell_struct[box_idxNbr_ch_B].cell_mass;
                ptr_ch_A->no_ptcl_full_node -= ptr_ch_B->ptr_cell_struct[box_idxNbr_ch_B].ptcl_size;

                // Transferring from end of child A to removed cells of child A
                ptr_ch_A->ptr_cell_idx_x[aux_int_A] = ptr_ch_A->ptr_cell_idx_x[aux_int_X];
                ptr_ch_A->ptr_cell_idx_y[aux_int_A] = ptr_ch_A->ptr_cell_idx_y[aux_int_X];
                ptr_ch_A->ptr_cell_idx_z[aux_int_A] = ptr_ch_A->ptr_cell_idx_z[aux_int_X];
                ptr_ch_A->ptr_box_idx[aux_int_A] = ptr_ch_A->ptr_box_idx[aux_int_X];
              }
            }
          }
          ptr_ch_B->cell_size += 8;
          no_cells_ch_A -= 8;
          cell_idx -= 8;
        }
        else if (bdry_cond_type == 0 && ptr_ch_A->pbc_correction_due_pbc_flag == true) //
        {
          box_real_dim_X_ch_A = ptr_ch_A->box_real_dim_x;
          box_real_dim_X_times_Y_ch_A = ptr_ch_A->box_real_dim_x * ptr_ch_A->box_real_dim_y;

          aux_idx_x = ptr_ch_A->ptr_cell_idx_x[cell_idx];
          aux_idx_y = ptr_ch_A->ptr_cell_idx_y[cell_idx];
          aux_idx_z = ptr_ch_A->ptr_cell_idx_z[cell_idx];

          box_idx_x_ch_A = aux_idx_x - ptr_ch_A->box_ts_x;
          box_idx_y_ch_A = aux_idx_y - ptr_ch_A->box_ts_y;
          box_idx_z_ch_A = aux_idx_z - ptr_ch_A->box_ts_z;

          if (ptr_ch_A->pbc_crosses_sim_box_bdry == true)
          {
            if (aux_idx_x > ptr_ch_A->box_max_x)
            {
              box_idx_x_ch_A -= (1 << (lv + 1));
            }

            if (aux_idx_y > ptr_ch_A->box_max_y)
            {
              box_idx_y_ch_A -= (1 << (lv + 1));
            }

            if (aux_idx_z > ptr_ch_A->box_max_z)
            {
              box_idx_z_ch_A -= (1 << (lv + 1));
            }
          }

          // This index corresponds to the corrected one
          box_idx_ch_A = box_idx_x_ch_A + box_idx_y_ch_A * box_real_dim_X_ch_A + box_idx_z_ch_A * box_real_dim_X_times_Y_ch_A;

          //* changing the box index and transfer cells
          if (ptr_ch_A->ptr_box[box_idx_ch_A] == -4)
          {
            // printf("correcting\n");
            for (int kk = 0; kk < 2; kk++)
            {
              for (int jj = 0; jj < 2; jj++)
              {
                for (int ii = 0; ii < 2; ii++)
                {
                  aux_int = cell_idx + ii + jj * 2 + kk * 4;
                  box_idxNbr_ch_A = box_idx_ch_A + ii + jj * box_real_dim_X_ch_A + kk * box_real_dim_X_times_Y_ch_A;
                  box_idxNbr_ch_A_old = ptr_ch_A->ptr_box_idx[aux_int];

                  ptr_ch_A->ptr_box[box_idxNbr_ch_A_old] = -4; // Putting the status of No-EXIST (-4) in the old cell
                  ptr_ch_A->ptr_box[box_idxNbr_ch_A] = -3;     // Putting the status of EXIST (-3) in the corrected cell

                  ptr_ch_A->ptr_box_idx[aux_int] = box_idxNbr_ch_A; // Changing the new box index of the child node A

                  cell_struct_aux = ptr_ch_A->ptr_cell_struct[box_idxNbr_ch_A];
                  ptr_ch_A->ptr_cell_struct[box_idxNbr_ch_A] = ptr_ch_A->ptr_cell_struct[box_idxNbr_ch_A_old];
                  ptr_ch_A->ptr_cell_struct[box_idxNbr_ch_A_old] = cell_struct_aux;
                }
              }
            }
          }
        }
      }
    }
    else // ptr_ch_A->box_check_fit == false
    {
      for (int cell_idx = 0; cell_idx < no_cells_ch_A; cell_idx += 8)
      {

        aux_idx_x = ptr_ch_A->ptr_cell_idx_x[cell_idx] >> 1;
        aux_idx_y = ptr_ch_A->ptr_cell_idx_y[cell_idx] >> 1;
        aux_idx_z = ptr_ch_A->ptr_cell_idx_z[cell_idx] >> 1;

        box_idx_x_node = aux_idx_x - ptr_node->box_ts_x;
        box_idx_y_node = aux_idx_y - ptr_node->box_ts_y;
        box_idx_z_node = aux_idx_z - ptr_node->box_ts_z;

        if (ptr_node->pbc_crosses_sim_box_bdry == true)
        {
          if (aux_idx_x > ptr_node->box_max_x)
          {
            box_idx_x_node -= (1 << lv);
          }

          if (aux_idx_y > ptr_node->box_max_y)
          {
            box_idx_y_node -= (1 << lv);
          }

          if (aux_idx_z > ptr_node->box_max_z)
          {
            box_idx_z_node -= (1 << lv);
          }
        }

        box_idx_node = box_idx_x_node + box_idx_y_node * box_real_dim_X_node + box_idx_z_node * box_real_dim_X_times_Y_node;

        if (ptr_node->ptr_box[box_idx_node] == new_zone_idx)
        {

          box_real_dim_X_ch_A = ptr_ch_A->box_real_dim_x;
          box_real_dim_X_times_Y_ch_A = ptr_ch_A->box_real_dim_x * ptr_ch_A->box_real_dim_y;

          aux_idx_x = ptr_ch_A->ptr_cell_idx_x[cell_idx];
          aux_idx_y = ptr_ch_A->ptr_cell_idx_y[cell_idx];
          aux_idx_z = ptr_ch_A->ptr_cell_idx_z[cell_idx];

          box_idx_x_ch_A = aux_idx_x - ptr_ch_A->box_ts_x;
          box_idx_y_ch_A = aux_idx_y - ptr_ch_A->box_ts_y;
          box_idx_z_ch_A = aux_idx_z - ptr_ch_A->box_ts_z;

          if (ptr_ch_A->pbc_crosses_sim_box_bdry == true)
          {
            if (aux_idx_x > ptr_ch_A->box_max_x)
            {
              box_idx_x_ch_A -= (1 << (lv + 1));
            }

            if (aux_idx_y > ptr_ch_A->box_max_y)
            {
              box_idx_y_ch_A -= (1 << (lv + 1));
            }

            if (aux_idx_z > ptr_ch_A->box_max_z)
            {
              box_idx_z_ch_A -= (1 << (lv + 1));
            }
          }

          box_idx_ch_A = box_idx_x_ch_A + box_idx_y_ch_A * box_real_dim_X_ch_A + box_idx_z_ch_A * box_real_dim_X_times_Y_ch_A;

          //* changing the box and transfer box_mass_aux to box_mass of the cells
          if (ptr_ch_A->ptr_box[box_idx_ch_A] == -4)
          {
            for (int kk = 0; kk < 2; kk++)
            {
              for (int jj = 0; jj < 2; jj++)
              {
                for (int ii = 0; ii < 2; ii++)
                {
                  aux_int = cell_idx + ii + jj * 2 + kk * 4;
                  box_idxNbr_ch_A = box_idx_ch_A + ii + jj * box_real_dim_X_ch_A + kk * box_real_dim_X_times_Y_ch_A;
                  box_idxNbr_ch_A_old = ptr_ch_A->ptr_box_idx[aux_int];
                  ptr_ch_A->ptr_box[box_idxNbr_ch_A] = -3;          // Putting the status of EXIST (-3) in the child node cell
                  ptr_ch_A->ptr_box_idx[aux_int] = box_idxNbr_ch_A; // Changing the new box index of the child node A

                  cell_struct_aux = ptr_ch_A->ptr_cell_struct[box_idxNbr_ch_A];
                  ptr_ch_A->ptr_cell_struct[box_idxNbr_ch_A] = ptr_ch_A->ptr_cell_struct_old[box_idxNbr_ch_A_old];
                  ptr_ch_A->ptr_cell_struct_old[box_idxNbr_ch_A_old] = cell_struct_aux;
                }
              }
            }
          }
        }
        else
        {
          ptr_ch_B = ptr_node->pptr_chn[ptr_node->ptr_links_old_ord_new[ptr_node->ptr_box[box_idx_node]]];

          //* >> CELLS *//
          no_cells_ch_B = ptr_ch_B->cell_size;

          box_real_dim_X_ch_B = ptr_ch_B->box_real_dim_x;
          box_real_dim_X_times_Y_ch_B = ptr_ch_B->box_real_dim_x * ptr_ch_B->box_real_dim_y;

          aux_idx_x = ptr_ch_A->ptr_cell_idx_x[cell_idx];
          aux_idx_y = ptr_ch_A->ptr_cell_idx_y[cell_idx];
          aux_idx_z = ptr_ch_A->ptr_cell_idx_z[cell_idx];

          box_idx_x_ch_B = aux_idx_x - ptr_ch_B->box_ts_x;
          box_idx_y_ch_B = aux_idx_y - ptr_ch_B->box_ts_y;
          box_idx_z_ch_B = aux_idx_z - ptr_ch_B->box_ts_z;

          if (ptr_ch_B->pbc_crosses_sim_box_bdry == true)
          {
            if (aux_idx_x > ptr_ch_B->box_max_x)
            {
              box_idx_x_ch_B -= (1 << (lv + 1));
            }

            if (aux_idx_y > ptr_ch_B->box_max_y)
            {
              box_idx_y_ch_B -= (1 << (lv + 1));
            }

            if (aux_idx_z > ptr_ch_B->box_max_z)
            {
              box_idx_z_ch_B -= (1 << (lv + 1));
            }
          }

          box_idx_ch_B = box_idx_x_ch_B + box_idx_y_ch_B * box_real_dim_X_ch_B + box_idx_z_ch_B * box_real_dim_X_times_Y_ch_B;

          // box_idx_x_ch_B = ptr_ch_A->ptr_cell_idx_x[cell_idx] - ptr_ch_B->box_ts_x;
          // box_idx_y_ch_B = ptr_ch_A->ptr_cell_idx_y[cell_idx] - ptr_ch_B->box_ts_y;
          // box_idx_z_ch_B = ptr_ch_A->ptr_cell_idx_z[cell_idx] - ptr_ch_B->box_ts_z;
          // box_idx_ch_B = box_idx_x_ch_B + box_idx_y_ch_B * box_real_dim_X_ch_B + box_idx_z_ch_B * box_real_dim_X_times_Y_ch_B;

          for (int kk = 0; kk < 2; kk++)
          {
            for (int jj = 0; jj < 2; jj++)
            {
              for (int ii = 0; ii < 2; ii++)
              {
                aux_int = ii + jj * 2 + kk * 4;
                aux_int_A = cell_idx + aux_int;
                aux_int_B = no_cells_ch_B + aux_int;
                aux_int_X = no_cells_ch_A - 8 + aux_int;

                // Transferring from child A to child B
                ptr_ch_B->ptr_cell_idx_x[aux_int_B] = ptr_ch_A->ptr_cell_idx_x[aux_int_A];
                ptr_ch_B->ptr_cell_idx_y[aux_int_B] = ptr_ch_A->ptr_cell_idx_y[aux_int_A];
                ptr_ch_B->ptr_cell_idx_z[aux_int_B] = ptr_ch_A->ptr_cell_idx_z[aux_int_A];

                box_idxNbr_ch_A_old = ptr_ch_A->ptr_box_idx[aux_int_A];
                box_idxNbr_ch_B = box_idx_ch_B + ii + jj * box_real_dim_X_ch_B + kk * box_real_dim_X_times_Y_ch_B;
                ptr_ch_B->ptr_box_idx[aux_int_B] = box_idxNbr_ch_B;

                // Updating the box status of child node B
                ptr_ch_B->ptr_box[box_idxNbr_ch_B] = -3; // Putting the status of EXIST (-3) in the child node cell

                //* >> Cell structure and local mass *//
                cell_struct_aux = ptr_ch_B->ptr_cell_struct[box_idxNbr_ch_B];
                ptr_ch_B->ptr_cell_struct[box_idxNbr_ch_B] = ptr_ch_A->ptr_cell_struct_old[box_idxNbr_ch_A_old];
                ptr_ch_A->ptr_cell_struct_old[box_idxNbr_ch_A_old] = cell_struct_aux;

                // Local mass
                ptr_ch_B->node_mass += ptr_ch_B->ptr_cell_struct[box_idxNbr_ch_B].cell_mass;
                ptr_ch_B->no_ptcl_full_node += ptr_ch_B->ptr_cell_struct[box_idxNbr_ch_B].ptcl_size;

                ptr_ch_A->node_mass -= ptr_ch_B->ptr_cell_struct[box_idxNbr_ch_B].cell_mass;
                ptr_ch_A->no_ptcl_full_node -= ptr_ch_B->ptr_cell_struct[box_idxNbr_ch_B].ptcl_size;

                // Transferring from end of child A to removed cells of child A
                ptr_ch_A->ptr_cell_idx_x[aux_int_A] = ptr_ch_A->ptr_cell_idx_x[aux_int_X];
                ptr_ch_A->ptr_cell_idx_y[aux_int_A] = ptr_ch_A->ptr_cell_idx_y[aux_int_X];
                ptr_ch_A->ptr_cell_idx_z[aux_int_A] = ptr_ch_A->ptr_cell_idx_z[aux_int_X];
                ptr_ch_A->ptr_box_idx[aux_int_A] = ptr_ch_A->ptr_box_idx[aux_int_X];
              }
            }
          }
          ptr_ch_B->cell_size += 8;
          no_cells_ch_A -= 8;
          cell_idx -= 8;
        }
      }
    }
    ptr_ch_A->cell_size = no_cells_ch_A;
    zone_idx++;
  }

  // for (int zone_idx = 0; zone_idx < ptr_node->zones_size; zone_idx++)
  // {
  //     if (ptr_node->ptr_links_old_ord_old[zone_idx] < ptr_node->chn_size)
  //     {

  //     }
  // }

  // Case of old child nodes that will not be used
  cntr_ch = 0;
  for (int ch = 0; ch < ptr_node->chn_size - ptr_node->zones_size; ch++)
  {
    while (ptr_node->ptr_links_old_ord_old[cntr_ch] == (cntr_ch + ch) && cntr_ch < ptr_node->zones_size)
    {
      cntr_ch++;
    }

    ptr_ch_A = ptr_node->pptr_chn[cntr_ch + ch];

    no_cells_ch_A = ptr_ch_A->cell_size;
    for (int cell_idx = 0; cell_idx < no_cells_ch_A; cell_idx += 8)
    {

      aux_idx_x = ptr_ch_A->ptr_cell_idx_x[cell_idx] >> 1;
      aux_idx_y = ptr_ch_A->ptr_cell_idx_y[cell_idx] >> 1;
      aux_idx_z = ptr_ch_A->ptr_cell_idx_z[cell_idx] >> 1;

      box_idx_x_node = aux_idx_x - ptr_node->box_ts_x;
      box_idx_y_node = aux_idx_y - ptr_node->box_ts_y;
      box_idx_z_node = aux_idx_z - ptr_node->box_ts_z;

      if (ptr_node->pbc_crosses_sim_box_bdry == true)
      {
        if (aux_idx_x > ptr_node->box_max_x)
        {
          box_idx_x_node -= (1 << lv);
        }

        if (aux_idx_y > ptr_node->box_max_y)
        {
          box_idx_y_node -= (1 << lv);
        }

        if (aux_idx_z > ptr_node->box_max_z)
        {
          box_idx_z_node -= (1 << lv);
        }
      }

      box_idx_node = box_idx_x_node + box_idx_y_node * box_real_dim_X_node + box_idx_z_node * box_real_dim_X_times_Y_node;

      // box_idx_x_node = (ptr_ch_A->ptr_cell_idx_x[cell_idx] >> 1) - ptr_node->box_ts_x;
      // box_idx_y_node = (ptr_ch_A->ptr_cell_idx_y[cell_idx] >> 1) - ptr_node->box_ts_y;
      // box_idx_z_node = (ptr_ch_A->ptr_cell_idx_z[cell_idx] >> 1) - ptr_node->box_ts_z;
      // box_idx_node = box_idx_x_node + box_idx_y_node * box_real_dim_X_node + box_idx_z_node * box_real_dim_X_times_Y_node;
      if (ptr_node->ptr_box[box_idx_node] >= 0)
      {
        ptr_ch_B = ptr_node->pptr_chn[ptr_node->ptr_links_old_ord_new[ptr_node->ptr_box[box_idx_node]]];

        //* >> CELLS *//
        no_cells_ch_B = ptr_ch_B->cell_size;

        box_real_dim_X_ch_B = ptr_ch_B->box_real_dim_x;
        box_real_dim_X_times_Y_ch_B = ptr_ch_B->box_real_dim_x * ptr_ch_B->box_real_dim_y;

        aux_idx_x = ptr_ch_A->ptr_cell_idx_x[cell_idx];
        aux_idx_y = ptr_ch_A->ptr_cell_idx_y[cell_idx];
        aux_idx_z = ptr_ch_A->ptr_cell_idx_z[cell_idx];

        box_idx_x_ch_B = aux_idx_x - ptr_ch_B->box_ts_x;
        box_idx_y_ch_B = aux_idx_y - ptr_ch_B->box_ts_y;
        box_idx_z_ch_B = aux_idx_z - ptr_ch_B->box_ts_z;

        if (ptr_ch_B->pbc_crosses_sim_box_bdry == true)
        {
          if (aux_idx_x > ptr_ch_B->box_max_x)
          {
            box_idx_x_ch_B -= (1 << (lv + 1));
          }

          if (aux_idx_y > ptr_ch_B->box_max_y)
          {
            box_idx_y_ch_B -= (1 << (lv + 1));
          }

          if (aux_idx_z > ptr_ch_B->box_max_z)
          {
            box_idx_z_ch_B -= (1 << (lv + 1));
          }
        }

        box_idx_ch_B = box_idx_x_ch_B + box_idx_y_ch_B * box_real_dim_X_ch_B + box_idx_z_ch_B * box_real_dim_X_times_Y_ch_B;

        // box_idx_x_ch_B = ptr_ch_A->ptr_cell_idx_x[cell_idx] - ptr_ch_B->box_ts_x;
        // box_idx_y_ch_B = ptr_ch_A->ptr_cell_idx_y[cell_idx] - ptr_ch_B->box_ts_y;
        // box_idx_z_ch_B = ptr_ch_A->ptr_cell_idx_z[cell_idx] - ptr_ch_B->box_ts_z;
        // box_idx_ch_B = box_idx_x_ch_B + box_idx_y_ch_B * box_real_dim_X_ch_B + box_idx_z_ch_B * box_real_dim_X_times_Y_ch_B;

        for (int kk = 0; kk < 2; kk++)
        {
          for (int jj = 0; jj < 2; jj++)
          {
            for (int ii = 0; ii < 2; ii++)
            {
              // Transferring from child A to child B
              aux_int = ii + jj * 2 + kk * 4;
              aux_int_A = cell_idx + aux_int;
              aux_int_B = no_cells_ch_B + aux_int;

              ptr_ch_B->ptr_cell_idx_x[aux_int_B] = ptr_ch_A->ptr_cell_idx_x[aux_int_A];
              ptr_ch_B->ptr_cell_idx_y[aux_int_B] = ptr_ch_A->ptr_cell_idx_y[aux_int_A];
              ptr_ch_B->ptr_cell_idx_z[aux_int_B] = ptr_ch_A->ptr_cell_idx_z[aux_int_A];

              box_idxNbr_ch_A = ptr_ch_A->ptr_box_idx[aux_int_A];
              box_idxNbr_ch_B = box_idx_ch_B + ii + jj * box_real_dim_X_ch_B + kk * box_real_dim_X_times_Y_ch_B;
              ptr_ch_B->ptr_box_idx[aux_int_B] = box_idxNbr_ch_B;

              // Updating the box status of child node B

              ptr_ch_B->ptr_box[box_idxNbr_ch_B] = -3; // Putting the status of EXIST (-3) in the child node cell

              //* >> Cell structure and local mass *//
              cell_struct_aux = ptr_ch_B->ptr_cell_struct[box_idxNbr_ch_B];
              ptr_ch_B->ptr_cell_struct[box_idxNbr_ch_B] = ptr_ch_A->ptr_cell_struct[box_idxNbr_ch_A];
              ptr_ch_A->ptr_cell_struct[box_idxNbr_ch_A] = cell_struct_aux;

              ptr_ch_B->node_mass += ptr_ch_B->ptr_cell_struct[box_idxNbr_ch_B].cell_mass;
              ptr_ch_B->no_ptcl_full_node += ptr_ch_B->ptr_cell_struct[box_idxNbr_ch_B].ptcl_size;
            }
          }
        }
        ptr_ch_B->cell_size += 8;
      }
    }
    ptr_ch_A->cell_size = no_cells_ch_A;
  }
  return _SUCCESS_;
}

/**
 * @brief Transferring new refined parent cells from the new refinement zone to 
 * the new child nodes
 *
 * **SHORT DESCRIPTION**: Transferring new parent cells that require be refined 
 * from the parent node to the new child nodes to which they will belong.
 *
 * **PREREQUISITES**: If there is at least a new refinement zone.
 *
 * @param[in,out] ptr_node Pointer to node structure
 *
 * **RETURN**: The error status.
 *
 * **LONG DESCRIPTION**:
 *
 * Transferring new parent cells that require be refined from the parent node to 
 * the new child nodes to which they will belong. Furthermore, the new box 
 * status is also updated with the new information of child nodes, thus also 
 * like the mass and number of particles of the child nodes and other child node 
 * parameters.
 *
 * Unlike to the function transferring_old_child_to_new_child(), here there is
 * no need to consider the box flag parameter \ref node.box_check_fit
 * "box_check_fit" because there is no sending child node, instead it is the
 * parent node that is in charge of sending the new information to the child 
 * nodes.
 *
 * The flux of this function can be seen in the figure (work in progress), and
 * it is explained below:
 *
 * - [0]  <b> THE transferring_new_zones_to_new_child() FUNCTION STARTS....</b>
 *
 * - [1]  Defining some internal useful parameters.
 *
 * - [2]  Run a "for" loop over every new refined zone.
 *
 * - [3]  Run a "while" loop over current the number of cells between the child 
 * node and the size of the associated new refinement zone in the parent.
 *
 * - [4]  The child box index associated with the parent cell of the "while"
 * loop above is computed.
 *
 * - [5]  For the child node, the box status is updated, and the box indices are
 * computed for the package of 8 cells associated with the parent node cell.
 *
 * - [6]  Finally, the particles are transferred from the parent node cell to 
 * the package of 8 cells of the child node. Moreover, other child node
 * parameters, like the cell masses, are also updated.
 *
 * - [7]  <b> THE transferring_new_zones_to_new_child() FUNCTION ENDS....</b>
 *
 * **ILLUSTRATIVE EXAMPLES**:
 * - [a]  Trivial.
 *
 * **RATIONALES**:
 * - [a]  aaa
 *
 * **NOTES**:
 * - [a]
 */

static int transferring_new_zones_to_new_child(struct node *ptr_node)
{
  struct node *ptr_ch; // child node

  int box_idx_node; // Box index of the node cell

  // int box_idx_x_ch; // Box index in X direcction of the node cell
  // int box_idx_y_ch; // Box index in Y direcction of the node cell
  // int box_idx_z_ch; // Box index in Z direcction of the node cell

  int box_idxNbr_ch; // Box index of the node cell

  int box_idx_ptcl_ch; // Box index of particles in the child node

  int aux_idx_x;
  int aux_idx_y;
  int aux_idx_z;

  int lv = ptr_node->lv;

  int box_idx_x_ch;
  int box_idx_y_ch;
  int box_idx_z_ch;

  int box_idx_ch; // Box index of the node cell

  int cell_idx_node;
  int cell_idx_ch;

  int ptcl_idx;

  int no_cells_ch; // Number of cell in the child node

  int zone_element;

  int box_real_dim_X_ch;
  int box_real_dim_X_times_Y_ch;

  for (int zone_idx = 0; zone_idx < ptr_node->zones_size; zone_idx++)
  {
    ptr_ch = ptr_node->pptr_chn[ptr_node->ptr_links_old_ord_new[zone_idx]];

    no_cells_ch = ptr_ch->cell_size;
    zone_element = 0;
    while (no_cells_ch < ptr_node->ptr_zone_size[zone_idx] * 8)
    {
      cell_idx_node = ptr_node->pptr_zones[zone_idx][zone_element];
      box_idx_node = ptr_node->ptr_box_idx[cell_idx_node];

      // Case cell is a new refinement cell. If it was not new, the past child nodes would have already added it to the corresponding child node
      if (ptr_node->ptr_box_old[box_idx_node] < 0)
      {
        box_real_dim_X_ch = ptr_ch->box_real_dim_x;
        box_real_dim_X_times_Y_ch = ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;

        //* >> CELLS, BOX STATUS AND CAPACITY OF CELL STRUCTURE *//
        aux_idx_x = ptr_node->ptr_cell_idx_x[cell_idx_node] * 2;
        aux_idx_y = ptr_node->ptr_cell_idx_y[cell_idx_node] * 2;
        aux_idx_z = ptr_node->ptr_cell_idx_z[cell_idx_node] * 2;

        box_idx_x_ch = aux_idx_x - ptr_ch->box_ts_x;
        box_idx_y_ch = aux_idx_y - ptr_ch->box_ts_y;
        box_idx_z_ch = aux_idx_z - ptr_ch->box_ts_z;

        if (ptr_ch->pbc_crosses_sim_box_bdry == true)
        {
          if (aux_idx_x > ptr_ch->box_max_x)
          {
            box_idx_x_ch -= (1 << (lv + 1));
          }

          if (aux_idx_y > ptr_ch->box_max_y)
          {
            box_idx_y_ch -= (1 << (lv + 1));
          }

          if (aux_idx_z > ptr_ch->box_max_z)
          {
            box_idx_z_ch -= (1 << (lv + 1));
          }
        }

        box_idx_ch = box_idx_x_ch + box_idx_y_ch * box_real_dim_X_ch + box_idx_z_ch * box_real_dim_X_times_Y_ch;

        for (int kk = 0; kk < 2; kk++)
        {
          for (int jj = 0; jj < 2; jj++)
          {
            for (int ii = 0; ii < 2; ii++)
            {
              cell_idx_ch = no_cells_ch + ii + jj * 2 + kk * 4;
              // Cells
              ptr_ch->ptr_cell_idx_x[cell_idx_ch] = aux_idx_x + ii;
              ptr_ch->ptr_cell_idx_y[cell_idx_ch] = aux_idx_y + jj;
              ptr_ch->ptr_cell_idx_z[cell_idx_ch] = aux_idx_z + kk;

              // Box status
              // box_idx_x_ch = ptr_ch->ptr_cell_idx_x[cell_idx_ch] - ptr_ch->box_ts_x;
              // box_idx_y_ch = ptr_ch->ptr_cell_idx_y[cell_idx_ch] - ptr_ch->box_ts_y;
              // box_idx_z_ch = ptr_ch->ptr_cell_idx_z[cell_idx_ch] - ptr_ch->box_ts_z;
              box_idxNbr_ch = box_idx_ch + ii + jj * box_real_dim_X_ch + kk * box_real_dim_X_times_Y_ch;

              ptr_ch->ptr_box_idx[cell_idx_ch] = box_idxNbr_ch;

              ptr_ch->ptr_box[box_idxNbr_ch] = -3;

              // Cell structure capacity
              //* >> Space checking of the capacity of the particles in the child cells *//
              if (space_check(&(ptr_ch->ptr_cell_struct[box_idxNbr_ch].ptcl_cap), ptr_node->ptr_cell_struct[box_idx_node].ptcl_size + 1, 1.0f, "p1i1", &(ptr_ch->ptr_cell_struct[box_idxNbr_ch].ptr_ptcl)) == _FAILURE_)
              {
                printf("Error, in space_check function\n");
                return _FAILURE_;
              }
            }
          }
        }
        no_cells_ch += 8;

        // Transfer particles from node to child
        for (int j = 0; j < ptr_node->ptr_cell_struct[box_idx_node].ptcl_size; j++)
        {
          ptcl_idx = ptr_node->ptr_cell_struct[box_idx_node].ptr_ptcl[j];
          box_idx_ptcl_ch = ptcl_idx_to_box_idx(ptr_ch, ptcl_idx);
          ptr_ch->ptr_cell_struct[box_idx_ptcl_ch].cell_mass += GL_ptcl_mass[ptcl_idx];
          // ptr_ch->node_mass += GL_ptcl_mass[ptcl_idx];
          // ptr_ch->no_ptcl_full_node += 1;
          ptr_ch->ptr_cell_struct[box_idx_ptcl_ch].ptr_ptcl[ptr_ch->ptr_cell_struct[box_idx_ptcl_ch].ptcl_size] = ptcl_idx;
          ptr_ch->ptr_cell_struct[box_idx_ptcl_ch].ptcl_size += 1;
        }
        ptr_ch->node_mass += ptr_node->ptr_cell_struct[box_idx_node].cell_mass;
        ptr_ch->no_ptcl_full_node += ptr_node->ptr_cell_struct[box_idx_node].ptcl_size;
      }
      zone_element++;
    }
    ptr_ch->cell_size = no_cells_ch;
  }
  return _SUCCESS_;
}

/**
 * @brief Properly sort child nodes into the parent node array \ref
 * node.pptr_chn "pptr_chn"
 *
 * **SHORT DESCRIPTION**: As the new child nodes can have different IDs of
 * identification, this function properly sorts them into the parent node array
 * \ref node.pptr_chn "pptr_chn"
 *
 * **PREREQUISITES**: If there is at least a new refinement zone.
 *
 * @param[in,out] ptr_node Pointer to node structure
 *
 * **RETURN**: No parameter is returned.
 *
 * **LONG DESCRIPTION**:
 *
 * At this point, the identification of all child nodes corresponds to the same
 * new refinement zone identification of the parent node. However, how the 
 * parent node has access to these child nodes can be different from this sort
 * This function is in charge of properly sorting the child node into the parent 
 * node array \ref node.pptr_chn "pptr_chn".
 *
 * The flux of this function can be seen in the figure (work in progress), and
 * it is explained below:
 *
 * - [0]  <b> THE reorganization_child_node() FUNCTION STARTS....</b>
 *
 * - [1]  Defining some internal useful parameters.
 *
 * - [2]  An auxiliary struct node pointer to the pointer array is allocated and
 * filled with the correct sort of child nodes. First, it is filled with the
 * linked child nodes, and then, if there are, it is filled
 * with the old child nodes which will be removed from the parent node.
 *
 * - [3]  Finally, the \ref node.pptr_chn "pptr_chn" is sorted using the 
 * already properly sorted auxiliary array.
 *
 * - [4]  <b> THE reorganization_child_node() FUNCTION ENDS....</b>
 *
 * **ILLUSTRATIVE EXAMPLES**:
 * - [a]  Trivial.
 *
 * **RATIONALES**:
 * - [a]  aaa
 *
 * **NOTES**:
 * - [a]  This function can be improved in the future by avoiding the use of the 
 * auxiliary pointer to the pointer array, and simply using a pointer an array with no
 * additional memory allocation or release of it. 
 */

static void reorganization_child_node(struct node *ptr_node)
{
  // struct node *ptr_ch;
  struct node **pptr_aux;

  int cntr_ch;

  pptr_aux = (struct node **)malloc((ptr_node->ptr_links_old_ord_old[ptr_node->zones_size - 1] + 1) * sizeof(struct node *));

  // // Updating the ID of the child nodes
  // for (int zone_idx = 0; zone_idx < ptr_node->zones_size; zone_idx++)
  // {
  //     ptr_ch = ptr_node->pptr_chn[ptr_node->ptr_links_old_ord_old[zone_idx]];
  //     ptr_ch->ID = ptr_node->ptr_links_new_ord_old[zone_idx];
  // }

  // Filling the pptr_aux array with child nodes linked with the new refinement zones
  for (int i = 0; i < ptr_node->zones_size; i++)
  {
    pptr_aux[i] = ptr_node->pptr_chn[ptr_node->ptr_links_old_ord_new[i]];
  }

  // Filling the pptr_aux array with child nodes without link
  cntr_ch = 0;
  for (int i = 0; i < ptr_node->ptr_links_old_ord_old[ptr_node->zones_size - 1] + 1 - ptr_node->zones_size; i++)
  {
    while (ptr_node->ptr_links_old_ord_old[cntr_ch] == cntr_ch + i)
    {
      cntr_ch++;
    }
    pptr_aux[ptr_node->zones_size + i] = ptr_node->pptr_chn[cntr_ch + i];
  }

  // Reorganization of child nodes
  for (int i = 0; i < ptr_node->ptr_links_old_ord_old[ptr_node->zones_size - 1] + 1; i++)
  {
    ptr_node->pptr_chn[i] = pptr_aux[i];
  }

  free(pptr_aux);
}

/**
 * @brief Properly sort every child node with the grandchild nodes
 *
 * **SHORT DESCRIPTION**: Equal to the function reorganization_child_node(), the
 * grandchild nodes can also need to be sorted properly relative to the child
 * nodes (their parent nodes).
 *
 * **PREREQUISITES**: If there is at least one new and one old refinement zone,
 * and if the level of refinement is lower than the maximum level of refinement
 * minus 1.
 *
 * @param[in,out] ptr_node Pointer to node structure
 *
 * **RETURN**: The error status.
 *
 * **LONG DESCRIPTION**:
 *
 * Equal than at the function reorganization_child_node(), the grandchild nodes
 * can also need to be sorted properly to be identified with the child nodes
 * (their parent nodes).
 *
 * The grandchild nodes do not change in their geometry, node properties, or
 * amount. They are preserved completely except their only modification can be
 * that they change their parents. How to proceed is very similar to the
 * previous function reorganization_child_node(), using an auxiliary pointer to
 * the struct node pointer array to store the correct sort of grandchild nodes.
 *
 * Here, it is important to see that since there is no direct parameter to know
 * how many grandchildren the parent node has, it is necessary to count them.
 *
 * The flux of this function can be seen in the figure (work in progress), and
 * it is explained below:
 *
 * - [0]  <b> THE reorganization_grandchild_node() FUNCTION STARTS....</b>
 *
 * - [1]  Defining some internal useful parameters.
 *
 * - [2]  Counting the total number of grandchild nodes of the parent node using
 * both child nodes, the linked and those who will be removed.
 *
 * - [3]  An auxiliary struct node pointer to the pointer array is allocated and
 * filled with all grandchild nodes in any order. Moreover, the number of child
 * nodes of all child nodes is reset to 0.
 *
 * - [4]  Run a "for" loop over the number of grandchild nodes.
 *
 * - [5]  The parent box index associated with the first cell of the grandchild
 * node is computed
 *
 * - [6]  Using this index in the parent box, the corresponding child node
 * identification for that grandchild node is localized.
 *
 * - [7]  The grandchild node is linked to this child node and the
 * identification parameter of the grandchild node is updated to correspond to
 * the child node.
 *
 * - [8]  Finally, the number of child nodes of the child node is increasing in
 * 1 and the "for" loop of step [4] begins again with another grandchild node.
 *
 * - [9]  <b> THE reorganization_grandchild_node() FUNCTION ENDS....</b>
 *
 * **ILLUSTRATIVE EXAMPLES**:
 * - [a]  Imagine the situation where in the old time step there is only one
 * child node A, with two children GC1 and GC2, which are the grandchild of the
 * parent node. Then, in the next time step, the new map of refinement of the
 * parent node breaks the child node A into two parts, now known as A1 and A2.
 * In the function adapt_child_nodes(), the old child node A has been adapted to
 * fit, for example, with the new zone A1. In the function
 * create_new_child_nodes(), a new child node is created to fit the other new
 * zone A2. At this point, it is not known if the grandchild nodes GC1 and GC2
 * belong to A1 or A2 or a combination of them.
 *
 * **RATIONALES**:
 * - [a]  aaa
 *
 * **NOTES**:
 * - [a]  This function can be improved in the future by avoiding the use of the
 * auxiliary pointer to the pointer array, and simply using a pointer to an
 * array with no additional memory allocation or release of it.
 *
 * - [b]  Maybe there is more efficient to store the total number of 
 * grandchildren of the parent node instead to compute it here. But, it is also 
 * important to note, that in every time step, the number of the grandchild 
 * nodes of the child nodes can change.
 */

static int reorganization_grandchild_node(struct node *ptr_node)
{

  struct node **pptr_aux;

  int no_grandchildren = 0;

  int box_idx_x_node; // Box index in X direcction of the node cell
  int box_idx_y_node; // Box index in Y direcction of the node cell
  int box_idx_z_node; // Box index in Z direcction of the node cell
  int box_idx_node;   // Box index of the node cell

  int aux_idx_x;
  int aux_idx_y;
  int aux_idx_z;

  int child_ID;

  int size;

  int box_real_dim_X_node = ptr_node->box_real_dim_x;
  int box_real_dim_X_times_Y_node = ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;

  int lv = ptr_node->lv;

  int cntr_grandch = 0; // Counter the number of grandchildren

  for (int ch = 0; ch < ptr_node->zones_size; ch++)
  {
    for (int grandch = 0; grandch < ptr_node->pptr_chn[ch]->chn_size; grandch++)
    {
      no_grandchildren++;
    }
  }
  for (int ch = ptr_node->zones_size; ch < ptr_node->chn_size; ch++)
  {
    for (int grandch = 0; grandch < ptr_node->pptr_chn[ch]->chn_size; grandch++)
    {
      no_grandchildren++;
    }
  }

  if (no_grandchildren > 0)
  {

    pptr_aux = (struct node **)malloc(no_grandchildren * sizeof(struct node *));

    for (int ch = 0; ch < ptr_node->zones_size; ch++)
    {
      for (int grandch = 0; grandch < ptr_node->pptr_chn[ch]->chn_size; grandch++)
      {
        pptr_aux[cntr_grandch] = ptr_node->pptr_chn[ch]->pptr_chn[grandch];
        ptr_node->pptr_chn[ch]->pptr_chn[grandch] = NULL;
        cntr_grandch++;
      }
      ptr_node->pptr_chn[ch]->chn_size = 0; // Reseting the number of grandchildren to 0
    }
    for (int ch = ptr_node->zones_size; ch < ptr_node->chn_size; ch++)
    {
      for (int grandch = 0; grandch < ptr_node->pptr_chn[ch]->chn_size; grandch++)
      {
        pptr_aux[cntr_grandch] = ptr_node->pptr_chn[ch]->pptr_chn[grandch];
        ptr_node->pptr_chn[ch]->pptr_chn[grandch] = NULL;
        cntr_grandch++;
      }
      ptr_node->pptr_chn[ch]->chn_size = 0; // Reseting the number of grandchildren to 0
    }

    for (int grandch = 0; grandch < no_grandchildren; grandch++)
    {

      aux_idx_x = pptr_aux[grandch]->ptr_cell_idx_x[0] >> 2;
      aux_idx_y = pptr_aux[grandch]->ptr_cell_idx_y[0] >> 2;
      aux_idx_z = pptr_aux[grandch]->ptr_cell_idx_z[0] >> 2;

      box_idx_x_node = aux_idx_x - ptr_node->box_ts_x;
      box_idx_y_node = aux_idx_y - ptr_node->box_ts_y;
      box_idx_z_node = aux_idx_z - ptr_node->box_ts_z;

      if (ptr_node->pbc_crosses_sim_box_bdry == true)
      {
        if (aux_idx_x > ptr_node->box_max_x)
        {
          box_idx_x_node -= (1 << lv);
        }

        if (aux_idx_y > ptr_node->box_max_y)
        {
          box_idx_y_node -= (1 << lv);
        }

        if (aux_idx_z > ptr_node->box_max_z)
        {
          box_idx_z_node -= (1 << lv);
        }
      }

      box_idx_node = box_idx_x_node + box_idx_y_node * box_real_dim_X_node + box_idx_z_node * box_real_dim_X_times_Y_node;

      // box_idx_x_node = (pptr_aux[grandch]->ptr_cell_idx_x[0] >> 2) - ptr_node->box_ts_x;
      // box_idx_y_node = (pptr_aux[grandch]->ptr_cell_idx_y[0] >> 2) - ptr_node->box_ts_y;
      // box_idx_z_node = (pptr_aux[grandch]->ptr_cell_idx_z[0] >> 2) - ptr_node->box_ts_z;
      // box_idx_node = box_idx_x_node + box_idx_y_node * box_real_dim_X_node + box_idx_z_node * box_real_dim_X_times_Y_node;
      child_ID = ptr_node->ptr_box[box_idx_node];
      size = ptr_node->pptr_chn[child_ID]->chn_size;

      //* >> Space checking of the child capacity of the child node "child_ID" *//

      if (space_check(&(ptr_node->pptr_chn[child_ID]->chn_cap), size + 1, 2.0f, "p1n2", &(ptr_node->pptr_chn[child_ID]->pptr_chn)) == _FAILURE_)
      {
        printf("Error, in space_check function\n");
        return _FAILURE_;
      }
      ptr_node->pptr_chn[child_ID]->pptr_chn[size] = pptr_aux[grandch];
      pptr_aux[grandch]->ptr_pt = ptr_node->pptr_chn[child_ID];
      pptr_aux[grandch]->ID = size;
      ptr_node->pptr_chn[child_ID]->chn_size = size + 1;
    }
    free(pptr_aux);
  }

  return _SUCCESS_;
}

/**
 * @brief Updating the refined zones identification of the child nodes.
 *
 * **SHORT DESCRIPTION**: Updating the child box status of the refined zones of the
 * child nodes.
 *
 * **PREREQUISITES**: If there is at least one new and one old refinement zone,
 * and if the level of refinement is lower than the maximum level of refinement
 * minus 1.
 *
 * @param[in,out] ptr_node Pointer to node structure
 *
 * **RETURN**: No parameter is returned.
 *
 * **LONG DESCRIPTION**:
 *
 * In the previous function reorganization_grandchild_node(), the child nodes
 * were filled with the proper grandchild nodes, and the identification of these
 * last was updated to match with the sort of the array \ref node.pptr_chn "pptr_chn" of the child
 * node. However, the box status of the child node has not been updated yet with
 * the information on its refined zones. This function is in charge to perform
 * that.
 *
 * The flux of this function can be seen in the figure (work in progress), and
 * it is explained below:
 *
 * - [0]  <b> THE updating_ref_zones_children() FUNCTION STARTS....</b>
 *
 * - [1]  Defining some internal useful parameters.
 *
 * - [2]  Run a "for" loop over the number of refinement zones, that at this
 * point of the module tree_adaptation.c is equal to the final number of child
 * nodes.
 *
 * - [3]  Run a "for" loop over the number of grandchild nodes.
 *
 * - [4]  Run a "for" loop over the cells in every grandchild node. The loop
 * runs each 8 cells which correspond to one cell in the corresponding child node.
 *
 * - [5]  The child box index associated with this package of 8 grandchild cells is
 * computed.
 *
 * - [6]  Using this index the status of the child node box is updated to be
 * equal to the grandchild identification.
 *
 * - [7]  <b> THE updating_ref_zones_children() FUNCTION ENDS....</b>
 *
 * **ILLUSTRATIVE EXAMPLES**:
 * - [a]  a
 *
 * **RATIONALES**:
 * - [a]  aaa
 *
 * **NOTES**:
 * - [a]  a
 */

static void updating_ref_zones_children(struct node *ptr_node)
{
  struct node *ptr_ch;
  struct node *ptr_grandch;

  int box_idx_x_ch; // Box index in X direcction of the child cell
  int box_idx_y_ch; // Box index in Y direcction of the child cell
  int box_idx_z_ch; // Box index in Z direcction of the child cell
  int box_idx_ch;   // Box index of the child cell

  int aux_idx_x;
  int aux_idx_y;
  int aux_idx_z;

  int lv = ptr_node->lv;

  int box_real_dim_X_ch;
  int box_real_dim_X_times_Y_ch;

  for (int ch = 0; ch < ptr_node->zones_size; ch++)
  {
    ptr_ch = ptr_node->pptr_chn[ch];

    box_real_dim_X_ch = ptr_ch->box_real_dim_x;
    box_real_dim_X_times_Y_ch = ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;

    for (int grandch = 0; grandch < ptr_ch->chn_size; grandch++)
    {
      ptr_grandch = ptr_ch->pptr_chn[grandch];
      for (int cell_idx = 0; cell_idx < ptr_grandch->cell_size; cell_idx += 8)
      {

        aux_idx_x = ptr_grandch->ptr_cell_idx_x[cell_idx] >> 1;
        aux_idx_y = ptr_grandch->ptr_cell_idx_y[cell_idx] >> 1;
        aux_idx_z = ptr_grandch->ptr_cell_idx_z[cell_idx] >> 1;

        box_idx_x_ch = aux_idx_x - ptr_ch->box_ts_x;
        box_idx_y_ch = aux_idx_y - ptr_ch->box_ts_y;
        box_idx_z_ch = aux_idx_z - ptr_ch->box_ts_z;

        if (ptr_ch->pbc_crosses_sim_box_bdry == true)
        {
          if (aux_idx_x > ptr_ch->box_max_x)
          {
            box_idx_x_ch -= (1 << (lv + 1));
          }

          if (aux_idx_y > ptr_ch->box_max_y)
          {
            box_idx_y_ch -= (1 << (lv + 1));
          }

          if (aux_idx_z > ptr_ch->box_max_z)
          {
            box_idx_z_ch -= (1 << (lv + 1));
          }
        }

        box_idx_ch = box_idx_x_ch + box_idx_y_ch * box_real_dim_X_ch + box_idx_z_ch * box_real_dim_X_times_Y_ch;

        // box_idx_x_ch = (ptr_grandch->ptr_cell_idx_x[cell_idx] >> 1) - ptr_ch->box_ts_x;
        // box_idx_y_ch = (ptr_grandch->ptr_cell_idx_y[cell_idx] >> 1) - ptr_ch->box_ts_y;
        // box_idx_z_ch = (ptr_grandch->ptr_cell_idx_z[cell_idx] >> 1) - ptr_ch->box_ts_z;
        // box_idx_ch = box_idx_x_ch + box_idx_y_ch * box_real_dim_X_ch + box_idx_z_ch * box_real_dim_X_times_Y_ch;
        ptr_ch->ptr_box[box_idx_ch] = grandch;
      }
    }
  }
}

/**
 * @brief Updating the child box status of the boundary simulation cells
 *
 * **SHORT DESCRIPTION**: Updating the child box status of the boundary 
 * simulation cells
 *
 * **PREREQUISITES**: If there is at least a new refinement zone.
 *
 * @param[in,out] ptr_node Pointer to node structure
 *
 * **RETURN**: No parameter is returned.
 *
 * **LONG DESCRIPTION**:
 *
 * If any of the child nodes touch the boundary simulation in the non-periodic
 * boundary condition (pbc), or crosses the whole simulation in the periodic 
 * boundary condition, then its box requires to be updated in the outer boundary 
 * cells. This function is in charge to perform this.
 *
 * The flux of this function can be seen in the figure (work in progress), and 
 * it is explained below:
 *
 * - [0]  <b> THE adding_boundary_simulation_box_status_to_children_nodes() FUNCTION STARTS....</b>
 *
 * - [1]  Defining some internal useful parameters.
 *
 * - [2]  Run a "for" loop over the number of refinement zones, that at this
 * point of the module tree_adaptation.c is equal to the final number of child
 * nodes.
 *
 * - [3]  If the periodic boundary conditions (pbc) are used, and the child node 
 * flag \ref node.pbc_crosses_whole_sim_box "pbc_crosses_whole_sim_box" is 
 * "true", then the box cell status from the beginning of the box to the 
 * boundary of the *Smallest Box* (see Key Concepts \ref 
 * Key_Concepts_Smallest_Box "Smallest Box"), and between the end of the 
 * "Smallest box" to the end of the box are updated to the value of -6. This is 
 * done in every direction if it corresponds.
 *
 * - [4]  If non-periodic boundary conditions (non-pbc) are used, and the child 
 * node flag \ref node.sim_bdry_contact "sim_bdry_contact" is "true", then the 
 * child box status is updated with the value of -5 when the outer of the 
 * *Smallest Box* (see Key Concepts \ref Key_Concepts_Smallest_Box 
 * "Smallest Box").
 *
 * - [5]  <b> THE adding_boundary_simulation_box_status_to_children_nodes() FUNCTION ENDS....</b>
 *
 * **ILLUSTRATIVE EXAMPLES**:
 * - [a]  a
 *
 * **RATIONALES**:
 * - [a]  aaa
 *
 * **NOTES**:
 * - [a]  a
 */

static void adding_boundary_simulation_box_status_to_children_nodes(struct node *ptr_node)
{

  struct node *ptr_ch;

  int box_real_dim_X_ch;
  int box_real_dim_X_times_Y_ch;

  int box_idx_ch;

  int lv = ptr_node->lv;

  for (int zone_idx = 0; zone_idx < ptr_node->zones_size; zone_idx++)
  {
    ptr_ch = ptr_node->pptr_chn[zone_idx];

    //* >> Adding -6 and -5 to the box boundary when it corresponds
    if (ptr_ch->sim_bdry_contact == true)
    {

      box_real_dim_X_ch = ptr_ch->box_real_dim_x;
      box_real_dim_X_times_Y_ch = ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;

      // Adding status of -6
      if (bdry_cond_type == 0 && ptr_ch->pbc_crosses_whole_sim_box == true)
      {
        // X axis
        if (ptr_ch->pbc_crosses_whole_sim_box_x == true)
        {
          for (int k = 0; k < ptr_ch->box_real_dim_z; k++)
          {
            for (int j = 0; j < ptr_ch->box_real_dim_y; j++)
            {
              // for (int i = 0; i < (ptr_ch->box_real_dim_x - ptr_ch->box_dim_x) / 2; i++)
              for (int i = 0; i < ptr_ch->box_min_x - ptr_ch->box_ts_x; i++)
              {
                box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
                ptr_ch->ptr_box[box_idx_ch] = -6;
              }

              // for (int i = (ptr_ch->box_real_dim_x + ptr_ch->box_dim_x) / 2; i < ptr_ch->box_real_dim_x; i++)
              for (int i = ptr_ch->box_max_x + 1 - ptr_ch->box_ts_x; i < ptr_ch->box_real_dim_x; i++)
              {
                box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
                ptr_ch->ptr_box[box_idx_ch] = -6;
              }
            }
          }
        }

        // Y axis
        if (ptr_ch->pbc_crosses_whole_sim_box_y == true)
        {
          for (int k = 0; k < ptr_ch->box_real_dim_z; k++)
          {
            for (int i = 0; i < ptr_ch->box_real_dim_x; i++)
            {
              // for (int j = 0; j < (ptr_ch->box_real_dim_y - ptr_ch->box_dim_y) / 2; j++)
              for (int j = 0; j < ptr_ch->box_min_y - ptr_ch->box_ts_y; j++)
              {
                box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
                ptr_ch->ptr_box[box_idx_ch] = -6;
              }

              // for (int j = (ptr_ch->box_real_dim_y + ptr_ch->box_dim_y) / 2; j < ptr_ch->box_real_dim_y; j++)
              for (int j = ptr_ch->box_max_y + 1 - ptr_ch->box_ts_y; j < ptr_ch->box_real_dim_y; j++)
              {
                box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
                ptr_ch->ptr_box[box_idx_ch] = -6;
              }
            }
          }
        }

        // Z axis
        if (ptr_ch->pbc_crosses_whole_sim_box_z == true)
        {
          for (int j = 0; j < ptr_ch->box_real_dim_y; j++)
          {
            for (int i = 0; i < ptr_ch->box_real_dim_x; i++)
            {
              // for (int k = 0; k < (ptr_ch->box_real_dim_z - ptr_ch->box_dim_z) / 2; k++)
              for (int k = 0; k < ptr_ch->box_min_z - ptr_ch->box_ts_z; k++)
              {
                box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
                ptr_ch->ptr_box[box_idx_ch] = -6;
              }

              // for (int k = (ptr_ch->box_real_dim_z + ptr_ch->box_dim_z) / 2; k < ptr_ch->box_real_dim_z; k++)
              for (int k = ptr_ch->box_max_z + 1 - ptr_ch->box_ts_z; k < ptr_ch->box_real_dim_z; k++)
              {
                box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
                ptr_ch->ptr_box[box_idx_ch] = -6;
              }
            }
          }
        }
      }
      // Removing status of -5
      else if (bdry_cond_type != 0)
      {
        // X axis
        if (ptr_ch->sim_bdry_contact_x == true)
        {
          if (ptr_ch->box_min_x == 0)
          {
            for (int k = 0; k < ptr_ch->box_real_dim_z; k++)
            {
              for (int j = 0; j < ptr_ch->box_real_dim_y; j++)
              {
                // for (int i = 0; i < (ptr_ch->box_real_dim_x - ptr_ch->box_dim_x) / 2; i++)
                for (int i = 0; i < ptr_ch->box_min_x - ptr_ch->box_ts_x; i++)
                {
                  box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
                  ptr_ch->ptr_box[box_idx_ch] = -5;
                }
              }
            }
          }

          if (ptr_ch->box_max_x == (1 << (lv + 1)) - 1)
          {
            for (int k = 0; k < ptr_ch->box_real_dim_z; k++)
            {
              for (int j = 0; j < ptr_ch->box_real_dim_y; j++)
              {
                // for (int i = (ptr_ch->box_real_dim_x + ptr_ch->box_dim_x) / 2; i < ptr_ch->box_real_dim_x; i++)
                for (int i = ptr_ch->box_max_x + 1 - ptr_ch->box_ts_x; i < ptr_ch->box_real_dim_x; i++)
                {
                  box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
                  ptr_ch->ptr_box[box_idx_ch] = -5;
                }
              }
            }
          }
        }

        // Y axis
        if (ptr_ch->sim_bdry_contact_y == true)
        {
          if (ptr_ch->box_min_y == 0)
          {
            for (int k = 0; k < ptr_ch->box_real_dim_z; k++)
            {
              for (int i = 0; i < ptr_ch->box_real_dim_x; i++)
              {
                // for (int j = 0; j < (ptr_ch->box_real_dim_y - ptr_ch->box_dim_y) / 2; j++)
                for (int j = 0; j < ptr_ch->box_min_y - ptr_ch->box_ts_y; j++)
                {
                  box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
                  ptr_ch->ptr_box[box_idx_ch] = -5;
                }
              }
            }
          }

          if (ptr_ch->box_max_y == (1 << (lv + 1)) - 1)
          {
            for (int k = 0; k < ptr_ch->box_real_dim_z; k++)
            {
              for (int i = 0; i < ptr_ch->box_real_dim_x; i++)
              {
                // for (int j = (ptr_ch->box_real_dim_y + ptr_ch->box_dim_y) / 2; j < ptr_ch->box_real_dim_y; j++)
                for (int j = ptr_ch->box_max_y + 1 - ptr_ch->box_ts_y; j < ptr_ch->box_real_dim_y; j++)
                {
                  box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
                  ptr_ch->ptr_box[box_idx_ch] = -5;
                }
              }
            }
          }
        }

        // Z axis
        if (ptr_ch->sim_bdry_contact_z == true)
        {
          if (ptr_ch->box_min_z == 0)
          {
            for (int j = 0; j < ptr_ch->box_real_dim_y; j++)
            {
              for (int i = 0; i < ptr_ch->box_real_dim_x; i++)
              {
                // for (int k = 0; k < (ptr_ch->box_real_dim_z - ptr_ch->box_dim_z) / 2; k++)
                for (int k = 0; k < ptr_ch->box_min_z - ptr_ch->box_ts_z; k++)
                {
                  box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
                  ptr_ch->ptr_box[box_idx_ch] = -5;
                }
              }
            }
          }

          if (ptr_ch->box_max_z == (1 << (lv + 1)) - 1)
          {
            for (int j = 0; j < ptr_ch->box_real_dim_y; j++)
            {
              for (int i = 0; i < ptr_ch->box_real_dim_x; i++)
              {
                // for (int k = (ptr_ch->box_real_dim_z + ptr_ch->box_dim_z) / 2; k < ptr_ch->box_real_dim_z; k++)
                for (int k = ptr_ch->box_max_z + 1 - ptr_ch->box_ts_z; k < ptr_ch->box_real_dim_z; k++)
                {
                  box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
                  ptr_ch->ptr_box[box_idx_ch] = -5;
                }
              }
            }
          }
        }
      }
    }
  }
}

/**
 * @brief Filling child node grid point arrays
 *
 * **SHORT DESCRIPTION**: Filling child node interior, boundary, and simulation
 * boundary grid point arrays. For example, some of the filled arrays are \ref
 * node.ptr_intr_grid_cell_idx_x "ptr_intr_grid_cell_idx_x", \ref
 * node.ptr_bdry_box_grid_idx "ptr_bdry_box_grid_idx", or \ref
 * node.ptr_sim_bdry_grid_cell_idx_z "ptr_sim_bdry_grid_cell_idx_z".
 *
 * **PREREQUISITES**: If there is at least a new refinement zone.
 *
 * @param[in,out] ptr_node Pointer to node structure
 *
 * **RETURN**: The error status.
 *
 * **LONG DESCRIPTION**:
 *
 * Filling child node interior, boundary, and simulation
 * boundary grid point arrays. This function fills 4 arrays per type of grid 
 * point. For example, for the interior grid points of the child nodes, the 
 * function fills the arrays \ref node.ptr_intr_grid_cell_idx_x 
 * "ptr_intr_grid_cell_idx_x" (\ref node.ptr_intr_grid_cell_idx_y 
 * "ptr_intr_grid_cell_idx_y", \ref node.ptr_intr_grid_cell_idx_z 
 * "ptr_intr_grid_cell_idx_z",), and \ref node.ptr_intr_box_grid_idx 
 * "ptr_intr_box_grid_idx". Moreover, other grid parameters of the child nodes 
 * are updated, the size (see Key Concepts \ref Key_Concepts_Size "Size") and
 * the capacity (see Key Concepts \ref Key_Concepts_Capacity "Capacity") of
 * every type of grid point.
 *
 * The flux of this function can be seen in the figure (work in progress), and
 * it is explained below:
 *
 * - [0]  <b> THE filling_child_grid_point_arrays() FUNCTION STARTS....</b>
 *
 * - [1]  Defining some internal useful parameters.
 *
 * - [2]  Run a "for" loop over the number of refinement zones, that at this
 * point of the module tree_adaptation.c is equal to the final number of child
 * nodes.
 *
 * - [3]  Run a triple "for" loop over the box grid points of the child nodes,
 * but only running over the *Smallest Box* (see Key Concepts \ref 
 * Key_Concepts_Smallest_Box "Smallest Box").
 *
 * - [4]  In every grid point of the child node, depending on the box status of
 * the neighboring cells of that grid point, the grid points are labeled 
 * no-exist, interior, boundary, and simulation boundary grid points using local 
 * parameters.
 *
 * - [5]  If the grid point is labeled no-exist, then the "for" loop continues 
 * with the next grid point. If the grid point is labeled interior, boundary, or
 * simulation boundary, the grid point is added to one of those categories 
 * putting its properties in the corresponding arrays.
 *
 * - [6]  <b> THE filling_child_grid_point_arrays() FUNCTION ENDS....</b>
 *
 * **ILLUSTRATIVE EXAMPLES**:
 * - [a]  a
 *
 * **RATIONALES**:
 * - [a]  aaa
 *
 * **NOTES**:
 * - [a]  It is possible to improve this function by looking for a better way to
 * run between the grid points, and also looking for a better way to find the 
 * type of grid point that it belongs to.
 */

static int filling_child_grid_point_arrays(struct node *ptr_node)
{
  struct node *ptr_ch;

  int box_idxNbr_i0_j0_k0_ch;    // Box index of the child node at x=0,y=0,z=0
  int box_idxNbr_im1_j0_k0_ch;   // Box index of the child node of the neighbor at x=-1,y=0,z=0
  int box_idxNbr_i0_jm1_k0_ch;   // Box index of the child node of the neighbor at x=0,y=-1,z=0
  int box_idxNbr_im1_jm1_k0_ch;  // Box index of the child node of the neighbor at x=-1,y=-1,z=0
  int box_idxNbr_i0_j0_km1_ch;   // Box index of the child node of the neighbor at x=0,y=0,z=-1
  int box_idxNbr_im1_j0_km1_ch;  // Box index of the child node of the neighbor at x=-1,y=0,z=-1
  int box_idxNbr_i0_jm1_km1_ch;  // Box index of the child node of the neighbor at x=0,y=-1,z=-1
  int box_idxNbr_im1_jm1_km1_ch; // Box index of the child node of the neighbor at x=-1,y=-1,z=-1

  int box_idx_ch; // Box index of the child node

  int box_grid_idx_ch; // Box grid index of the child node

  bool is_bder_grid_point; // Ask if the grid point is interior
  bool grid_point_exist;

  int box_real_dim_X_ch;
  int box_real_dim_X_times_Y_ch;

  int grid_box_real_dim_X_ch;
  int grid_box_real_dim_X_times_Y_ch;

  int size_bder_grid_points;
  int size_intr_grid_points;
  int size_SIMULATION_BOUNDARY_grid_points;

  bool grid_point_simulation_boundary;

  int aux_idx_x;
  int aux_idx_y;
  int aux_idx_z;

  int lv = ptr_node->lv;

  for (int ch = 0; ch < ptr_node->zones_size; ch++)
  {
    ptr_ch = ptr_node->pptr_chn[ch];

    box_real_dim_X_ch = ptr_ch->box_real_dim_x;
    box_real_dim_X_times_Y_ch = ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;

    grid_box_real_dim_X_ch = ptr_ch->box_real_dim_x + 1;
    grid_box_real_dim_X_times_Y_ch = (ptr_ch->box_real_dim_x + 1) * (ptr_ch->box_real_dim_y + 1);

    // 27 * N - 9 * (N-1) = 18N+9 >= Total of grid points; N = number of cubics with 8 cells and 27 grid points
    // Maximum of boundary grid points = 9*(N*2+1) - (N*2-1) = 16N+10 (Imagine a row of cubics)
    // Maximum of interior grid points = (2k-1)^3 < 8k^3, where k^3 = N
    size_bder_grid_points = ptr_ch->cell_size / 8 * 16 + 10;
    size_intr_grid_points = ptr_ch->cell_size;
    size_SIMULATION_BOUNDARY_grid_points = size_bder_grid_points;

    //* >> Space checking of border grid points array*//
    if (space_check(&(ptr_ch->grid_bdry_cap), size_bder_grid_points, 1.0f, "p4i1i1i1i1", &(ptr_ch->ptr_bdry_grid_cell_idx_x), &(ptr_ch->ptr_bdry_grid_cell_idx_y), &(ptr_ch->ptr_bdry_grid_cell_idx_z), &(ptr_ch->ptr_bdry_box_grid_idx)) == _FAILURE_)
    {
      printf("Error, in space_check function\n");
      return _FAILURE_;
    }

    //* >> Space checking of interior grid points array*//
    if (space_check(&(ptr_ch->grid_intr_cap), size_intr_grid_points, 1.0f, "p4i1i1i1i1", &(ptr_ch->ptr_intr_grid_cell_idx_x), &(ptr_ch->ptr_intr_grid_cell_idx_y), &(ptr_ch->ptr_intr_grid_cell_idx_z), &(ptr_ch->ptr_intr_box_grid_idx)) == _FAILURE_)
    {
      printf("Error, in space_check function\n");
      return _FAILURE_;
    }

    if (ptr_ch->sim_bdry_contact == true)
    {
      //* >> Space checking of simulation boundary grid points array*//
      if (space_check(&(ptr_ch->grid_sim_bdry_cap), size_SIMULATION_BOUNDARY_grid_points, 1.0f, "p4i1i1i1i1", &(ptr_ch->ptr_sim_bdry_grid_cell_idx_x), &(ptr_ch->ptr_sim_bdry_grid_cell_idx_y), &(ptr_ch->ptr_sim_bdry_grid_cell_idx_z), &(ptr_ch->ptr_sim_bdry_box_grid_idx)) == _FAILURE_)
      {
        printf("Error, in space_check function\n");
        return _FAILURE_;
      }
    }

    //* >> Grid points *//
    for (int kk = ptr_ch->box_min_z - ptr_ch->box_ts_z; kk < ptr_ch->box_max_z - ptr_ch->box_ts_z + 2; kk++)
    {
      for (int jj = ptr_ch->box_min_y - ptr_ch->box_ts_y; jj < ptr_ch->box_max_y - ptr_ch->box_ts_y + 2; jj++)
      {
        for (int ii = ptr_ch->box_min_x - ptr_ch->box_ts_x; ii < ptr_ch->box_max_x - ptr_ch->box_ts_x + 2; ii++)
        {
          box_idx_ch = ii + jj * ptr_ch->box_real_dim_x + kk * ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;

          grid_point_exist = false;
          grid_point_simulation_boundary = false;

          box_idxNbr_i0_j0_k0_ch = box_idx_ch;
          box_idxNbr_im1_j0_k0_ch = box_idx_ch - 1;
          box_idxNbr_i0_jm1_k0_ch = box_idx_ch - box_real_dim_X_ch;
          box_idxNbr_im1_jm1_k0_ch = box_idx_ch - 1 - box_real_dim_X_ch;
          box_idxNbr_i0_j0_km1_ch = box_idx_ch - box_real_dim_X_times_Y_ch;
          box_idxNbr_im1_j0_km1_ch = box_idx_ch - 1 - box_real_dim_X_times_Y_ch;
          box_idxNbr_i0_jm1_km1_ch = box_idx_ch - box_real_dim_X_ch - box_real_dim_X_times_Y_ch;
          box_idxNbr_im1_jm1_km1_ch = box_idx_ch - 1 - box_real_dim_X_ch - box_real_dim_X_times_Y_ch;

          // First
          if (ptr_ch->ptr_box[box_idxNbr_i0_j0_k0_ch] > -4)
          {
            grid_point_exist = true;
            is_bder_grid_point = false;
            //* Connection to the right  *//

            //* Connection to the left  *//
            if (ptr_ch->ptr_box[box_idxNbr_im1_j0_k0_ch] < -3 &&
                ptr_ch->ptr_box[box_idxNbr_im1_jm1_k0_ch] < -3 &&
                ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] < -3 &&
                ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] < -3)
            {
              is_bder_grid_point = true;
            }
            //* Backward connection   *//

            //* Forward connection   *//
            else if (ptr_ch->ptr_box[box_idxNbr_i0_jm1_k0_ch] < -3 &&
                     ptr_ch->ptr_box[box_idxNbr_im1_jm1_k0_ch] < -3 &&
                     ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] < -3 &&
                     ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] < -3)
            {
              is_bder_grid_point = true;
            }
            //* Upward connection *//

            //* Down connection *//
            else if (ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] < -3 &&
                     ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] < -3 &&
                     ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] < -3 &&
                     ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] < -3)
            {
              is_bder_grid_point = true;
            }
          }
          // Second
          else if (ptr_ch->ptr_box[box_idxNbr_im1_j0_k0_ch] > -4)
          {
            grid_point_exist = true;
            is_bder_grid_point = false;
            //* Connection to the right  *//
            if (ptr_ch->ptr_box[box_idxNbr_i0_jm1_k0_ch] < -3 &&
                ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] < -3 &&
                ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] < -3)
            {
              is_bder_grid_point = true;
            }
            //* Connection to the left  *//

            //* Backward connection   *//

            //* Forward connection   *//
            else if (ptr_ch->ptr_box[box_idxNbr_i0_jm1_k0_ch] < -3 &&
                     ptr_ch->ptr_box[box_idxNbr_im1_jm1_k0_ch] < -3 &&
                     ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] < -3 &&
                     ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] < -3)
            {
              is_bder_grid_point = true;
            }
            //* Upward connection *//

            //* Down connection *//
            else if (ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] < -3 &&
                     ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] < -3 &&
                     ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] < -3 &&
                     ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] < -3)
            {
              is_bder_grid_point = true;
            }
          }
          // Third
          else if (ptr_ch->ptr_box[box_idxNbr_i0_jm1_k0_ch] > -4)
          {
            grid_point_exist = true;
            is_bder_grid_point = false;
            //* Connection to the right  *//

            //* Connection to the left  *//
            if (ptr_ch->ptr_box[box_idxNbr_im1_jm1_k0_ch] < -3 &&
                ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] < -3 &&
                ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] < -3)
            {
              is_bder_grid_point = true;
            }
            //* Backward connection   *//
            else if (ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] < -3 &&
                     ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] < -3)
            {
              is_bder_grid_point = true;
            }
            //* Forward connection   *//

            //* Upward connection *//

            //* Down connection *//
            else if (ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] < -3 &&
                     ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] < -3 &&
                     ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] < -3 &&
                     ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] < -3)
            {
              is_bder_grid_point = true;
            }
          }
          // Fourth
          else if (ptr_ch->ptr_box[box_idxNbr_im1_jm1_k0_ch] > -4)
          {
            grid_point_exist = true;
            is_bder_grid_point = false;
            //* Connection to the right  *//
            if (ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] < -3 &&
                ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] < -3)
            {
              is_bder_grid_point = true;
            }
            //* Connection to the left  *//

            //* Backward connection   *//
            else if (ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] < -3 &&
                     ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] < -3)
            {
              is_bder_grid_point = true;
            }
            //* Forward connection   *//

            //* Upward connection *//

            //* Down connection *//
            else if (ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] < -3 &&
                     ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] < -3 &&
                     ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] < -3 &&
                     ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] < -3)
            {
              is_bder_grid_point = true;
            }
          }
          // Others
          else if (ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] > -4 ||
                   ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] > -4 ||
                   ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] > -4 ||
                   ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] > -4)
          {
            grid_point_exist = true;
            is_bder_grid_point = true;
          }

          if (ptr_ch->sim_bdry_contact == true && grid_point_exist == true)
          {
            if (bdry_cond_type == 0 || (bdry_cond_type != 0 && is_bder_grid_point == true))
            {
              aux_idx_x = ii + ptr_ch->box_ts_x;
              aux_idx_y = jj + ptr_ch->box_ts_y;
              aux_idx_z = kk + ptr_ch->box_ts_z;

              // Notes that here is not necessary to do the transformation for grid points with negative indices values
              // The reason is because if it is negative, then, it never will be a boundary simulation box grid point.

              if (aux_idx_x == 0 || aux_idx_x == (1 << (lv + 1)) ||
                  aux_idx_y == 0 || aux_idx_y == (1 << (lv + 1)) ||
                  aux_idx_z == 0 || aux_idx_z == (1 << (lv + 1)))
              {
                grid_point_simulation_boundary = true;
              }
            }
          }

          if (grid_point_exist == true)
          {

            aux_idx_x = ii + ptr_ch->box_ts_x < 0 ? ii + ptr_ch->box_ts_x + (1 << (lv + 1)) : ii + ptr_ch->box_ts_x;
            aux_idx_y = jj + ptr_ch->box_ts_y < 0 ? jj + ptr_ch->box_ts_y + (1 << (lv + 1)) : jj + ptr_ch->box_ts_y;
            aux_idx_z = kk + ptr_ch->box_ts_z < 0 ? kk + ptr_ch->box_ts_z + (1 << (lv + 1)) : kk + ptr_ch->box_ts_z;

            box_grid_idx_ch = ii + jj * grid_box_real_dim_X_ch + kk * grid_box_real_dim_X_times_Y_ch;

            //* >> Adding the grid point *//

            //* >> Simulation Boundary grid point*//
            if (grid_point_simulation_boundary == true)
            {
              ptr_ch->ptr_sim_bdry_grid_cell_idx_x[ptr_ch->grid_sim_bdry_size] = aux_idx_x;
              ptr_ch->ptr_sim_bdry_grid_cell_idx_y[ptr_ch->grid_sim_bdry_size] = aux_idx_y;
              ptr_ch->ptr_sim_bdry_grid_cell_idx_z[ptr_ch->grid_sim_bdry_size] = aux_idx_z;
              ptr_ch->ptr_sim_bdry_box_grid_idx[ptr_ch->grid_sim_bdry_size] = box_grid_idx_ch;
              ptr_ch->grid_sim_bdry_size += 1; // Increasing the number of boundary grid points in the array
            }
            //* >> Boundary grid point*//
            else if (is_bder_grid_point == true)
            {
              //* >> Adding the grid point to the Boundary array *//
              ptr_ch->ptr_bdry_grid_cell_idx_x[ptr_ch->grid_bdry_size] = aux_idx_x;
              ptr_ch->ptr_bdry_grid_cell_idx_y[ptr_ch->grid_bdry_size] = aux_idx_y;
              ptr_ch->ptr_bdry_grid_cell_idx_z[ptr_ch->grid_bdry_size] = aux_idx_z;
              ptr_ch->ptr_bdry_box_grid_idx[ptr_ch->grid_bdry_size] = box_grid_idx_ch;
              ptr_ch->grid_bdry_size += 1; // Increasing the number of Boundary grid points in the array
            }
            //* Interior grid point *//
            else
            {

              //* >> Adding the grid point to the interior array *//
              ptr_ch->ptr_intr_grid_cell_idx_x[ptr_ch->grid_intr_size] = aux_idx_x;
              ptr_ch->ptr_intr_grid_cell_idx_y[ptr_ch->grid_intr_size] = aux_idx_y;
              ptr_ch->ptr_intr_grid_cell_idx_z[ptr_ch->grid_intr_size] = aux_idx_z;
              ptr_ch->ptr_intr_box_grid_idx[ptr_ch->grid_intr_size] = box_grid_idx_ch;
              ptr_ch->grid_intr_size += 1; // Increasing the number of interior grid points in the array
            }
          }
        }
      }
    }
    // for (int kk = ptr_ch->box_min_z - ptr_ch->box_ts_z; kk < ptr_ch->box_max_z - ptr_ch->box_ts_z + 2; kk++)
    // {
    //     for (int jj = ptr_ch->box_min_y - ptr_ch->box_ts_y; jj < ptr_ch->box_max_y - ptr_ch->box_ts_y + 2; jj++)
    //     {
    //         for (int ii = ptr_ch->box_min_x - ptr_ch->box_ts_x; ii < ptr_ch->box_max_x - ptr_ch->box_ts_x + 2; ii++)
    //         {
    //             box_idx_ch = ii + jj * box_real_dim_X_ch + kk * box_real_dim_X_times_Y_ch;

    //             grid_point_exist = false;

    //             box_idxNbr_i0_j0_k0_ch = box_idx_ch;
    //             box_idxNbr_im1_j0_k0_ch = box_idx_ch - 1;
    //             box_idxNbr_i0_jm1_k0_ch = box_idx_ch - box_real_dim_X_ch;
    //             box_idxNbr_im1_jm1_k0_ch = box_idx_ch - 1 - box_real_dim_X_ch;
    //             box_idxNbr_i0_j0_km1_ch = box_idx_ch - box_real_dim_X_times_Y_ch;
    //             box_idxNbr_im1_j0_km1_ch = box_idx_ch - 1 - box_real_dim_X_times_Y_ch;
    //             box_idxNbr_i0_jm1_km1_ch = box_idx_ch - box_real_dim_X_ch - box_real_dim_X_times_Y_ch;
    //             box_idxNbr_im1_jm1_km1_ch = box_idx_ch - 1 - box_real_dim_X_ch - box_real_dim_X_times_Y_ch;

    //             //First
    //             if(ptr_ch->ptr_box[box_idxNbr_i0_j0_k0_ch] > -4)
    //             {
    //                 grid_point_exist = true;
    //                 is_bder_grid_point = false;
    //                 //* Connection to the right  *//

    //                 //* Connection to the left  *//
    //                 if (ptr_ch->ptr_box[box_idxNbr_im1_j0_k0_ch] < -3 &&
    //                          ptr_ch->ptr_box[box_idxNbr_im1_jm1_k0_ch] < -3 &&
    //                          ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] < -3 &&
    //                          ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] < -3)
    //                 {
    //                     is_bder_grid_point = true;
    //                 }
    //                 //* Backward connection   *//

    //                 //* Forward connection   *//
    //                 else if (ptr_ch->ptr_box[box_idxNbr_i0_jm1_k0_ch] < -3 &&
    //                          ptr_ch->ptr_box[box_idxNbr_im1_jm1_k0_ch] < -3 &&
    //                          ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] < -3 &&
    //                          ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] < -3)
    //                 {
    //                     is_bder_grid_point = true;
    //                 }
    //                 //* Upward connection *//

    //                 //* Down connection *//
    //                 else if (ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] < -3 &&
    //                          ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] < -3 &&
    //                          ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] < -3 &&
    //                          ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] < -3)
    //                 {
    //                     is_bder_grid_point = true;
    //                 }
    //             }
    //             // Second
    //             else if (ptr_ch->ptr_box[box_idxNbr_im1_j0_k0_ch] > -4)
    //             {
    //                 grid_point_exist = true;
    //                 is_bder_grid_point = false;
    //                 //* Connection to the right  *//
    //                 if (ptr_ch->ptr_box[box_idxNbr_i0_jm1_k0_ch] < -3 &&
    //                     ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] < -3 &&
    //                     ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] < -3)
    //                 {
    //                     is_bder_grid_point = true;
    //                 }
    //                 //* Connection to the left  *//

    //                 //* Backward connection   *//

    //                 //* Forward connection   *//
    //                 else if (ptr_ch->ptr_box[box_idxNbr_i0_jm1_k0_ch] < -3 &&
    //                          ptr_ch->ptr_box[box_idxNbr_im1_jm1_k0_ch] < -3 &&
    //                          ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] < -3 &&
    //                          ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] < -3)
    //                 {
    //                     is_bder_grid_point = true;
    //                 }
    //                 //* Upward connection *//

    //                 //* Down connection *//
    //                 else if (ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] < -3 &&
    //                          ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] < -3 &&
    //                          ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] < -3 &&
    //                          ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] < -3)
    //                 {
    //                     is_bder_grid_point = true;
    //                 }
    //             }
    //             // Third
    //             else if(ptr_ch->ptr_box[box_idxNbr_i0_jm1_k0_ch] > -4)
    //             {
    //                 grid_point_exist = true;
    //                 is_bder_grid_point = false;
    //                 //* Connection to the right  *//

    //                 //* Connection to the left  *//
    //                 if (ptr_ch->ptr_box[box_idxNbr_im1_jm1_k0_ch] < -3 &&
    //                          ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] < -3 &&
    //                          ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] < -3)
    //                 {
    //                     is_bder_grid_point = true;
    //                 }
    //                 //* Backward connection   *//
    //                 else if (ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] < -3 &&
    //                          ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] < -3)
    //                 {
    //                     is_bder_grid_point = true;
    //                 }
    //                 //* Forward connection   *//

    //                 //* Upward connection *//

    //                 //* Down connection *//
    //                 else if (ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] < -3 &&
    //                          ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] < -3 &&
    //                          ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] < -3 &&
    //                          ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] < -3)
    //                 {
    //                     is_bder_grid_point = true;
    //                 }
    //             }
    //             // Fourth
    //             else if (ptr_ch->ptr_box[box_idxNbr_im1_jm1_k0_ch] > -4)
    //             {
    //                 grid_point_exist = true;
    //                 is_bder_grid_point = false;
    //                 //* Connection to the right  *//
    //                 if (ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] < -3 &&
    //                     ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] < -3)
    //                 {
    //                     is_bder_grid_point = true;
    //                 }
    //                 //* Connection to the left  *//

    //                 //* Backward connection   *//
    //                 else if (ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] < -3 &&
    //                          ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] < -3)
    //                 {
    //                     is_bder_grid_point = true;
    //                 }
    //                 //* Forward connection   *//

    //                 //* Upward connection *//

    //                 //* Down connection *//
    //                 else if (ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] < -3 &&
    //                          ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] < -3 &&
    //                          ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] < -3 &&
    //                          ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] < -3)
    //                 {
    //                     is_bder_grid_point = true;
    //                 }
    //             }
    //             // Others
    //             else if (ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] > -4 ||
    //                      ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] > -4 ||
    //                      ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] > -4 ||
    //                      ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] > -4)
    //             {
    //                 grid_point_exist = true;
    //                 is_bder_grid_point = true;
    //             }

    //             // //* >> The grid point exist *//
    //             // if (ptr_ch->ptr_box[box_idxNbr_i0_j0_k0_ch] > -4 ||
    //             //     ptr_ch->ptr_box[box_idxNbr_im1_j0_k0_ch] > -4 ||
    //             //     ptr_ch->ptr_box[box_idxNbr_i0_jm1_k0_ch] > -4 ||
    //             //     ptr_ch->ptr_box[box_idxNbr_im1_jm1_k0_ch] > -4 ||
    //             //     ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] > -4 ||
    //             //     ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] > -4 ||
    //             //     ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] > -4 ||
    //             //     ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] > -4)
    //             // {
    //             //     grid_point_exist = true;
    //             //     is_bder_grid_point = false;
    //             //     //* Connection to the right  *//
    //             //     if (ptr_ch->ptr_box[box_idxNbr_i0_j0_k0_ch] < -3 &&
    //             //         ptr_ch->ptr_box[box_idxNbr_i0_jm1_k0_ch] < -3 &&
    //             //         ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] < -3 &&
    //             //         ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] < -3)
    //             //     {
    //             //         is_bder_grid_point = true;
    //             //     }
    //             //     //* Connection to the left  *//
    //             //     else if (ptr_ch->ptr_box[box_idxNbr_im1_j0_k0_ch] < -3 &&
    //             //              ptr_ch->ptr_box[box_idxNbr_im1_jm1_k0_ch] < -3 &&
    //             //              ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] < -3 &&
    //             //              ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] < -3)
    //             //     {
    //             //         is_bder_grid_point = true;
    //             //     }
    //             //     //* Backward connection   *//
    //             //     else if (ptr_ch->ptr_box[box_idxNbr_i0_j0_k0_ch] < -3 &&
    //             //              ptr_ch->ptr_box[box_idxNbr_im1_j0_k0_ch] < -3 &&
    //             //              ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] < -3 &&
    //             //              ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] < -3)
    //             //     {
    //             //         is_bder_grid_point = true;
    //             //     }
    //             //     //* Forward connection   *//
    //             //     else if (ptr_ch->ptr_box[box_idxNbr_i0_jm1_k0_ch] < -3 &&
    //             //              ptr_ch->ptr_box[box_idxNbr_im1_jm1_k0_ch] < -3 &&
    //             //              ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] < -3 &&
    //             //              ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] < -3)
    //             //     {
    //             //         is_bder_grid_point = true;
    //             //     }
    //             //     //* Upward connection *//
    //             //     else if (ptr_ch->ptr_box[box_idxNbr_i0_j0_k0_ch] < -3 &&
    //             //              ptr_ch->ptr_box[box_idxNbr_im1_j0_k0_ch] < -3 &&
    //             //              ptr_ch->ptr_box[box_idxNbr_i0_jm1_k0_ch] < -3 &&
    //             //              ptr_ch->ptr_box[box_idxNbr_im1_jm1_k0_ch] < -3)
    //             //     {
    //             //         is_bder_grid_point = true;
    //             //     }
    //             //     //* Down connection *//
    //             //     else if (ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] < -3 &&
    //             //              ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] < -3 &&
    //             //              ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] < -3 &&
    //             //              ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] < -3)
    //             //     {
    //             //         is_bder_grid_point = true;
    //             //     }
    //             // }

    //             if (grid_point_exist == true)
    //             {
    //                 box_grid_idx_ch = ii + jj * grid_box_real_dim_X_ch + kk * grid_box_real_dim_X_times_Y_ch;
    //                 //* >> Adding the grid point *//

    //                 //* >> Border grid point*//
    //                 if (is_bder_grid_point == true)
    //                 {
    //                     //* >> Adding the grid point to the border array *//
    //                     ptr_ch->ptr_bdry_grid_cell_idx_x[ptr_ch->grid_bdry_size] = ii + ptr_ch->box_ts_x;
    //                     ptr_ch->ptr_bdry_grid_cell_idx_y[ptr_ch->grid_bdry_size] = jj + ptr_ch->box_ts_y;
    //                     ptr_ch->ptr_bdry_grid_cell_idx_z[ptr_ch->grid_bdry_size] = kk + ptr_ch->box_ts_z;
    //                     ptr_ch->ptr_bdry_box_grid_idx[ptr_ch->grid_bdry_size] = box_grid_idx_ch;
    //                     ptr_ch->grid_bdry_size += 1; // Increasing the number of border grid points in the array
    //                 }
    //                 //* Interior grid point *//
    //                 else
    //                 {
    //                     //* >> Adding the grid point to the interior array *//
    //                     ptr_ch->ptr_intr_grid_cell_idx_x[ptr_ch->grid_intr_size] = ii + ptr_ch->box_ts_x;
    //                     ptr_ch->ptr_intr_grid_cell_idx_y[ptr_ch->grid_intr_size] = jj + ptr_ch->box_ts_y;
    //                     ptr_ch->ptr_intr_grid_cell_idx_z[ptr_ch->grid_intr_size] = kk + ptr_ch->box_ts_z;
    //                     ptr_ch->ptr_intr_box_grid_idx[ptr_ch->grid_intr_size] = box_grid_idx_ch;
    //                     ptr_ch->grid_intr_size += 1; // Increasing the number of interior grid points in the array
    //                 }
    //             }
    //         }
    //     }
    // }
  }
  return _SUCCESS_;
}

/**
 * @brief Computing the number of particles outside of the refinement zones of
 * the node used.
 *
 * **SHORT DESCRIPTION**: Computing the number of particles outside of the
 * refinement zones of the node used. This number is stored at the node 
 * parameter \ref node.no_ptcl_outs_ref_zones "no_ptcl_outs_ref_zones".
 *
 * **PREREQUISITES**: If there is at least a new refinement zone, or if the node
 * used is the Head node.
 *
 * @param[in,out] ptr_node Pointer to node structure
 *
 * **RETURN**: The error status.
 *
 * **LONG DESCRIPTION**:
 *
 * Computing the number of particles outside of the refinement zones of the node 
 * used. This number is stored at the node parameter \ref 
 * node.no_ptcl_outs_ref_zones "no_ptcl_outs_ref_zones".
 *
 * This function is initially called by the "for" loop over the level of 
 * refinement at the function tree_adaptation() with the child nodes as input, 
 * and also it is called one time at the end of the tree_adaptation() function 
 * with the Head node as the input.
 *
 * The flux of this function can be seen in the figure (work in progress), and
 * it is explained below:
 *
 * - [0]  <b> THE computing_no_ptcl_outside_refinement_zones() FUNCTION STARTS....</b>
 *
 * - [1]  Defining some internal useful parameters. The value of \ref
 * node.no_ptcl_outs_ref_zones "no_ptcl_outs_ref_zones" is initialized.
 *
 * - [2]  Run a "for" loop over the child nodes.
 *
 * - [3]  Per every child node, the value of the \ref 
 * node.no_ptcl_outs_ref_zones "no_ptcl_outs_ref_zones" is initialized is 
 * updated.
 *
 * - [6]  <b> THE computing_no_ptcl_outside_refinement_zones() FUNCTION ENDS....</b>
 *
 * **ILLUSTRATIVE EXAMPLES**:
 * - [a]  a
 *
 * **RATIONALES**:
 * - [a]  aaa
 *
 * **NOTES**:
 * - [a]  a
 */

static void computing_no_ptcl_outside_refinement_zones(struct node *ptr_node)
{
  // The input can be a child node or the head node

  struct node *ptr_ch;

  ptr_node->no_ptcl_outs_ref_zones = ptr_node->no_ptcl_full_node;

  for (int i = 0; i < ptr_node->chn_size; i++)
  {
    ptr_ch = ptr_node->pptr_chn[i];
    ptr_node->no_ptcl_outs_ref_zones -= ptr_ch->no_ptcl_full_node;
  }
}

/**
 * @brief Transferring unused child nodes to the memory pool
 *
 * **SHORT DESCRIPTION**: Transferring unused child nodes to the stack of memory 
 * pool to be reused in the future by any other parent node at any level of 
 * refinement.
 *
 * **PREREQUISITES**: Always used.
 *
 * @param[in,out] ptr_node Pointer to node structure
 *
 * **RETURN**: No parameter is returned.
 *
 * **LONG DESCRIPTION**:
 *
 * Transferring unused child nodes to the stack of memory pool to be reused in 
 * the future by any other parent node at any level of refinement.
 *
 * Moreover, a recent implementation lets remove refinement zones of the parent 
 * node which contains fewer particles than the minimum particles required to 
 * refine a cell. So, this zone is removed, and the corresponding child node is 
 * also removed and transferred to the stack of the memory pool. The gap in the 
 * parent node is filled by the last node in the list of child nodes, and also 
 * this surrogate child node is updated in its parameters.
 *
 * The flux of this function can be seen in the figure (work in progress), and 
 * it is explained below:
 *
 * - [0]  <b> THE transferring_unused_child_node_to_memory_pool() FUNCTION STARTS....</b>
 *
 * - [1]  Defining some internal useful parameters.
 *
 * - [2]  Run a "for" loop over the excess between the number of refinement 
 * zones and the total number of child nodes including linked and unlinked nodes
 *
 * - [3]  The nodes are removed. Remember that at this point the child nodes are
 * sorted in the parent parameter array \ref node.pptr_chn "pptr_chn" such that
 * the first elements correspond exactly with the first refinement zones, so the
 * excess of child nodes in that parent parameter array corresponds to the nodes
 * that will not be used by the parent node.
 *
 * - [4]  The second part of the function consists in to remove zones and child
 * which contains fewer particles than the minimum required to refine a cell. To
 * do that, a "for" loop over the new refinement zones is performed.
 *
 * - [5]  For every new refinement zone, if the child node corresponding has 
 * fewer particles than the refinement cell criteria, then it will be replaced 
 * by the last node in the list.
 *
 * - [6]  The surrogating node is updated in its identification, and the box 
 * cells of the parent node are updated with the information of the surrogating 
 * node and the removed node.
 *
 * - [7]  The removed node is sent to the stack of the memory pool.
 *
 * - [8]  Finally, some other important parameters are updated.
 *
 * - [9]  <b> THE transferring_unused_child_node_to_memory_pool() FUNCTION ENDS....</b>
 *
 * **ILLUSTRATIVE EXAMPLES**:
 * - [a]  a
 *
 * **RATIONALES**:
 * - [a]  aaa
 *
 * **NOTES**:
 * - [a]  a
 */

static void transferring_unused_child_node_to_memory_pool(struct node *ptr_node)
{
  int box_idx_node;

  int aux_int;
  int *aux_ptr_int;

  //* >> Putting unused child nodes to the stack of memory pool *//
  for (int i = ptr_node->zones_size; i < ptr_node->chn_size; i++)
  {
    add_node_to_stack(ptr_node->pptr_chn[i]);
  }

  //* >> Removing zones (children nodes) with 0 particles
  for (int zone_idx = 0; zone_idx < ptr_node->zones_size; zone_idx++)
  {
    if (ptr_node->pptr_chn[zone_idx]->no_ptcl_full_node < ref_criterion_ptcl)
    {
      //* >> Changing the ID of the replacement zone
      ptr_node->pptr_chn[ptr_node->zones_size - 1]->ID = zone_idx;

      //* >> Changing the box status of the replacement zone from the old zone value "ptr_node->zones_size - 1" to the new one = ch
      for (int cell_idx = 0; cell_idx < ptr_node->ptr_zone_size[ptr_node->zones_size - 1]; cell_idx++)
      {
        box_idx_node = ptr_node->ptr_box_idx[ptr_node->pptr_zones[ptr_node->zones_size - 1][cell_idx]];
        ptr_node->ptr_box[box_idx_node] = zone_idx;
      }

      //* >> Changing the box status of the removed zone from the old zone value ch to the new one = -3
      for (int cell_idx = 0; cell_idx < ptr_node->ptr_zone_size[zone_idx]; cell_idx++)
      {
        box_idx_node = ptr_node->ptr_box_idx[ptr_node->pptr_zones[zone_idx][cell_idx]];
        ptr_node->ptr_box[box_idx_node] = -3;
      }

      //* >> Removing the zone ch without particles from the children of the node ptr_node
      add_node_to_stack(ptr_node->pptr_chn[zone_idx]);

      //* >> Replacing the ch values of the last children
      ptr_node->pptr_chn[zone_idx] = ptr_node->pptr_chn[ptr_node->zones_size - 1];
      ptr_node->pptr_chn[ptr_node->zones_size - 1] = NULL;

      //* >> Updating the zone values
      // Size of the zones
      aux_int = ptr_node->ptr_zone_size[zone_idx];
      ptr_node->ptr_zone_size[zone_idx] = ptr_node->ptr_zone_size[ptr_node->zones_size - 1];
      ptr_node->ptr_zone_size[ptr_node->zones_size - 1] = aux_int;
      // Cap of the zones
      aux_int = ptr_node->ptr_zone_cap[zone_idx];
      ptr_node->ptr_zone_cap[zone_idx] = ptr_node->ptr_zone_cap[ptr_node->zones_size - 1];
      ptr_node->ptr_zone_cap[ptr_node->zones_size - 1] = aux_int;
      // Zone array
      aux_ptr_int = ptr_node->pptr_zones[zone_idx];
      ptr_node->pptr_zones[zone_idx] = ptr_node->pptr_zones[ptr_node->zones_size - 1];
      ptr_node->pptr_zones[ptr_node->zones_size - 1] = aux_ptr_int;

      // Total number of active zones
      ptr_node->zones_size -= 1; //-1 the the total number of children in the node
                                 // ptr_node->zones_size -= 1;	//-1 the the total number of zones in the node

      zone_idx--;
    }
  }
}

/**
 * @brief Updating of the \ref tree_adaptation__REMINDER__Tentacles "Tentacles" 
 * at the child nodes level of refinement
 *
 * **SHORT DESCRIPTION**: Updating of the \ref 
 * tree_adaptation__REMINDER__Tentacles "Tentacles" at the child nodes level of 
 * refinement. tree_adaptation() function reset the size of the \ref 
 * tree_adaptation__REMINDER__Tentacles "Tentacles" at every level of 
 * refinement, and call this function at every parent node of that level of 
 * refinement but using the level of the child nodes as input. So, in the end, 
 * the \ref tree_adaptation__REMINDER__Tentacles "Tentacles" update accordingly 
 * the new configuration of child nodes.
 *
 * **PREREQUISITES**: If there is at least a new refinement zone.
 *
 * @param[in] ptr_node Pointer to node structure
 * @param[in] tentacle_lv Level of refinement of the \ref 
 * tree_adaptation__REMINDER__Tentacles "Tentacles"
 *
 * **RETURN**: The error status.
 *
 * **LONG DESCRIPTION**:
 *
 * Updating of the \ref tree_adaptation__REMINDER__Tentacles "Tentacles" at the 
 * child nodes level of refinement. tree_adaptation() function reset the size of 
 * the \ref tree_adaptation__REMINDER__Tentacles "Tentacles" at every level of 
 * refinement, and call this function at every parent node of that level of 
 * refinement but using the level of the child nodes as input. So, in the end, 
 * the \ref tree_adaptation__REMINDER__Tentacles "Tentacles" update accordingly 
 * the new configuration of child nodes.
 *
 * Here the only global parameters associated with the \ref 
 * tree_adaptation__REMINDER__Tentacles "Tentacles" are going to be updated.
 *
 * The flux of this function can be seen in the figure (work in progress), and
 * it is explained below:
 *
 * - [0]  <b> THE tentacles_updating() FUNCTION STARTS....</b>
 *
 * - [1]  Defining some internal useful parameters.
 *
 * - [2]  Run a "for" loop over the number of refinement zones, that at this
 * point of the module tree_adaptation.c is equal to the final number of child
 * nodes.
 *
 * - [3]  The child node is pointed by the tentacle at the corresponding 
 * position and level of refinement.
 *
 * - [4]  Finally, the tentacles Size (see Key Concepts \ref Key_Concepts_Size 
 * "Size") of the global parameter \link GL_tentacles_size \endlink is updated.
 *
 * - [5]  <b> THE tentacles_updating() FUNCTION ENDS....</b>
 *
 * **ILLUSTRATIVE EXAMPLES**:
 * - [a]  a
 *
 * **RATIONALES**:
 * - [a]  aaa
 *
 * **NOTES**:
 * - [a]  a
 */

static int tentacles_updating(const struct node *ptr_node, int tentacle_lv)
{
  int no_tentacles = GL_tentacles_size[tentacle_lv];
  int size = no_tentacles + ptr_node->zones_size;
  int no_zones = ptr_node->zones_size;

  //* >> Space checking of the capacity of the refined cells *//
  if (space_check(&(GL_tentacles_cap[tentacle_lv]), size, 4.0f, "p1n2", &(GL_tentacles[tentacle_lv])) == _FAILURE_)
  {
    printf("Error, in space_check function\n");
    return _FAILURE_;
  }

  for (int i = 0; i < no_zones; i++)
  {
    GL_tentacles[tentacle_lv][no_tentacles + i] = ptr_node->pptr_chn[i];
  }

  GL_tentacles_size[tentacle_lv] = size;
  return _SUCCESS_;
}

/**
 * @brief Updating the maximum level parameter \link GL_tentacles_level_max 
 * \endlink of the \ref tree_adaptation__REMINDER__Tentacles "Tentacles"
 *
 * **SHORT DESCRIPTION**: Updating the maximum level parameter \link 
 * GL_tentacles_level_max \endlink of the \ref 
 * tree_adaptation__REMINDER__Tentacles "Tentacles".
 *
 * **PREREQUISITES**: Always used.
 *
 * **RETURN**: No parameter is returned.
 *
 * **LONG DESCRIPTION**:
 *
 * Updating the maximum level parameter \link GL_tentacles_level_max \endlink of 
 * the \ref tree_adaptation__REMINDER__Tentacles "Tentacles". this function is 
 * used only once for each call of function tree_adaptation().
 *
 * The flux of this function can be seen in the figure (work in progress), and
 * it is explained below:
 *
 * - [0]  <b> THE updating_tentacles_max_lv() FUNCTION STARTS....</b>
 *
 * - [1]  Defining some internal useful parameters.
 *
 * - [2]  Run a "for" loop over the refinement levels but avoiding the finest
 * refinement level. The loop goes from a finer refinement level to a coarser
 * one.
 *
 * - [3]  The child node is pointed by the tentacle at the corresponding
 * position and level of refinement.
 *
 * - [4]  If the number of  \ref tree_adaptation__REMINDER__Tentacles 
 * "Tentacles" at that level is bigger than cero, then the maimum level 
 * parameter \link GL_tentacles_level_max \endlink is obtained
 *
 * - [5]  <b> THE updating_tentacles_max_lv() FUNCTION ENDS....</b>
 *
 * **ILLUSTRATIVE EXAMPLES**:
 * - [a]  a
 *
 * **RATIONALES**:
 * - [a]  aaa
 *
 * **NOTES**:
 * - [a]  a
 */

static void updating_tentacles_max_lv(void)
{
  int no_lvs = GL_tentacles_level_max < (lmax - lmin) ? GL_tentacles_level_max : (GL_tentacles_level_max - 1);

  for (int lv = no_lvs + 1; lv > -1; lv--)
  {
    if (GL_tentacles_size[lv] > 0)
    {
      GL_tentacles_level_max = lv;
      lv = -1;
    }
  }
}

/**
 * @brief Make the calls of all the local functions of the tree_adaptation.c 
 * module
 *
 * **SHORT DESCRIPTION**: Make the calls of all the local functions of the 
 * tree_adaptation.c module
 *
 * **PREREQUISITES**: Always used.
 *
 * **RETURN**: The error status.
 *
 * **ILLUSTRATIVE EXAMPLES**:
 * - [a]  a
 *
 * **RATIONALES**:
 * - [a]  aaa
 *
 * **NOTES**:
 * - [a]  a
 */

int tree_adaptation(void)
{
  //* >> Working in the refinement zones *//
  if (lmin < lmax)
  {

    clock_t aux_clock;

    struct node *ptr_node;
    int no_pts; // Number of parents in the loop
    int no_lvs; // Number of level of refinement to adapt

    no_lvs = GL_tentacles_level_max < (lmax - lmin) ? GL_tentacles_level_max : (GL_tentacles_level_max - 1);

    for (int lv = no_lvs; lv > -1; lv--)
    {
      GL_tentacles_size[lv + 1] = 0;
      no_pts = GL_tentacles_size[lv];
      //* >> "for" loop over parent nodes *//
      for (int i = 0; i < no_pts; i++)
      {
        ptr_node = GL_tentacles[lv][i];

        // printf("lv = %d, pt = %d, chn_size = %d\n", ptr_node->lv, i, ptr_node->chn_size);

        //* Updating the box mass information *//
        // printf("\n\nupdating_cell_struct\n\n");
        aux_clock = clock();
        if (updating_cell_struct(ptr_node) == _FAILURE_)
        {
          printf("Error at function updating_cell_struct()\n");
          return _FAILURE_;
        }
        GL_times[30] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

        aux_clock = clock();
        check_error(ptr_node, 1);
        GL_times[31] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;
        // printf("pre-check ends\n");

        //* Initialization of node boxes *//
        // printf("\n\ninitialization_node_boxes\n\n");
        aux_clock = clock();
        initialization_node_boxes(ptr_node);
        GL_times[32] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

        //* Initialization of the auiliary refinement arrays*//
        // printf("\n\nInitialization ref aux\n\n\n");
        // aux_clock = clock();
        // initialization_ref_aux(ptr_node);
        // GL_times[32] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

        //* >> Filling the refinement cells array *//
        // printf("\n\nFill cell ref\n\n\n");
        aux_clock = clock();
        if (fill_cell_ref(ptr_node) == _FAILURE_)
        {
          printf("Error at function fill_cell_ref()\n");
          return _FAILURE_;
        }
        GL_times[33] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

        //* >> Filling the different zones of refinement *//
        // printf("\n\nFill zones ref\n\n");
        aux_clock = clock();
        if (fill_zones_ref(ptr_node) == _FAILURE_)
        {
          printf("Error at function fill_zones_ref()\n");
          return _FAILURE_;
        }
        GL_times[34] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

        if (ptr_node->zones_size > 0)
        {
          //* >> Create links *//
          // printf("\n\nCreate links\n\n");
          aux_clock = clock();
          if (create_links(ptr_node) == _FAILURE_)
          {
            printf("Error at function create_links()\n");
            return _FAILURE_;
          }
          GL_times[35] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

          aux_clock = clock();
          if (ptr_node->chn_size > 0)
          {
            //* >> Removing cells that no longer require refinement *//
            // printf("\n\nRemoving cells nolonger requiere refinement\n\n");
            remov_cells_nolonger_require_refinement(ptr_node);
          }
          GL_times[36] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

          //* >> Adapting child boxes to the new space *//
          // printf("\n\nAdapt child box\n\n");
          aux_clock = clock();
          if (ptr_node->chn_size > 0)
          {
            if (adapt_child_nodes(ptr_node) == _FAILURE_)
            {
              printf("Error at function adapt_child_nodes()\n");
              return _FAILURE_;
            }
          }
          GL_times[37] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

          //* >> Creating and defining new children nodes for excess in refinement zones and and linking them to the parent node ptr_node *//
          // printf("\n\nCreate new child nodes\n\n");
          aux_clock = clock();
          if (ptr_node->zones_size > ptr_node->chn_size)
          {

            if (create_new_child_nodes(ptr_node) == _FAILURE_)
            {
              printf("Error, in new_child_nodes function\n");
              return _FAILURE_;
            }
          }
          GL_times[38] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

          //* >> Adapting the information from old child nodes to new child nodes *//
          // printf("\n\Transferring old child to new child\n\n");
          aux_clock = clock();
          if (ptr_node->chn_size > 0)
          {
            if (transferring_old_child_to_new_child(ptr_node) == _FAILURE_)
            {
              printf("Error at function transferring_old_child_to_new_child()\n");
              return _FAILURE_;
            }
          }
          GL_times[39] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

          //* >> Transferring new zones of refienemnt information to all child nodes *//
          // printf("\n\Transferring new zones to new child\n\n");
          aux_clock = clock();
          if (transferring_new_zones_to_new_child(ptr_node) == _FAILURE_)
          {
            printf("Error at function transferring_new_zones_to_new_child()\n");
            return _FAILURE_;
          }
          GL_times[40] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

          //* >> Reorganization child nodes *//
          // printf("\n\nReorganization child nodes\n\n");
          aux_clock = clock();
          reorganization_child_node(ptr_node);
          GL_times[42] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

          //* >> Reorganization grandchild nodes *//
          // printf("\n\nReorganization grandchild nodes\n\n");
          aux_clock = clock();
          if (ptr_node->chn_size > 0 && lv < no_lvs)
          {
            if (reorganization_grandchild_node(ptr_node) == _FAILURE_)
            {
              printf("Error at function reorganization_grandchild_node()\n");
              return _FAILURE_;
            }
          }
          GL_times[43] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

          //* >> Updating refinement zones of the grandchildren *//
          // printf("\n\nUpdating refinement zones of the grandchildren\n\n");
          aux_clock = clock();
          if (ptr_node->chn_size > 0 && lv < no_lvs)
          {
            updating_ref_zones_children(ptr_node);
          }
          GL_times[44] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

          //* >> Adding boundary simulation box information
          // printf("\n\nAdding boundary simulation box information\n\n");
          aux_clock = clock();
          adding_boundary_simulation_box_status_to_children_nodes(ptr_node);
          GL_times[41] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

          //* >> Updating children grid points *//
          // printf("\n\nUpdating children grid points\n\n");
          aux_clock = clock();
          if (filling_child_grid_point_arrays(ptr_node) == _FAILURE_)
          {
            printf("Error at function filling_child_grid_point_arrays()\n");
            return _FAILURE_;
          }
          GL_times[45] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

          // Filling number of particles to be updated outside of the refinement zones in the child nodes
          aux_clock = clock();
          for (int ch = 0; ch < ptr_node->zones_size; ch++)
          {
            computing_no_ptcl_outside_refinement_zones(ptr_node->pptr_chn[ch]);
          }
          GL_times[46] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;
        }

        //* >> Moved Unused child node to the stack of memory pool *//
        // printf("\n\nMoved Unused child node to the stack of memory pool\n\n");
        aux_clock = clock();
        transferring_unused_child_node_to_memory_pool(ptr_node);
        GL_times[47] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

        //* >> Tentacles updating *//
        // printf("\n\nTentacles Updating\n\n");
        aux_clock = clock();
        if (0 < ptr_node->zones_size)
        {
          if (tentacles_updating(ptr_node, lv + 1) == _FAILURE_)
          {
            printf("Error at function tentacles_updating()\n");
            return _FAILURE_;
          }
        }
        GL_times[48] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

        //* >> Updated new value for the chn size
        ptr_node->chn_size = ptr_node->zones_size;

        aux_clock = clock();
        check_error(ptr_node, 2);
        GL_times[31] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;
        // printf("post-check ends\n");
      }
    }

    // Filling number of particles to be updated outside of the refinement zones in the head node
    aux_clock = clock();
    computing_no_ptcl_outside_refinement_zones(GL_tentacles[0][0]);
    GL_times[46] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

    //* >> Tentacles Updating lv max *//
    // printf("\n\nTentacles updating lv max\n\n");
    aux_clock = clock();
    updating_tentacles_max_lv();
    GL_times[49] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;
  }

  //* >> Reseting

  return _SUCCESS_;
}