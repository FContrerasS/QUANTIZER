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
 * @file tree_construction.c
 *
 * \f[{\color{magenta} \mathbf{ DOCUMENTED\ ``tree\_construction.c"\ MODULE}}\f]
 *
 * @brief Construction of the tree to the initial particle configuration.
 *
 * **VERSION INFORMATION**: Felipe Contreras, 2022-10-01, version 1.0.
 *
 * **SHORT DESCRIPTION**: Construction of the tree of refinement zones to the
 * initial particle configuration. The nodes are created according to the
 * refinement zones. The structure of \ref tree_adaptation__REMINDER__Tentacles
 * "Tentacles" is modified to pointer to the initial configuration of nodes of
 * the tree.
 *
 * **PREREQUISITES**: Always called, but only used if the user uses more than
 * one level of refinement, i.e. \f$ l_{min} < l_{max} \f$.
 *
 * **RETURN**: The error status.
 *
 * **LONG DESCRIPTION**:
 *
 * The goal of the tree_construction.c module is to perform the construction of
 * the tree of refinements creating the nodes. In this module, the tentacles
 * structure is also modified to the initial node configuration. The
 * modifications go from the coarsest level of refinement to one level before
 * the finest level, running over every node of those levels. How this process
 * should be understood is by considering the parent node which is the current
 * node of the loop, and its child nodes. In every loop over a node, the parent
 * node localize the refinement zones and creates the new child nodes. At the
 * end of every loop the function modified the \ref
 * tree_adaptation__REMINDER__Tentacles "Tentacles" structure to be according to
 * the new child nodes.
 *
 * The flux of this function can be seen in the figure (work in progress), and
 * it is explained below:
 *
 * - [0]  <b> THE tree_construction.c MODULE STARTS....</b>
 *
 * - [1]  The tree_construction() begins to work.
 *
 * - [2]  Some useful parameters are defined.
 *
 * - [3]  Run a \c "while" loop over the refinement levels but avoiding the
 *   finest refinement level. The loop goes from a coarser refinement level to a
 *   finer one, starting from the coarsest level.
 *
 * - [4]  Run a \c "for" loop over all the nodes in the level of refinement
 *   begins. From now on we will call these nodes \f${\color{red} \mathbf{
 *   parent\
 *   nodes}}\f$.
 *
 * - [5]  With the initial configuration of the particles given by the user,
 *   cells in the parent node are labeled as refinement cells (-1 in the code
 *   language) in the box array \ref node.ptr_box "ptr_box" through the function
 *   fill_cell_ref().
 *
 * - [6]  Using the map of refinement cells obtained in the previous step, the
 *   different refinement zones (blocks in code language) in the parent node are
 *   built through the function fill_zones_ref(). Every one of these zones is
 *   independent of the other.
 *
 * - [7] Using the function create_child_nodes(), new child nodes are created to
 *   fit with the refinement zones of the parent node.
 *
 * - [8] Using the function transfer_tiny_child_nodes_to_memory_pool() the child
 *   nodes with fewer particles than the minimum particles required to refine a
 *   cell are removed from the parent node and put in the stack of the memory
 *   pool to be used in the future for any other parent node as its child node.
 *
 * - [9] Using the function fill_tentacles() the \ref
 *   tree_adaptation__REMINDER__Tentacles "Tentacles" structure is updated with
 *   the new child node distribution.
 *
 * - [10] <b> THE tree_construction.c MODULE ENDS....</b>
 *
 * **ILLUSTRATIVE EXAMPLES**:
 * - [a]  Trivial.
 *
 * **RATIONALES**:
 * - [a]  It is possible to join several functions in only one function, for
 *   example, the functions fill_cell_ref() and fill_zones_ref(), which are in
 *   charge of labeling the new box of the parent node with the cells to be
 *   refined, and create the refinement zones using this new map of refinement
 *   respectively. However, we decide to separate the functions as much as
 *   possible to increase the modularity and the cleanness of the code.
 *
 * **NOTES**:
 * - [a]  ll
 *
 */

#include "tree_construction.h"

//* >> Local Functions
static int find_min_max_subzones_ref_PERIODIC_BOUNDARY(struct node *ptr_node, int zone_idx);
static int fill_cell_ref(struct node *ptr_node);
static int fill_zones_ref(struct node *ptr_node);
static int create_child_nodes(struct node *ptr_node);
static void transfer_tiny_child_nodes_to_memory_pool(struct node *ptr_node);
static int fill_tentacles(const struct node *ptr_node);

/**
 * @brief Find the minimum and maximum cell indices of the subzone \e zone_idx.
 *
 * **SHORT DESCRIPTION**: Find the minimum and maximum cell indices in the *Code
 * Space* (see Key Concepts \ref Key_Concepts_Code_Space "Code Space") that
 * defines the *Smallest Box* (see Key Concepts \ref Key_Concepts_Smallest_Box
 * "Smallest Box") associated with the corresponding subzone \e zone_idx.
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
 * *Smallest Box* (see Key Concepts \ref Key_Concepts_Smallest_Box "Smallest
 * Box") ssociated with the corresponding subzone \e zone_idx., storing them in
 * the parent node structure parameters \ref node::ptr_pbc_min_subzones_x
 * "ptr_pbc_min_subzones_x" (\ref node::ptr_pbc_min_subzones_y "y", \ref
 * node::ptr_pbc_min_subzones_z "z"), and \ref node::ptr_pbc_max_subzones_x
 * "ptr_pbc_max_subzones_x" (\ref node::ptr_pbc_max_subzones_y "y", \ref
 * node::ptr_pbc_max_subzones_z "z").
 *
 * The flux of this function can be seen in the figure (work in progress), and
 * it is explained below:
 *
 * - [0]  <b> THE find_min_max_subzones_ref_PERIODIC_BOUNDARY() FUNCTION
 *   STARTS....</b>
 *
 * - [1]  Defining some internal useful parameters.
 *
 * - [2]  The current status of the new box contains the new map of refinement,
 *   which will correspond to the final one. In this step, to find the subzones,
 *   the cells of the zone of refinement with ID \e zone_idx are returned to the
 *   previous status to \e refinement_required (-1) in the new box.
 *
 * - [3]  How we are going to create the subzones is very similar to the
 *   creation of zones of refinement. However, there are some differences, for
 *   example, there is a difference in the criterion of the definition of zone
 *   and subzone of refinement. The first one accepts that a refinement zone
 *   crosses the box of the simulation, while the second one does not.
 *
 * - [4]  A \c "while" loop runs over all cells in the zone of refinement of ID
 *   \e zone_idx that have not been tagged to any subzone. A counter is used to
 *   perform this task.
 *
 * - [5]  In every step of the \c "while" loop over the cells, we ask if the
 *   cell belongs to a subzone of refinement (box value \f$ \geq 0 \f$). If it
 *   belongs, we continue asking to the next cell, if it does not belong (box
 *   value \f$ = -1 \f$), we change its box satus to "zone_idx" and pass to the
 *   next step.
 *
 * - [6]  Having found a cell with no subzone, we are going to create an entire
 *   new subzone starting with this cell as the foundation stone. To perform
 *   this, a new \c "while" loop is executed running until there are no more
 *   cells without analysis in the subzone, i.e. the final state of the block
 *   found is a solid piece completely isolated from the other subzones, labeled
 *   with the "subzone_idx" value in the parent node box.
 *
 * - [7]  At this point both \c "while" loops of the steps [4] and [6] end. Now,
 *   using the box with this new information the minimum and maximum of every
 *   subzone of the zone \e zone_idx are found.
 *
 * - [8]  Finally, the box status of the parent node is returned to its initial
 *   state, putting the value of \e zone_idx in the corresponding cell of the
 *   new box.
 *
 * - [9]  <b> THE find_min_max_subzones_ref_PERIODIC_BOUNDARY() FUNCTION
 *   ENDS....</b>
 *
 * **ILLUSTRATIVE EXAMPLES**:
 * - [a]  Trivial.
 *
 * **RATIONALES**:
 * - [a]  In step [2], we decide to use the same box to store the information
 *   about the new refinement zones. It requires overwriting the current
 *   information and then returning the old one. We choose this option instead
 *   of creating a copy of the box because we are avoiding to use malloc,
 *   memcpy, and free functions.
 *
 * **NOTES**:
 * - [a]  If the number of subzones found is equal to 1, then it is not
 *   necessary to find the minimum and maximum cell indices, because they  are
 *   equal to the minimum and maximum of the refinement level. It is  also valid
 *   when more than one dimension crosses the box of the simulation.
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
    // printf("Inside of the \c "while" loop\n");
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
      }                                              // End \c "while" loop, now the box contains the the information about all cells of the subzone "subzone_idx"
      subzone_idx++;                                 // Increasing the subzone number
    }                                                // subZone defined in the box
    cell_ref_idx++;
  } // At this point the box contains the information of all refinement subzones

  if (subzone_idx > 1)
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

  // ptr_node->pbc_subzones_cap = subzone_idx_max; // Maximum amount of subzones
  ptr_node->pbc_subzones_size = subzone_idx; // Total amount of subzones

  return _SUCCESS_;
}

/**
 * @brief Cells in the parent node are labeled as refinement cells in the \ref
 * node.ptr_box "ptr_box", and added to the auxiliary array \ref
 * node.ptr_cell_ref "ptr_cell_ref"
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
 * arrays (for example \ref node.ptr_cell_idx_x "ptr_cell_idx_x"). The criterion
 * to decide if a cell requires refinement is given by the user. If that
 * happens, the neighboring cells of the refined cell should be also refined
 * according to the \link n_exp \endlink parameter.
 *
 * The flux of this function can be seen in the figure (work in progress), and
 * it is explained below:
 *
 * - [0]  <b> THE fill_cell_ref() FUNCTION STARTS....</b>
 *
 * - [1]  Defining some internal useful parameters.
 *
 * - [2]  A \c "for" loop over the parent node cell is performed.
 *
 * - [3]  If the cell satisfies the refinement criteria, and its status hasn't
 *   been modified yet, it changes from "Exist" (-3) to "Requires-Refinement"
 *   (-1) in the box index
 *
 * - [4]  The neighboring cells at the distance of \link n_exp \endlink are also
 *   modified in their box indices if they haven't been modified yet.
 *
 * - [5]  At this point the \c "for" loop ends, and every cell in the node which
 *   will require refinement has been labeled to do that.
 *
 * - [6]  The next step is to add those cells to the auxiliary array of
 *   positional cell indices \ref node.ptr_cell_ref "ptr_cell_ref". To do that,
 *   a \c "for" loop running over the parent node cells is performed.
 *
 * - [7]  If the box index associated with the cell has the status of
 *   "Requires-Refinement" (-1), the positional index of the cell is added to
 *   the array \ref node.ptr_cell_ref "ptr_cell_ref".
 *
 * - [8]  The \c "for" loop over the parent node cells ends
 *
 * - [9]  <b> THE fill_cell_ref() FUNCTION ENDS....</b>
 *
 * **ILLUSTRATIVE EXAMPLES**:
 * - [a]  Trivial.
 *
 * **RATIONALES**:
 * - [a]  The reason we decided to add the elements to the \ref
 *   node.ptr_cell_ref "ptr_cell_ref" array, is because we want to avoid asking
 *   if there is enough capacity (see Key Concepts \ref Key_Concepts_Capacity
 *   "Capacity") in the array. in the array.
 *
 * **NOTES**:
 * - [a]
 */

static int fill_cell_ref(struct node *ptr_node)
{
  //* >> Adding cells which satisfy the refinement criteria to the array ptr_cell_ref and chaning the box the status of refinement -1 *//
  int box_idx_node;    // Box index
  int box_idxNbr_node; // Box index in the neigborhood

  int cell_ref_idx = 0; // Index of the position in the cell refined array

  int cntr; // Counter used to add cell index to the ptr_cell_ref

  int lv = ptr_node->lv;

  int box_real_dim_X_node = ptr_node->box_real_dim_x;
  int box_real_dim_X_times_Y_node = ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;

  //* >> Changing the box status from EXIST (-3) to REFINEMENT REQUIRED (-1) *//
  for (int i = 0; i < ptr_node->cell_size; i++)
  {
    box_idx_node = ptr_node->ptr_box_idx[i];

    // Refinement criterion in the box_mass in no border box points
    if (ptr_node->ptr_cell_struct[box_idx_node].cell_mass >= ref_criterion_mass || ptr_node->ptr_cell_struct[box_idx_node].ptcl_size >= ref_criterion_ptcl)
    //if (ptr_node->ptr_cell_struct[box_idx_node].cell_mass >= ref_criterion_mass && ptr_node->ptr_cell_struct[box_idx_node].ptcl_size >= ref_criterion_ptcl)
    {
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

  //* >> Adding the infomation about size of the ptr_cell_ref array *//
  ptr_node->cell_ref_size = cell_ref_idx;
  ptr_node->cell_ref_cap = cell_ref_idx;
  ptr_node->ptr_cell_ref = (int *)malloc(ptr_node->cell_ref_cap * sizeof(int));

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
 * refinement are separated in the refinement zones of the parent node.
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
 * groups which will be the refinement zones of the parent node. This zones are
 * created and stored in the pointer to the arrays \ref node.pptr_zones
 * "pptr_zones".
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
 * - [2]  A \c "while" loop runs over all cells which require refinement, i.e.
 *   cells belonging to the array \ref node.ptr_cell_ref "ptr_cell_ref", that
 *   have not been added to any refinement zone yet. A counter is used to
 *   perform this task.
 *
 * - [3]  In every step of the \c "while" loop over the cells, we ask if the
 *   cell belongs to a subzone of refinement (box value \f$ \geq 0 \f$). If it
 *   belongs, we continue asking to the next cell, if it does not belong (box
 *   value \f$ = -1 \f$), we change its box status to "zone_idx" and pass to the
 *   next step.
 *
 * - [4]  Having found a cell with no subzone, we are going to create an
 *   entirely new zone starting with this cell as the foundation stone. To
 *   perform this, a new \c "while" loop is executed running until there are no
 *   more cells without analysis in the zone, i.e. the final state of the block
 *   found is a solid piece completely isolated from the other zones, labeled
 *   with the "zone_idx" value in the parent node box.
 *
 * - [5]  At this point both \c "while" loops of the steps [2] and [3] end. Now,
 *   using the box with this information the zones array \ref node.pptr_zones
 *   "pptr_zones" is filled with the zones of refinement. To do this, a \c "for"
 *   loop over the refined cells array \ref node.ptr_cell_ref "ptr_cell_ref" is
 *   performed.
 *
 * - [6] <b> THE fill_zones_ref() FUNCTION ENDS....</b>
 *
 * **ILLUSTRATIVE EXAMPLES**:
 * - [a]  Trivial.
 *
 * **RATIONALES**:
 * - [a]  This design allows the contact of edges and corks between the
 *   different refinement zones, it also allows particles inside of a child node
 *   to jump to another child node of the same parent node in one time-step of
 *   the simulation. This leads to a computational spend (see
 *   particle_updating.c module) can be avoided  if we decided to change the
 *   design of the refinement zone considering now to ask to more cells than the
 *   closest 6 neighboring. However, this modification also implies the
 *   computational spend to ask to these new neighboring and also can increase
 *   the size of the refinement zone which in turn makes it difficult to load
 *   balance in multiple cores.
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

  int box_idx_node; // Box index

  int box_idxNbr_x_plus;  // Box index in the neigborhood on the right
  int box_idxNbr_x_minus; // Box index in the neigborhood on the left
  int box_idxNbr_y_plus;  // Box index in the neigborhood behind
  int box_idxNbr_y_minus; // Box index in the neigborhood in front
  int box_idxNbr_z_plus;  // Box index in the neigborhood up
  int box_idxNbr_z_minus; // Box index in the neigborhood down

  int lv = ptr_node->lv;

  int cell_idx;         // The cell index is simply i of the for loop
  int cell_ref_idx = 0; // The index in the cell refined array ptr_cell_ref

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
    box_idx_node = ptr_node->ptr_box_idx[cell_idx];

    if (ptr_node->ptr_box[box_idx_node] == -1) // A cell without zone has been founded
    {
      zone_size = 0; // Initial number of element in the zone

      //* >> Including the first element of the box to the auxiliary array ptr_aux_idx *//
      ptr_node->ptr_aux_idx[0] = box_idx_node;

      //* >>  Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0) *//
      ptr_node->ptr_box[box_idx_node] = zone_idx;

      zone_size++;               // +1 to the number of cells in the zone
      cntr_cell_add_all_zones++; // +1 to the number of cells added in total

      //* >> Building of the zone of ID = zone_idx *//
      cntr_insp = 0;                                                                     // Counter for the number of elements inspectioned in the current zone array
      while (zone_size > cntr_insp && ptr_node->cell_ref_size > cntr_cell_add_all_zones) // The cycle ends when all the elements of the zone have been inspected or when the number of cells refined is equal to the total number of cells added.
      {
        // Note that the number of elements in the zone increases if the neighbors of the inspected cell must be added to the zone

        cntr_cell_add = 0; // Counter number of cells added per cell inspected
        box_idx_node = ptr_node->ptr_aux_idx[cntr_insp];

        box_idxNbr_x_plus  = box_idx_node + 1;
        box_idxNbr_x_minus = box_idx_node - 1;
        box_idxNbr_y_plus  = box_idx_node + box_real_dim_X_node;
        box_idxNbr_y_minus = box_idx_node - box_real_dim_X_node;
        box_idxNbr_z_plus  = box_idx_node + box_real_dim_X_times_Y_node;
        box_idxNbr_z_minus = box_idx_node - box_real_dim_X_times_Y_node;

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
              if (ptr_node->ptr_box[box_idxNbr_z_plus] == -1 || ptr_node->ptr_box[box_idxNbr_z_plus] == zone_idx)
              {
                ptr_node->ptr_pbc_bool_bdry_anomalies_z[zone_idx] = true;
              }
            }
            else if (ptr_node->ptr_box[box_idxNbr_z_minus] == -6)
            {
              box_idxNbr_z_minus += (1 << lv) * box_real_dim_X_times_Y_node;
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
      }                                           // End \c "while" cycle, now the box contains the the information about all cells of the zone "zone_idx"

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
  } // At this point the box contains the information of all refinement zones

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
    box_idx_node = ptr_node->ptr_box_idx[cell_idx];
    zone_idx = ptr_node->ptr_box[box_idx_node];
    cntr_cell_add = ptr_node->ptr_aux_idx[zone_idx];          // Counter the element in the zone "zone_idx"
    ptr_node->pptr_zones[zone_idx][cntr_cell_add] = cell_idx; // Adding the index of the cell array in the block to the zone
    ptr_node->ptr_aux_idx[zone_idx] += 1;                     // Counter the number of elements added in the zone "zone_idx"
  }

  return _SUCCESS_;
} // End function fill_zones_ref

/**
 * @brief The child nodes are created to fit with their refinement zones.
 *
 * **SHORT DESCRIPTION**: The child nodes are created to fit with their
 * corresponding refinement zones.
 *
 * **PREREQUISITES**: If there is at least one refinement zone.
 *
 * @param[in,out] ptr_node Pointer to node structure
 *
 * **RETURN**: The error status.
 *
 * **LONG DESCRIPTION**:
 *
 * The child nodes are created to fit with their corresponding refinement zones.
 *
 * The flux of this function can be seen in the figure (work in progress), and
 * it is explained below:
 *
 * - [0]  <b> THE create_new_child_nodes() FUNCTION STARTS....</b>
 *
 * - [1]  Defining some internal useful parameters, and some parent node
 *   parameters are allocated.
 *
 * - [2]  A \c "for" loop over the refinement zones is performed.
 *
 * - [3]  The child node is obtained and initialized through the external
 *   function new_node(), and its level \ref node.lv "lv" and its identification
 *   \ref node.ID "ID" are assigned.
 *
 * - [4]  If periodic boundary conditions (pbc) is activated, and if any
 *   boundary flag anomalies \ref node.ptr_pbc_bool_bdry_anomalies_x
 *   "ptr_pbc_bool_bdry_anomalies_x" (\ref node.ptr_pbc_bool_bdry_anomalies_y
 *   "y", \ref node.ptr_pbc_bool_bdry_anomalies_z "z")  is activated, then the
 *   minimum and maximum of the subzones is computed through the function
 *   find_min_max_subzones_ref_PERIODIC_BOUNDARY().
 *
 * - [5]  For every dirección X, Y, and Z, the minimum and maximum, and the
 *   boundary flags are computed and stored in local parameters.
 *
 * - [6]  The minimum and maximum of the previous point were computed in the
 *   parent node, so then, they are transformed to the next level of refinement
 *   to be accord with the respective child node.
 *
 * - [7]  The box dimensions parameters and traslation constants (see \ref
 *   node.box_ts_x "box_ts_x" (\ref node.box_ts_y "y", \ref node.box_ts_z "z")
 *   are computed.
 *
 * - [8]  The cells in the refinement zone of the parent node are refined at
 *   packages of 8 cells, and added to the child node arrays \ref
 *   node.ptr_cell_idx_x "ptr_cell_idx_x" (\ref node.ptr_cell_idx_y
 *   "ptr_cell_idx_y", \ref node.ptr_cell_idx_z "ptr_cell_idx_z"), and \ref
 *   node.ptr_box_idx "ptr_box_idx".
 *
 * - [9]  The child box array is initialized at "No-Exist" status (-4), and
 *   "Exist" status (-3) for all child cells as appropriate.
 *
 * - [10] If the child node touches the boundary simulation in the non-periodic
 *   boundary condition (non-pbc), or crosses the whole simulation in the
 *   periodic boundary condition (pbc), then its box is updated in the outer
 *   boundary cells with the values of (-5) for non-pbc, and (-6) for pbc.
 *
 * - [11] The child cell structure is allocated, and filled with the particles
 *   and the corresponding parameters using the parent node cells information,
 *   and the external function ptcl_idx_to_box_idx() to know in which child cell
 *   the particle is located.
 *
 * - [12] Filling child node interior, boundary, and simulation boundary grid
 *   point arrays. Here are filled 4 arrays per type of grid point. For example,
 *   for the interior grid points of the child nodes, the function fills the
 *   arrays \ref node.ptr_intr_grid_cell_idx_x "ptr_intr_grid_cell_idx_x" (\ref
 *   node.ptr_intr_grid_cell_idx_y "ptr_intr_grid_cell_idx_y", \ref
 *   node.ptr_intr_grid_cell_idx_z "ptr_intr_grid_cell_idx_z",), and \ref
 *   node.ptr_intr_box_grid_idx "ptr_intr_box_grid_idx". Moreover, other grid
 *   parameters of the child node are updated, the size (see Key Concepts \ref
 *   Key_Concepts_Size "Size") and the capacity (see Key Concepts \ref
 *   Key_Concepts_Capacity "Capacity") of every type of grid point.
 *
 * - [13] The grid potentials \ref node.ptr_pot "ptr_pot", \ref node.ptr_pot_old
 *   "ptr_pot_old", accelerations \ref node.ptr_ax "ptr_ax", \ref node.ptr_ay
 *   "ptr_ay", \ref node.ptr_az "ptr_az", density \ref node.ptr_d "ptr_d", and
 *   the grid Capacity (see Key Concepts \ref Key_Concepts_Capacity "Capacity")
 *   parameter \ref node.grid_properties_cap "grid_properties_cap" are defined.
 *
 * - [14]  The child node is linked to the parent node through the struct node
 *   pointers parameters \ref node.pptr_chn "pptr_chn" and \ref node.ptr_pt
 *   "ptr_pt".
 *
 * - [15] Finally, the number of particles outside of the refinement zones is
 *   \ref node.no_ptcl_outs_ref_zones "no_ptcl_outs_ref_zones" is defined in the
 *   child node and updated in the parent node.
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

static int create_child_nodes(struct node *ptr_node)
{
  struct node *ptr_ch; // Pointer to the child node

  int cell_idx;  // Cell index
  int aux_idx_x; // Auxiliary indices
  int aux_idx_y;
  int aux_idx_z;

  int box_idx_ch;    // Box index
  int box_idxNbr_ch; // Box index of the neighboring child cell

  int box_idx_node; // Box index

  int box_idxNbr_i0_j0_k0_ch;    // Box index at x=0,y=0,z=0
  int box_idxNbr_im1_j0_k0_ch;   // Box index of the neighbor at x=-1,y=0,z=0
  int box_idxNbr_i0_jm1_k0_ch;   // Box index of the neighbor at x=0,y=-1,z=0
  int box_idxNbr_im1_jm1_k0_ch;  // Box index of the neighbor at x=-1,y=-1,z=0
  int box_idxNbr_i0_j0_km1_ch;   // Box index of the neighbor at x=0,y=0,z=-1
  int box_idxNbr_im1_j0_km1_ch;  // Box index of the neighbor at x=-1,y=0,z=-1
  int box_idxNbr_i0_jm1_km1_ch;  // Box index of the neighbor at x=0,y=-1,z=-1
  int box_idxNbr_im1_jm1_km1_ch; // Box index of the neighbor at x=-1,y=-1,z=-1

  int box_idx_ptcl_ch;

  int box_grid_idx_ch; // Box grid index

  bool is_bder_grid_point; // Ask if the grid point is interior
  bool grid_point_exist;
  bool grid_point_sim_bdry;

  int ptcl_idx;

  int pos_x; // Distance between the real box and the "Smallest Box" when the last is localized in the middle of the real one
  int pos_y;
  int pos_z;

  int cap;  // Capacity of arrays in child nodes
  int size; // Size or number of elements in some array in child nodes

  int lv = ptr_node->lv;

  int box_real_dim_X_ch;
  int box_real_dim_X_times_Y_ch;

  int grid_box_real_dim_X_ch;
  int grid_box_real_dim_X_times_Y_ch;

  int size_bder_grid_points_ch;
  int size_intr_grid_points_ch;

  bool check; // check subzone analized

  int aux_subzones_analized;

  int box_idx_x_ch;
  int box_idx_y_ch;
  int box_idx_z_ch;

  //* >> Creating child nodes in parent node *//
  ptr_node->chn_cap = ptr_node->zones_cap; // Same amount than refinement zones
  ptr_node->chn_size = ptr_node->zones_size;

  // if (ptr_node->chn_cap > 0)
  // {
  //   ptr_node->pptr_chn = (struct node **)malloc(ptr_node->chn_cap * sizeof(struct node *)); // Allocating children pointers
  // }
  ptr_node->pptr_chn = (struct node **)malloc(ptr_node->chn_cap * sizeof(struct node *)); // Allocating children pointers
  for (int i = 0; i < ptr_node->chn_cap; i++)
  {
    ptr_node->pptr_chn[i] = NULL;
  }

  //* >> Filling child nodes *//
  for (int zone_idx = 0; zone_idx < ptr_node->zones_size; zone_idx++)
  {
    //* >> Create and Initialize child node *//
    ptr_ch = (struct node *)malloc(sizeof(struct node));
    initialize_node(ptr_ch);

    //* >> Global properties *//
    ptr_ch->ID = zone_idx;
    ptr_ch->lv = lv + 1;

    //* >> Boxes *//
    // MIN and MAX cell indices values of the node.
    if (bdry_cond_type == 0 && ptr_node->pbc_crosses_whole_sim_box == true)
    {
      if (ptr_node->ptr_pbc_bool_bdry_anomalies_x[zone_idx] == true ||
          ptr_node->ptr_pbc_bool_bdry_anomalies_y[zone_idx] == true ||
          ptr_node->ptr_pbc_bool_bdry_anomalies_z[zone_idx] == true)
      {
        if (find_min_max_subzones_ref_PERIODIC_BOUNDARY(ptr_node, zone_idx) == _FAILURE_)
        {
          printf("Error at function fill_subzones_ref_PERIODIC_BOUNDARY()\n");
          return _FAILURE_;
        }
      }
    }

    //* >> X axis
    if (ptr_node->pbc_crosses_whole_sim_box_x == true)
    {
      if (ptr_node->ptr_pbc_bool_bdry_anomalies_x[zone_idx] == false)
      {
        for (int j = 0; j < ptr_node->ptr_zone_size[zone_idx]; j++)
        {
          cell_idx = ptr_node->pptr_zones[zone_idx][j]; // Cell index in the parent zone
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
      if (ptr_ch->box_max_x - ptr_ch->box_min_x >= (1 << lv) - 1)
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
      for (int j = 0; j < ptr_node->ptr_zone_size[zone_idx]; j++)
      {
        cell_idx = ptr_node->pptr_zones[zone_idx][j]; // Cell index in the parent zone

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
    } // End else if(ptr_node->pbc_crosses_sim_box_bdry_x = true)
    // Case Parent x axis do not cross the simulation box
    else
    {
      for (int j = 0; j < ptr_node->ptr_zone_size[zone_idx]; j++)
      {
        cell_idx = ptr_node->pptr_zones[zone_idx][j]; // Cell index in the parent zone
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
      if (ptr_node->ptr_pbc_bool_bdry_anomalies_y[zone_idx] == false)
      {
        for (int j = 0; j < ptr_node->ptr_zone_size[zone_idx]; j++)
        {
          cell_idx = ptr_node->pptr_zones[zone_idx][j]; // Cell index in the parent zone
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
      if (ptr_ch->box_max_y - ptr_ch->box_min_y >= (1 << lv) - 1)
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
      for (int j = 0; j < ptr_node->ptr_zone_size[zone_idx]; j++)
      {
        cell_idx = ptr_node->pptr_zones[zone_idx][j]; // Cell index in the parent zone

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
    } // End else if(ptr_node->pbc_crosses_sim_box_bdry_y = true)
    // Case Parent x axis do not cross the simulation box
    else
    {
      for (int j = 0; j < ptr_node->ptr_zone_size[zone_idx]; j++)
      {
        cell_idx = ptr_node->pptr_zones[zone_idx][j]; // Cell index in the parent zone
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
      if (ptr_node->ptr_pbc_bool_bdry_anomalies_z[zone_idx] == false)
      {
        for (int j = 0; j < ptr_node->ptr_zone_size[zone_idx]; j++)
        {
          cell_idx = ptr_node->pptr_zones[zone_idx][j]; // Cell index in the parent zone
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
      if (ptr_ch->box_max_z - ptr_ch->box_min_z >= (1 << lv) - 1)
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
      for (int j = 0; j < ptr_node->ptr_zone_size[zone_idx]; j++)
      {
        cell_idx = ptr_node->pptr_zones[zone_idx][j]; // Cell index in the parent zone

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
    } // End else if(ptr_node->pbc_crosses_sim_box_bdry_z = true)
    // Case Parent x axis do not cross the simulation box
    else
    {
      for (int j = 0; j < ptr_node->ptr_zone_size[zone_idx]; j++)
      {
        cell_idx = ptr_node->pptr_zones[zone_idx][j]; // Cell index in the parent zone
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

    // Changing the min and max of the ""Smallest Box"" from parent units to child units
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
    // ptr_ch->box_real_dim_x = 5 > (n_exp - 1) ? (ptr_ch->box_dim_x + 10) : (ptr_ch->box_dim_x + 2 * n_exp - 2);
    // ptr_ch->box_real_dim_y = 5 > (n_exp - 1) ? (ptr_ch->box_dim_y + 10) : (ptr_ch->box_dim_y + 2 * n_exp - 2);
    // ptr_ch->box_real_dim_z = 5 > (n_exp - 1) ? (ptr_ch->box_dim_z + 10) : (ptr_ch->box_dim_z + 2 * n_exp - 2);

    ptr_ch->box_real_dim_x = ptr_ch->box_dim_x + min_box_extra_size_per_side_real;
    ptr_ch->box_real_dim_y = ptr_ch->box_dim_y + min_box_extra_size_per_side_real;
    ptr_ch->box_real_dim_z = ptr_ch->box_dim_z + min_box_extra_size_per_side_real;

    // ptr_ch->box_real_dim_x_old = ptr_ch->box_real_dim_x;
    // ptr_ch->box_real_dim_y_old = ptr_ch->box_real_dim_y;
    // ptr_ch->box_real_dim_z_old = ptr_ch->box_real_dim_z;

    // Translations between cell array and box
    pos_x = (ptr_ch->box_real_dim_x - ptr_ch->box_dim_x) / 2; // Half of the distance of the box side less the ""Smallest Box"" side
    ptr_ch->box_ts_x = ptr_ch->box_min_x - pos_x;             // Every cell in the level l in the box must be subtracted this value to obtain the box index
    pos_y = (ptr_ch->box_real_dim_y - ptr_ch->box_dim_y) / 2;
    ptr_ch->box_ts_y = ptr_ch->box_min_y - pos_y;
    pos_z = (ptr_ch->box_real_dim_z - ptr_ch->box_dim_z) / 2;
    ptr_ch->box_ts_z = ptr_ch->box_min_z - pos_z;

    // ptr_ch->box_ts_x_old = ptr_ch->box_ts_x;
    // ptr_ch->box_ts_y_old = ptr_ch->box_ts_y;
    // ptr_ch->box_ts_z_old = ptr_ch->box_ts_z;

    //* >> Cells in the node *//
    cap = 2 * 8 * ptr_node->ptr_zone_size[zone_idx]; // Capacity is the double of total child cells
    ptr_ch->cell_cap = cap;
    size = 8 * ptr_node->ptr_zone_size[zone_idx];
    ptr_ch->cell_size = size;
    ptr_ch->ptr_cell_idx_x = (int *)malloc(ptr_ch->cell_cap * sizeof(int));
    ptr_ch->ptr_cell_idx_y = (int *)malloc(ptr_ch->cell_cap * sizeof(int));
    ptr_ch->ptr_cell_idx_z = (int *)malloc(ptr_ch->cell_cap * sizeof(int));
    ptr_ch->ptr_box_idx = (int *)malloc(ptr_ch->cell_cap * sizeof(int));

    box_real_dim_X_ch = ptr_ch->box_real_dim_x;
    box_real_dim_X_times_Y_ch = ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;

    for (int j = 0; j < ptr_node->ptr_zone_size[zone_idx]; j++)
    {
      aux_idx_x = ptr_node->ptr_cell_idx_x[ptr_node->pptr_zones[zone_idx][j]] * 2;
      aux_idx_y = ptr_node->ptr_cell_idx_y[ptr_node->pptr_zones[zone_idx][j]] * 2;
      aux_idx_z = ptr_node->ptr_cell_idx_z[ptr_node->pptr_zones[zone_idx][j]] * 2;

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
            cell_idx = ii + jj * 2 + kk * 4 + j * 8; // Child cell index
            ptr_ch->ptr_cell_idx_x[cell_idx] = aux_idx_x + ii;
            ptr_ch->ptr_cell_idx_y[cell_idx] = aux_idx_y + jj;
            ptr_ch->ptr_cell_idx_z[cell_idx] = aux_idx_z + kk;
            ptr_ch->ptr_box_idx[cell_idx] = box_idx_ch + ii + jj * box_real_dim_X_ch + kk * box_real_dim_X_times_Y_ch;
          }
        }
      }
    }

    // Filling the box status
    cap = ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y * ptr_ch->box_real_dim_z; // In general, the size of each side must be 5 times bigger than the same side of the ""Smallest Box""
    ptr_ch->box_cap = cap;
    ptr_ch->ptr_box = (int *)malloc(cap * sizeof(int));
    ptr_ch->ptr_box_old = (int *)malloc(cap * sizeof(int));
    // Putting the value of NO-EXIST (-4) in every box index
    for (int j = 0; j < cap; j++)
    {
      ptr_ch->ptr_box[j] = -4;
    }
    // Changing from NO-EXIST (-4) to EXIST (-3) to all cells in the block of the child node
    for (int j = 0; j < ptr_ch->cell_size; j++)
    {
      box_idx_ch = ptr_ch->ptr_box_idx[j];
      ptr_ch->ptr_box[box_idx_ch] = -3;
    }

    // Adding -6 to the boundary when it corresponds
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
    // Adding -5 to the boundary when it corresponds
    else if (bdry_cond_type != 0 && ptr_ch->sim_bdry_contact == true)
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

    //* >> Struct of cells (Particles and cell mass)
    //* >> Total mass in the node *//

    ptr_ch->ptr_cell_struct = (struct cell_struct *)malloc(ptr_ch->box_cap * sizeof(struct cell_struct));
    // ptr_ch->ptr_cell_struct_old = (struct cell_struct *)malloc(ptr_ch->box_cap * sizeof(struct cell_struct));
    for (int j = 0; j < ptr_ch->box_cap; j++)
    {
      initialize_cell_struct(&(ptr_ch->ptr_cell_struct[j]));
      // initialize_cell_struct(&(ptr_ch->ptr_cell_struct_old[j]));
    }

    for (int j = 0; j < ptr_node->ptr_zone_size[zone_idx]; j++)
    {
      box_idx_node = ptr_node->ptr_box_idx[ptr_node->pptr_zones[zone_idx][j]];
      // box_idx_x_node = ptr_node->ptr_cell_idx_x[ptr_node->pptr_zones[i][j]] - ptr_node->box_ts_x;
      // box_idx_y_node = ptr_node->ptr_cell_idx_y[ptr_node->pptr_zones[i][j]] - ptr_node->box_ts_y;
      // box_idx_z_node = ptr_node->ptr_cell_idx_z[ptr_node->pptr_zones[i][j]] - ptr_node->box_ts_z;
      // box_idx_node = box_idx_x_node + box_idx_y_node * ptr_node->box_real_dim_x + box_idx_z_node * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;
      if (ptr_node->ptr_cell_struct[box_idx_node].ptcl_size > 0)
      {
        // Allocating space fot the child particles in the 8  cells
        aux_idx_x = ptr_node->ptr_cell_idx_x[ptr_node->pptr_zones[zone_idx][j]] * 2;
        aux_idx_y = ptr_node->ptr_cell_idx_y[ptr_node->pptr_zones[zone_idx][j]] * 2;
        aux_idx_z = ptr_node->ptr_cell_idx_z[ptr_node->pptr_zones[zone_idx][j]] * 2;

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

        cap = ptr_node->ptr_cell_struct[box_idx_node].ptcl_size;
        for (int kk = 0; kk < 2; kk++)
        {
          for (int jj = 0; jj < 2; jj++)
          {
            for (int ii = 0; ii < 2; ii++)
            {
              box_idxNbr_ch = box_idx_ch + ii + jj * box_real_dim_X_ch + kk * box_real_dim_X_times_Y_ch;
              ptr_ch->ptr_cell_struct[box_idxNbr_ch].ptcl_cap = cap;
              ptr_ch->ptr_cell_struct[box_idxNbr_ch].ptr_ptcl = (int *)malloc(cap * sizeof(int));
            }
          }
        }

        // Putting parent cell particles into the child cell
        for (int k = 0; k < ptr_node->ptr_cell_struct[box_idx_node].ptcl_size; k++)
        {
          ptcl_idx = ptr_node->ptr_cell_struct[box_idx_node].ptr_ptcl[k];

          box_idx_ptcl_ch = ptcl_idx_to_box_idx(ptr_ch, ptcl_idx);

          ptr_ch->ptr_cell_struct[box_idx_ptcl_ch].ptr_ptcl[ptr_ch->ptr_cell_struct[box_idx_ptcl_ch].ptcl_size] = ptcl_idx;
          ptr_ch->ptr_cell_struct[box_idx_ptcl_ch].ptcl_size += 1;
          ptr_ch->ptr_cell_struct[box_idx_ptcl_ch].cell_mass += GL_ptcl_mass[ptcl_idx];
        }

        ptr_ch->node_mass += ptr_node->ptr_cell_struct[box_idx_node].cell_mass;              // Local mass
        ptr_ch->no_ptcl_full_node += ptr_node->ptr_cell_struct[box_idx_node].ptcl_size; // Local number of particles
      }
    }

    // size = ptr_ch->cell_size / 8 * 18 + 9; // 27 * N - 9 * (N-1)

    // 27 * N - 9 * (N-1) = 18N+9 >= Total of grid points; N = number of cubics with 8 cells and 27 grid points
    // Maximum of border grid points = 9*(N*2+1) - (N*2-1) = 16N+10 (Imagine a row of cubics)
    // Maximum of interior grid points = (2k-1)^3 < 8k^3, where k^3 = N
    size_bder_grid_points_ch = ptr_ch->cell_size / 8 * 16 + 10;
    size_intr_grid_points_ch = ptr_ch->cell_size;

    //* >> Space checking of border grid points array*//
    ptr_ch->grid_intr_cap = size_intr_grid_points_ch;
    ptr_ch->grid_bdry_cap = size_bder_grid_points_ch;
    ptr_ch->ptr_intr_grid_cell_idx_x = (int *)malloc(ptr_ch->grid_intr_cap * sizeof(int));
    ptr_ch->ptr_intr_grid_cell_idx_y = (int *)malloc(ptr_ch->grid_intr_cap * sizeof(int));
    ptr_ch->ptr_intr_grid_cell_idx_z = (int *)malloc(ptr_ch->grid_intr_cap * sizeof(int));
    ptr_ch->ptr_intr_box_grid_idx = (int *)malloc(ptr_ch->grid_intr_cap * sizeof(int));
    ptr_ch->ptr_bdry_grid_cell_idx_x = (int *)malloc(ptr_ch->grid_bdry_cap * sizeof(int));
    ptr_ch->ptr_bdry_grid_cell_idx_y = (int *)malloc(ptr_ch->grid_bdry_cap * sizeof(int));
    ptr_ch->ptr_bdry_grid_cell_idx_z = (int *)malloc(ptr_ch->grid_bdry_cap * sizeof(int));
    ptr_ch->ptr_bdry_box_grid_idx = (int *)malloc(ptr_ch->grid_bdry_cap * sizeof(int));

    if (ptr_ch->sim_bdry_contact == true)
    {
      ptr_ch->grid_sim_bdry_cap = size_bder_grid_points_ch;

      ptr_ch->ptr_sim_bdry_grid_cell_idx_x = (int *)malloc(ptr_ch->grid_sim_bdry_cap * sizeof(int));
      ptr_ch->ptr_sim_bdry_grid_cell_idx_y = (int *)malloc(ptr_ch->grid_sim_bdry_cap * sizeof(int));
      ptr_ch->ptr_sim_bdry_grid_cell_idx_z = (int *)malloc(ptr_ch->grid_sim_bdry_cap * sizeof(int));
      ptr_ch->ptr_sim_bdry_box_grid_idx = (int *)malloc(ptr_ch->grid_sim_bdry_cap * sizeof(int));
    }

    //* >> Grid points *//
    grid_box_real_dim_X_ch = ptr_ch->box_real_dim_x + 1;
    grid_box_real_dim_X_times_Y_ch = (ptr_ch->box_real_dim_x + 1) * (ptr_ch->box_real_dim_y + 1);

    for (int kk = ptr_ch->box_min_z - ptr_ch->box_ts_z; kk < ptr_ch->box_max_z - ptr_ch->box_ts_z + 2; kk++)
    {
      for (int jj = ptr_ch->box_min_y - ptr_ch->box_ts_y; jj < ptr_ch->box_max_y - ptr_ch->box_ts_y + 2; jj++)
      {
        for (int ii = ptr_ch->box_min_x - ptr_ch->box_ts_x; ii < ptr_ch->box_max_x - ptr_ch->box_ts_x + 2; ii++)
        {
          box_idx_ch = ii + jj * ptr_ch->box_real_dim_x + kk * ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;

          grid_point_exist = false;
          grid_point_sim_bdry = false;

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
                grid_point_sim_bdry = true;
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
            if (grid_point_sim_bdry == true)
            {
              ptr_ch->ptr_sim_bdry_grid_cell_idx_x[ptr_ch->grid_sim_bdry_size] = aux_idx_x;
              ptr_ch->ptr_sim_bdry_grid_cell_idx_y[ptr_ch->grid_sim_bdry_size] = aux_idx_y;
              ptr_ch->ptr_sim_bdry_grid_cell_idx_z[ptr_ch->grid_sim_bdry_size] = aux_idx_z;
              ptr_ch->ptr_sim_bdry_box_grid_idx[ptr_ch->grid_sim_bdry_size] = box_grid_idx_ch;
              ptr_ch->grid_sim_bdry_size += 1; // Increasing the number of border grid points in the array
            }
            //* >> Border grid point*//
            else if (is_bder_grid_point == true)
            {
              //* >> Adding the grid point to the border array *//
              ptr_ch->ptr_bdry_grid_cell_idx_x[ptr_ch->grid_bdry_size] = aux_idx_x;
              ptr_ch->ptr_bdry_grid_cell_idx_y[ptr_ch->grid_bdry_size] = aux_idx_y;
              ptr_ch->ptr_bdry_grid_cell_idx_z[ptr_ch->grid_bdry_size] = aux_idx_z;
              ptr_ch->ptr_bdry_box_grid_idx[ptr_ch->grid_bdry_size] = box_grid_idx_ch;
              ptr_ch->grid_bdry_size += 1; // Increasing the number of border grid points in the array
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

    //* Potential, Acceleration and density of the grid *//
    cap = (ptr_ch->box_real_dim_x + 1) * (ptr_ch->box_real_dim_y + 1) * (ptr_ch->box_real_dim_z + 1);
    ptr_ch->grid_properties_cap = cap;
    ptr_ch->ptr_pot = (vtype *)calloc(cap, sizeof(vtype));
    ptr_ch->ptr_pot_old = (vtype *)calloc(cap, sizeof(vtype));
    ptr_ch->ptr_ax = (vtype *)calloc(cap, sizeof(vtype));
    ptr_ch->ptr_ay = (vtype *)calloc(cap, sizeof(vtype));
    ptr_ch->ptr_az = (vtype *)calloc(cap, sizeof(vtype));
    ptr_ch->ptr_d = (vtype *)calloc(cap, sizeof(vtype));

    //* >> Tree structure *//
    ptr_node->pptr_chn[zone_idx] = ptr_ch; // Parent node pointing to child node i
    ptr_ch->ptr_pt = ptr_node;             // Child node i pointig to parent node

    // Filling number of particles to be updated outside of the refinement zones
    ptr_ch->no_ptcl_outs_ref_zones = ptr_ch->no_ptcl_full_node;
    ptr_node->no_ptcl_outs_ref_zones -= ptr_ch->no_ptcl_full_node;

  } // End filling child nodes

  return _SUCCESS_;
} // end function create_child_nodes

/**
 * @brief Transferring child nodes with a lower number of particles to the stack
 * of the memory pool
 *
 * **SHORT DESCRIPTION**: Transferring child nodes with fewer particles than the
 * minimum particles required to refine a cell to the stack of the memory pool
 * to be used in the future for any other parent node as its child node.
 *
 * **PREREQUISITES**: If there is at least one refinement zone.
 *
 * @param[in,out] ptr_node Pointer to node structure
 *
 * **RETURN**: No parameter is returned.
 *
 * **LONG DESCRIPTION**:
 *
 * Remove refinement zones of the parent node which contain fewer particles than
 * the minimum particles required to refine a cell. So, this zone is removed,
 * and the corresponding child node is also removed and transferred to the stack
 * of the memory pool. The gap in the parent node is filled by the last node in
 * the list of child nodes, and also this surrogate child node is updated in its
 * parameters.
 *
 * The flux of this function can be seen in the figure (work in progress), and
 * it is explained below:
 *
 * - [0]  <b> THE transfer_tiny_child_nodes_to_memory_pool() FUNCTION
 *   STARTS....</b>
 *
 * - [1]  Defining some internal useful parameters.
 *
 * - [4]  A \c "for" loop over the child nodes is performed.
 *
 * - [5]  For every child node, it has fewer particles than the refinement cell
 *   criteria, it is replaced by the last node in the list.
 *
 * - [6]  The surrogating node is updated in its identification, and the box
 *   cells of the parent node are updated with the information of the
 *   surrogating node and the removed node.
 *
 * - [7]  The removed node is sent to the stack of the memory pool.
 *
 * - [8]  Finally, some other important parameters are updated.
 *
 * - [9]  <b> THE transfer_tiny_child_nodes_to_memory_pool() FUNCTION
 *   ENDS....</b>
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

static void transfer_tiny_child_nodes_to_memory_pool(struct node *ptr_node)
{
  int box_idx_node;

  int aux_int;
  int *aux_ptr_int;

  //* >> Removing zones (children nodes) with 0 particles
  for (int ch = 0; ch < ptr_node->chn_size; ch++)
  {
    if (ptr_node->pptr_chn[ch]->no_ptcl_full_node < ref_criterion_ptcl)
    {
      //* >> Changing the ID of the replacement zone
      ptr_node->pptr_chn[ptr_node->chn_size - 1]->ID = ch;

      //* >> Changing the box status of the replacement zone from the old zone value "ptr_node->chn_size - 1" to the new one = ch
      for (int cell_idx = 0; cell_idx < ptr_node->ptr_zone_size[ptr_node->chn_size - 1]; cell_idx++)
      {
        box_idx_node = ptr_node->ptr_box_idx[ptr_node->pptr_zones[ptr_node->chn_size - 1][cell_idx]];
        ptr_node->ptr_box[box_idx_node] = ch;
      }

      //* >> Changing the box status of the removed zone from the old zone value ch to the new one = -3
      for (int cell_idx = 0; cell_idx < ptr_node->ptr_zone_size[ch]; cell_idx++)
      {
        box_idx_node = ptr_node->ptr_box_idx[ptr_node->pptr_zones[ch][cell_idx]];
        ptr_node->ptr_box[box_idx_node] = -3;
      }

      //* >> Removing the zone ch without particles from the children of the node ptr_node
      add_node_to_stack(ptr_node->pptr_chn[ch]);
      
      //* >> Replacing the ch zone for the last children
      ptr_node->pptr_chn[ch] = ptr_node->pptr_chn[ptr_node->chn_size - 1];
      ptr_node->pptr_chn[ptr_node->zones_size - 1] = NULL;

      //* >> Updating the zone values
      // Size of the zones
      aux_int = ptr_node->ptr_zone_size[ch];
      ptr_node->ptr_zone_size[ch] = ptr_node->ptr_zone_size[ptr_node->zones_size - 1];
      ptr_node->ptr_zone_size[ptr_node->zones_size - 1] = aux_int;

      // Cap of the zones
      aux_int = ptr_node->ptr_zone_cap[ch];
      ptr_node->ptr_zone_cap[ch] = ptr_node->ptr_zone_cap[ptr_node->zones_size - 1];
      ptr_node->ptr_zone_cap[ptr_node->zones_size - 1] = aux_int;

      // Zone array
      aux_ptr_int = ptr_node->pptr_zones[ch];
      ptr_node->pptr_zones[ch] = ptr_node->pptr_zones[ptr_node->zones_size - 1];
      ptr_node->pptr_zones[ptr_node->zones_size - 1] = aux_ptr_int;

      // Total number of active zones
      ptr_node->chn_size -= 1; //-1 the the total number of children in the node
                                // ptr_node->zones_size -= 1;	//-1 the the total number of zones in the node

      ch--;
    }
  }
}

/**
 * @brief Filling of the \ref tree_adaptation__REMINDER__Tentacles "Tentacles"
 * at the child nodes level of refinement
 *
 * **SHORT DESCRIPTION**: Filling of the \ref
 * tree_adaptation__REMINDER__Tentacles "Tentacles" at the child nodes level of
 * refinement
 *
 * **PREREQUISITES**: If there is at least one refinement zone.
 *
 * @param[in] ptr_node Pointer to node structure
 *
 * **RETURN**: The error status.
 *
 * **LONG DESCRIPTION**:
 *
 * Filling of the \ref tree_adaptation__REMINDER__Tentacles "Tentacles" at the
 * child nodes level of refinement. Moreover, the maximum level parameter \link
 * GL_tentacles_level_max \endlink of the \ref
 * tree_adaptation__REMINDER__Tentacles "Tentacles" is updated in every call.
 *
 * Here the only global parameters associated with the \ref
 * tree_adaptation__REMINDER__Tentacles "Tentacles" are going to be updated.
 *
 * The flux of this function can be seen in the figure (work in progress), and
 * it is explained below:
 *
 * - [0]  <b> THE fill_tentacles() FUNCTION STARTS....</b>
 *
 * - [1]  Defining some internal useful parameters.
 *
 * - [2]  Run a \c "for" loop over the number of refinement zones, that at this
 *   point of the module tree_adaptation.c is equal to the final number of child
 *   nodes.
 *
 * - [3]  The child node is pointed by the tentacle at the corresponding
 *   position and level of refinement.
 *
 * - [4]  Finally, the tentacles Size (see Key Concepts \ref Key_Concepts_Size
 *   "Size") of the global parameter \link GL_tentacles_size \endlink, and the
 *   \link GL_tentacles_level_max \endlink are updated.
 *
 * - [5]  <b> THE fill_tentacles() FUNCTION ENDS....</b>
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

static int fill_tentacles(const struct node *ptr_node)
{
  int tentacle_lv = ptr_node->lv - lmin + 1; // Children level
  int no_tentacles = GL_tentacles_size[tentacle_lv];
  int no_chn = ptr_node->chn_size;
  int size = no_tentacles + no_chn;
  

  if (space_check(&(GL_tentacles_cap[tentacle_lv]), size, 4.0f, "p1n2", &(GL_tentacles[tentacle_lv])) == _FAILURE_)
  {
    printf("Error, in space_check function\n");
    return _FAILURE_;
  }

  for (int i = 0; i < no_chn; i++)
  {
    //* >> Putting elements in the new tentacles *//
    GL_tentacles[tentacle_lv][no_tentacles + i] = ptr_node->pptr_chn[i];
  }

  //* Increasing the number of structs in the level lv *//
  GL_tentacles_size[tentacle_lv] = size;

  //* >> Increasing the deepth of levels of the tentacles *//
  if (GL_tentacles_level_max < tentacle_lv)
  {
    GL_tentacles_level_max++;
  }

  return _SUCCESS_;
} // end function fill_tentacles

/**
 * @brief Make the calls of all the local functions of the tree_construction.c
 * module
 *
 * **SHORT DESCRIPTION**: Make the calls of all the local functions of the
 * tree_construction.c module
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
 * - [a]  Example of a Image:
 *
 * 
 *
 */

int tree_construction(void)
{
  struct node *ptr_node; // Parent node

  int no_pts; // Number of parents in the cycle

  int lv = 0;

  while (lv < lmax - lmin && lv <= GL_tentacles_level_max)
  {
    no_pts = GL_tentacles_size[lv];
    //* >> For cycle over parent nodes *//
    for (int i = 0; i < no_pts; i++)
    {
      ptr_node = GL_tentacles[lv][i];

      //* >> Filling the refinement cells array *//
      if (fill_cell_ref(ptr_node) == _FAILURE_)
      {
        printf("Error at function fill_cell_ref()\n");
        return _FAILURE_;
      }

      //* >> Filling the different zones of refinement *//
      if (fill_zones_ref(ptr_node) == _FAILURE_)
      {
        printf("Error at function fill_zones_ref()\n");
        return _FAILURE_;
      }

      if (ptr_node->zones_size > 0)
      {
        //* >> Creating child nodes *//
        if (create_child_nodes(ptr_node) == _FAILURE_)
        {
          printf("Error at function create_child_nodes()\n");
          return _FAILURE_;
        }

        //* >> Removing child nodes with 0 particles
        transfer_tiny_child_nodes_to_memory_pool(ptr_node);

        //* >> Filling Tentacles for the next cycle at level of refinement
        if (fill_tentacles(ptr_node) == _FAILURE_)
        {
          printf("Error at function fill_tentacles()\n");
          return _FAILURE_;
        }
      }

    } // End of cycle over parent nodes
    lv++;

  } // End of cycle over refinement levels

  return _SUCCESS_;
}