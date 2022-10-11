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
 * @file initialize_node.h ******************** Documented \e "initialize_node.h" header ******************** \n
 *
 * @brief This is the header file of the initialize_node.c script.
 *
 * **VERSION INFORMATION**: Felipe Contreras, 2022-10-01, version 1.0.
 *
 * **DESCRIPTION**: This is the header file of the initialize_node.c script.
 *
 * **PREREQUISITES**: Always used.
 */

#ifndef __INITIALIZENODE__
#define __INITIALIZENODE__

#include "common.h"

/**
 * @brief Stand-alone unit representing a refinement zone of its parent node.
 *
 * **SHORT DESCRIPTION**: Stand-alone unit to other nodes in the same
 * level of refinement. Nodes represent a refinement zone of their parent node.
 *
 * **Key Concepts**:
 * - \anchor Key_Concepts_Capacity [a] **Capacity**:  Maximum amount of elements 
 * that can be held in an array.
 *
 * - \anchor Key_Concepts_Size [b] **Size**: Current number of elements in an 
 * array.
 *
 * - \anchor Key_Concepts_Smallest_Box [c] **Smallest Box**:  Logical array with 
 * box geometry containing all the existing cells of the node, that is, cells 
 * with status \f$ > -4 \f$, using the smallest possible space.
 *
 * - \anchor Key_Concepts_Code_Space [d] **Code Space**: In a typical 
 * simulation, the user chooses the coarsest level of refinement \f$ l_{min}\f$, 
 * the maximum level \f$ l_{max}\f$, and the simulation box with any length 
 * unit, and any coordinate system. But, the code always transforms the user box 
 * in the coarsest level of refinement \f$ l_{min}\f$ to fit in a cube of side 
 * equal to \f$ 2^{l_{min}} \f$ localized at the position (0,0,0) in a 
 * coordinate system, so the cube can be described as the set of points which 
 * belong to \f$ [0,2^{l_{min}})\times [0,2^{l_{min}})\times [0,2^{l_{min}})\f$. 
 * The 3-Dimensional space which goes from \f$ (-\infty,\infty)\times (-\infty,
 * \infty)\times (-\infty,\infty)\f$ contains this cube of side \f$ 2^{l_{min}} 
 * \f$ localized at (0,0,0) is called the \f${\color{red} \mathbf{ Code\ Space\ 
 * of\ refinement\ l_{min} }}\f$. For any level of refinement *l*, the box is 
 * localized in the same coordinate (0,0,0) but using a space equal to \f$ 
 * [0,2^l)^3\f$.
 *
 * - [e] \anchor Key_Concepts_Logical_Space **Logical Space**: Because of the 
 * enormous memory space required to store a full level of refinement *l* of 
 * size \f$ (2^{l})^3 \f$, the concept of boxes is implemented in every node. 
 * Boxes are logical entities that represent only a small piece of space of the 
 * *Code Space of refinement level l*, where *l* is the level of refinement of 
 * the node. This representation of the space is called the \f${\color{red} 
 * \mathbf{ Logical\ Space}}\f$. For every node at every level of refinement, 
 * this *Logical Space* has dimensions of \f$ [0,A)\times[0,B)\times[0,C) \f$, 
 * where *A, B, C* are the dimensions of the box of the node.
 *
 * **LONG DESCRIPTION**:
 *
 * The node structure is a Stand-alone unit to other nodes in the same level of
 * refinement. They can be almost completely updated without the use of any
 * other exterior information except for the boundary conditions which come from
 * external agents, such as its parent node or the user's initial conditions.
 * The basic idea of the nodes is to represent a refinement zone of the parent
 * node. So, as these zones live in the parent node, its refinement cells live
 * in the child node as if they were normal cells.
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

struct node
{
  //* >> Global node properties *//

  int ID; /**< Number \f$(\geq 0)\f$ of identification*/
  int lv; /**< Level of refinement*/

  //* >> Boxes *//
  int *ptr_box;       /**< Array containing the Box of the cells status */
  int *ptr_box_old;   /**< Array containing the old box of the cell status of a previous time-step, and used to adapt the box ::ptr_box to a new time-step */
  int box_cap;        /**< *Capacity* (see Key Concepts \ref Key_Concepts_Capacity "Capacity") of the boxes ::ptr_box and ::ptr_box_old */
  int box_real_dim_x; /**< Dimension X of the boxes ::ptr_box and ::ptr_box_old */
  int box_real_dim_y; /**< Same as ::box_real_dim_x but in the Y direction */
  int box_real_dim_z; /**< Same as ::box_real_dim_x but in the Z direction */
  // int box_real_dim_x_old; // Auxiliary real dimension X of the box
  // int box_real_dim_y_old; // Auxiliary real dimension X of the box
  // int box_real_dim_z_old; // Auxiliary real dimension X of the box

  int box_dim_x; /**< Dimension X of the *Smallest Box* (see Key Concepts \ref Key_Concepts_Smallest_Box "Smallest Box") of the boxes ::ptr_box and ::ptr_box_old */
  int box_dim_y; /**< Same as ::box_dim_x but in the Y direction */
  int box_dim_z; /**< Same as ::box_dim_x but in the Z direction */

  int box_ts_x;  /**< Translation constant between the cell index in the *Code Space* and the box index in the *Logical Space* (see Key Concepts \ref Key_Concepts_Code_Space "Code Space", and \ref Key_Concepts_Logical_Space "Logical Space") in the X direction */
  int box_ts_y;  /**< Same as ::box_ts_x but in the Y direction */
  int box_ts_z;  /**< Same as ::box_ts_x but in the Z direction */
  // int box_ts_x_old;       // Auxiliary index translation from real local index cell to box index at dimension X
  // int box_ts_y_old;       // Auxiliary ndex translation from real local index cell to box index at dimension Y
  // int box_ts_z_old;       // Auxiliary ndex translation from real local index cell to box index at dimension Z
  int box_min_x;      /**< Minimum cell index at the dimension X in the *Code Space* (see Key Concepts \ref Key_Concepts_Code_Space "Code Space") that defines the *Smallest Box* (see Key Concepts \ref Key_Concepts_Smallest_Box "Smallest Box") associated with the boxes ::ptr_box and ::ptr_box_old */
  int box_min_y;      /**< Same as ::box_min_x but in the Y direction */
  int box_min_z;      /**< Same as ::box_min_x but in the Z direction */
  int box_max_x;      /**< Maximum cell index at the dimension X in the *Code Space* (see Key Concepts \ref Key_Concepts_Code_Space "Code Space") that defines the *Smallest Box* (see Key Concepts \ref Key_Concepts_Smallest_Box "Smallest Box") associated with the boxes ::ptr_box and ::ptr_box_old */
  int box_max_y;      /**< Same as ::box_max_x but in the Y direction */
  int box_max_z;      /**< Same as ::box_max_x but in the Z direction */
  bool box_check_fit; /**< Check if the new box fits in the old one */

  //* >> Cells in the node *//
  int *ptr_cell_idx_x; /**< Cell index at X in the *Code Space* (see Key Concepts \ref Key_Concepts_Code_Space "Code Space") of the cells in the node */
  int *ptr_cell_idx_y; /**< Same as ::ptr_cell_idx_x but in the Y direction */
  int *ptr_cell_idx_z; /**< Same as ::ptr_cell_idx_x but in the Z direction */
  int *ptr_box_idx;    /**< Box index in the *Logical Space* (see Key Concepts \ref Key_Concepts_Logical_Space "Logical space") of the cells in the node  */
  int cell_cap;        /**< *Capacity* (see Key Concepts \ref Key_Concepts_Capacity "Capacity") of cell arrays ::ptr_cell_idx_x (\ref ::ptr_cell_idx_y "y", \ref ::ptr_cell_idx_z "z"), and ::ptr_box_idx   */
  int cell_size;       /**< *Size* (see Key Concepts \ref Key_Concepts_Size "Size") of cell arrays ::ptr_cell_idx_x (\ref ::ptr_cell_idx_y "y", \ref ::ptr_cell_idx_z "z"), and ::ptr_box_idx */

  //* >> Struct of cells (Particles and cell mass)
  struct cell_struct *ptr_cell_struct; /**< Array containing the cell structures of the node */
  struct cell_struct *ptr_cell_struct_old; /**< Array containing the old cell structures of the node, and used to adapt the cell structure array ::ptr_cell_struct to a new time-step */
  int cell_struct_old_cap;                 /**< *Capacity* (see Key Concepts \ref Key_Concepts_Capacity "Capacity") of the cell structure arrays ::ptr_cell_struct and ::ptr_cell_struct_old */

  //* >> Total mass in the node
  vtype node_mass; /**< Total mass of the particles inside of the node */

  //* >> Total number of particles in the node
  int no_ptcl_full_node;                       /**< Total number of particles in the node*/
  int no_ptcl_outs_ref_zones; /**< Total number of particles outisde of the refinement zones of the node */

  //* >> Grid points *//
  int *ptr_intr_grid_cell_idx_x; /**< Grid cell index at X in the *Code Space* (see Key Concepts \ref Key_Concepts_Code_Space "Code Space") of the interior grid points in the node */
  int *ptr_intr_grid_cell_idx_y; /**< Same as ::ptr_intr_grid_cell_idx_x but in the Y direction */
  int *ptr_intr_grid_cell_idx_z; /**< Same as ::ptr_intr_grid_cell_idx_x but in the Z direction */
  int *ptr_intr_box_grid_idx;    /**< Box grid index in the *Logical Space* (see Key Concepts \ref Key_Concepts_Logical_Space "Logical space") of the interior grid points in the node */
  int grid_intr_cap;             /**< *Capacity* (see Key Concepts \ref Key_Concepts_Capacity "Capacity") of the interior grid points arrays ::ptr_intr_grid_cell_idx_x (\ref ::ptr_bdry_grid_cell_idx_y "y", \ref ::ptr_bdry_grid_cell_idx_z "z"), and ::ptr_intr_box_grid_idx */
  int grid_intr_size;            /**< *Size* (see Key Concepts \ref Key_Concepts_Size "Size") of the interior grid points arrays ::ptr_intr_grid_cell_idx_x (\ref ::ptr_bdry_grid_cell_idx_y "y", \ref ::ptr_bdry_grid_cell_idx_z "z"), and ::ptr_intr_box_grid_idx */

  int *ptr_bdry_grid_cell_idx_x; /**< Grid cell index at X in the *Code Space* (see Key Concepts \ref Key_Concepts_Code_Space "Code Space") of the boundary grid points in the node */
  int *ptr_bdry_grid_cell_idx_y; /**< Same as ::ptr_bdry_grid_cell_idx_x but in the Y direction */
  int *ptr_bdry_grid_cell_idx_z; /**< Same as ::ptr_bdry_grid_cell_idx_x but in the Z direction */
  int *ptr_bdry_box_grid_idx;    /**< Box grid index in the *Logical Space* (see Key Concepts \ref Key_Concepts_Logical_Space "Logical space") of the boundary grid points in the node */
  int grid_bdry_cap;             /**< *Capacity* (see Key Concepts \ref Key_Concepts_Capacity "Capacity") of the boundary points arrays ::ptr_bdry_grid_cell_idx_x (\ref ::ptr_bdry_grid_cell_idx_y "y", \ref ::ptr_bdry_grid_cell_idx_z "z"), and ::ptr_bdry_box_grid_idx */
  int grid_bdry_size;            /**< *Size* (see Key Concepts \ref Key_Concepts_Size "Size") of the boundary points arrays ::ptr_bdry_grid_cell_idx_x (\ref ::ptr_bdry_grid_cell_idx_y "y", \ref ::ptr_bdry_grid_cell_idx_z "z"), and ::ptr_bdry_box_grid_idx */

  int *ptr_sim_bdry_grid_cell_idx_x; /**< Grid cell index at X in the *Code Space* (see Key Concepts \ref Key_Concepts_Code_Space "Code Space") of the *boundary simulation* grid points in the node */
  int *ptr_sim_bdry_grid_cell_idx_y; /**< Same as ::ptr_sim_bdry_grid_cell_idx_x but in the Y direction */
  int *ptr_sim_bdry_grid_cell_idx_z; /**< Same as ::ptr_sim_bdry_grid_cell_idx_x but in the Z direction */
  int *ptr_sim_bdry_box_grid_idx;    /**< Box grid index in the *Logical Space* (see Key Concepts \ref Key_Concepts_Logical_Space "Logical space") of the *boundary simulation* grid points in the node */
  int grid_sim_bdry_cap;             /**< *Capacity* (see Key Concepts \ref Key_Concepts_Capacity "Capacity") of the *boundary simulation* points arrays ::ptr_sim_bdry_grid_cell_idx_x (\ref ::ptr_sim_bdry_grid_cell_idx_y "y", \ref ::ptr_sim_bdry_grid_cell_idx_z "z"), and ::ptr_sim_bdry_box_grid_idx */
  int grid_sim_bdry_size;            /**< *Size* (see Key Concepts \ref Key_Concepts_Size "Size") of the *boundary simulation* points arrays ::ptr_sim_bdry_grid_cell_idx_x (\ref ::ptr_sim_bdry_grid_cell_idx_y "y", \ref ::ptr_sim_bdry_grid_cell_idx_z "z"), and ::ptr_sim_bdry_box_grid_idx */

  //* Potential, Acceleration, and density of the grid *//
  vtype *ptr_pot; /**< Array containing the potential of the box grid points in the node */
  vtype *ptr_pot_old;      /**< Array containing the previous time-step potential of the box grid points in the node, used to compute the error of the solution of the Poisson equation */
  vtype *ptr_ax;           /**< Array containing the acceleration in the X direction of the box grid points in the node  */
  vtype *ptr_ay;           /**< Same as ::ptr_ax but in the Y direction */
  vtype *ptr_az;           /**< Same as ::ptr_ax but in the Z direction */
  vtype *ptr_d;            /**< Array containing the density of the box grid points in the node */
  int grid_properties_cap; /**< *Capacity* (see Key Concepts \ref Key_Concepts_Capacity "Capacity") of the grid potentials ::ptr_pot, ::ptr_pot_old, accelerations ::ptr_ax, ::ptr_ay, ::ptr_az, and density ptr_d */

  //* >> Tree structure *//
  struct node **pptr_chn; /**< Pointer to child nodes pointer array */
  struct node *ptr_pt;    /**< Pointer to parent node */
  int chn_cap;            /**< *Capacity* (see Key Concepts \ref Key_Concepts_Capacity "Capacity") of the child nodes array ::pptr_chn */
  int chn_size;           /**< *Size* (see Key Concepts \ref Key_Concepts_Size "Size") of the child nodes array ::pptr_chn  */

  //* >> Auxiliary arrays to go from old box to new box *//
  int *ptr_cell_ref;  /**< Auxiliary array containing the positional indices of the node's cell arrays, for example ::ptr_cell_idx_x, which will require refinement */
  int cell_ref_cap;   /**< *Capacity* (see Key Concepts \ref Key_Concepts_Capacity "Capacity") of the auxiliary array of positional indices of cell to refine ::ptr_cell_ref */
  int cell_ref_size;  /**< *Size* (see Key Concepts \ref Key_Concepts_Size "Size") of the auxiliary array of positional indices of cell to refine ::ptr_cell_ref */
  int **pptr_zones;   /**< Auxiliary pointer to the array of pointers to the refined zones. Each zone contains the same type of positional indies as the auxiliary array ::ptr_cell_ref */
  int zones_cap;      /**< *Capacity* (see Key Concepts \ref Key_Concepts_Capacity "Capacity") of the refined zones auxiliary pointer array ::pptr_zones */
  int zones_size;     /**< *Size* (see Key Concepts \ref Key_Concepts_Size "Size") of the refined zones auxiliary pointer array ::pptr_zones */
  int *ptr_zone_cap;  /**< Auxiliary array of the *Capacities* (see Key Concepts \ref Key_Concepts_Capacity "Capacity") of each refined zone pointed by the array ::pptr_zones */
  int *ptr_zone_size; /**< Auxiliary array of the *Sizes* (see Key Concepts \ref Key_Concepts_Size "Size") of each refined zone pointed by the array ::pptr_zones */

  int *ptr_aux_idx;                    /**< Auxiliary array, multiple uses */
  int aux_idx_cap;                     /**< *Capacity* (see Key Concepts \ref Key_Concepts_Capacity "Capacity") of the auxiliary array ::ptr_aux_idx */
  bool *ptr_pbc_bool_bdry_anomalies_x; /**< Only for periodic boundary conditions (pbc). Flag array saying if the corresponding refinement zone crosses the simulation box at the X axis. */
  bool *ptr_pbc_bool_bdry_anomalies_y; /**< Same as ::ptr_pbc_bool_bdry_anomalies_x but in the Y direction */
  bool *ptr_pbc_bool_bdry_anomalies_z; /**< Same as ::ptr_pbc_bool_bdry_anomalies_x but in the Z direction */
  int pbc_bool_bdry_anomalies_cap;     /**< *Capacity* (see Key Concepts \ref Key_Concepts_Capacity "Capacity") of the boolean anomaly arrays ::ptr_pbc_bool_bdry_anomalies_x (\ref ::ptr_pbc_bool_bdry_anomalies_y "y", \ref ::ptr_pbc_bool_bdry_anomalies_z "z") */

  // Sub zones for periodic boundary conditions
  // int **pptr_subzones;    // Pointer to refined subzones in the node
  int *ptr_pbc_min_subzones_x; /**< Only for periodic boundary conditions (pbc). Minimum cell index at the dimension X in the *Code Space* (see Key Concepts \ref Key_Concepts_Code_Space "Code Space") that defines the *Smallest Box* (see Key Concepts \ref Key_Concepts_Smallest_Box "Smallest Box") associated with the corresponding subzone */
  int *ptr_pbc_min_subzones_y; /**< Same as ::ptr_pbc_min_subzones_x but in the Y direction */
  int *ptr_pbc_min_subzones_z; /**< Same as ::ptr_pbc_min_subzones_x but in the Z direction */
  int *ptr_pbc_max_subzones_x; /**< Only for periodic boundary conditions (pbc). Maximum cell index at the dimension X in the *Code Space* (see Key Concepts \ref Key_Concepts_Code_Space "Code Space") that defines the *Smallest Box* (see Key Concepts \ref Key_Concepts_Smallest_Box "Smallest Box") associated with the corresponding subzone */
  int *ptr_pbc_max_subzones_y; /**< Same as ::ptr_pbc_max_subzones_x but in the Y direction */
  int *ptr_pbc_max_subzones_z; /**< Same as ::ptr_pbc_max_subzones_x but in the Z direction */
  int pbc_subzones_cap;        /**< *Capacity* (see Key Concepts \ref Key_Concepts_Capacity "Capacity") of minimum and maximum cell indices arrays ::ptr_pbc_min_subzones_x (\ref ::ptr_pbc_min_subzones_y "y", \ref ::ptr_pbc_min_subzones_z "z"), and ::ptr_pbc_max_subzones_x (\ref ::ptr_pbc_max_subzones_y "y", \ref ::ptr_pbc_max_subzones_z "z") */
  int pbc_subzones_size;       /**< *Size* (see Key Concepts \ref Key_Concepts_Size "Size") of minimum and maximum cell indices arrays ::ptr_pbc_min_subzones_x (\ref ::ptr_pbc_min_subzones_y "y", \ref ::ptr_pbc_min_subzones_z "z"), and ::ptr_pbc_max_subzones_x (\ref ::ptr_pbc_max_subzones_y "y", \ref ::ptr_pbc_max_subzones_z "z")*/

  //* >> Links in Tree adaptation *//
  int *ptr_links_old_ord_old; /**< Array of non-negative integer values sorted from least to greatest, representing the identification ID of the old zones to be matched with the new zones in the node */
  int *ptr_links_new_ord_new; /**< Array of non-negative integer values sorted from least to greatest, representing the identification ID of the new zones to be matched with the old zones in the node */
  int *ptr_links_new_ord_old; /**< Array of non-negative integer values sorted by the old zone order ::ptr_links_old_ord_old, representing the identification ID of the new zones to be matched with the old zones in the node */
  int *ptr_links_old_ord_new; /**< Array of non-negative integer values sorted by the new zone order ::ptr_links_new_ord_new, representing the identification ID of the old zones to be matched with the new zones in the node */
  int links_cap;              /**< *Capacity* (see Key Concepts \ref Key_Concepts_Capacity "Capacity") of the links arrays ::ptr_links_old_ord_old, ::ptr_links_new_ord_new, ::ptr_links_new_ord_old, and ::ptr_links_old_ord_new  */

  //* >> Boundary of the simulation box *//
  bool sim_bdry_contact;   /**< Flag saying if the node touch the simulation box boundary */
  bool sim_bdry_contact_x; /**< Flag saying if the node touch the simulation box boundary at the X axis */
  bool sim_bdry_contact_y; /**< Same as ::sim_bdry_contact_x but in the Y direction */
  bool sim_bdry_contact_z; /**< Same as ::sim_bdry_contact_x but in the Z direction */

  // The following parameters are only used for the Periodic boundary conditions (boundary_type = 0), pbc = periodic boundary conditions
  bool pbc_crosses_sim_box_bdry;   /**< Only for periodic boundary conditions (pbc). Flag saying if the node crosses the simulation box boundary, i.e. the node has at least one cell inside of the box and another outside of it */
  bool pbc_crosses_sim_box_bdry_x; /**< Only for periodic boundary conditions (pbc). Flag saying if the node crosses the simulation box boundary at the X axis, i.e. the node has at least one cell inside of the box and another outside of it in the X axis */
  bool pbc_crosses_sim_box_bdry_y; /**< Same as ::pbc_crosses_sim_box_bdry_x but in the Y direction */
  bool pbc_crosses_sim_box_bdry_z; /**< Same as ::pbc_crosses_sim_box_bdry_x but in the Z direction */

  bool pbc_crosses_whole_sim_box;   /**< Only for periodic boundary conditions (pbc). Flag saying if the node has at least one of its dimensions the same size as the simulation box at that refinement level */
  bool pbc_crosses_whole_sim_box_x; /**< Only for periodic boundary conditions (pbc). Flag saying if the node has the same size X dimension than the simulation box at that refinement level */
  bool pbc_crosses_whole_sim_box_y; /**< Same as ::pbc_crosses_whole_sim_box_x but in the Y direction */
  bool pbc_crosses_whole_sim_box_z; /**< Same as ::pbc_crosses_whole_sim_box_x but in the Z direction */

  /*
   * The following parameter pbc_correction_due_pbc_flag is a flag used to 
   * correct box indices used only in periodic boundary conditions when the box 
   * pass from a crosses the whole simulation box to just cross the boundary 
   * of the simulation and vice versa.
   */

  bool pbc_correction_due_pbc_flag; /**< Only for periodic boundary conditions (pbc). Flag saying if required to do corrections in the function moving_old_child_to_new_child() of the tree_adaptation.c module to the node cells because of changes in the pbc flags */
};

void initialize_node(struct node *ptr_head);

#endif
