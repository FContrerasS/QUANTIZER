/*
 * initialize_node.c
 *
 * Initialize basic parameters of a node
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

#include "initialize_node.h"

void initialize_node(struct node *ptr_node)
{
	//** >> Global node properties **/
	ptr_node->ID = -1; // Node ID
	ptr_node->lv = -1;  // Level of refinement

	//** >> Boxes **/
	ptr_node->ptr_box = NULL;	   // Box contaning the cells status of the minimal box cells and more
	ptr_node->ptr_box_old = NULL;  // Auxiliary box contaning used to adatp the box to a new time-step
	ptr_node->box_cap = 0;	//Maximum capacity of the box
	ptr_node->box_real_dim_x = 0;  // Real dimension X of the box
	ptr_node->box_real_dim_y = 0;  // Real dimension Y of the box
	ptr_node->box_real_dim_z = 0;  // Real dimension Z of the box
	ptr_node->box_real_dim_x_old = 0;  // Auxiliary real dimension X of the box
	ptr_node->box_real_dim_y_old = 0;  // Auxiliary real dimension Y of the box
	//ptr_node->box_real_dim_z_old = 0;  // Auxiliary real dimension Z of the box
	ptr_node->box_dim_x = 0;	   // Dimension X of the box (new and old)
	ptr_node->box_dim_y = 0;	   // Dimension Y of the box (new and old)
	ptr_node->box_dim_z = 0;	   // Dimension Z of the box (new and old)
	ptr_node->box_ts_x = 0;		   // Index translation from real local index cell to box index at dimension X
	ptr_node->box_ts_y = 0;		   // Index translation from real local index cell to box index at dimension Y
	ptr_node->box_ts_z = 0;		   // Index translation from real local index cell to box index at dimension Z
	ptr_node->box_ts_x_old = 0;		   // Auxiliary index translation from real local index cell to box index at dimension X
	ptr_node->box_ts_y_old = 0;		   // Auxiliary index translation from real local index cell to box index at dimension Y
	ptr_node->box_ts_z_old = 0;		   // Auxiliary index translation from real local index cell to box index at dimension Z
	ptr_node->box_min_x = INT_MAX; // Already minimal box value index in the real local space at the dimension X
	ptr_node->box_min_y = INT_MAX; // Already minimal box value index in the real local space at the dimension Y
	ptr_node->box_min_z = INT_MAX; // Already minimal box value index in the real local space at the dimension Z
	ptr_node->box_max_x = INT_MIN;	   // Already maximum box value index in the real local space at the dimension X
	ptr_node->box_max_y = INT_MIN;	   // Already maximum box value index in the real local space at the dimension Y
	ptr_node->box_max_z = INT_MIN;	   // Already maximum box value index in the real local space at the dimension Z
	ptr_node->box_check_fit = true; // Check if the new box will fit in the old one

	//** >> Cells in the node **/
	ptr_node->ptr_cell_idx_x = NULL; // X index of the cells in the node
	ptr_node->ptr_cell_idx_y = NULL; // Y index of the cells in the node
	ptr_node->ptr_cell_idx_z = NULL; // Z index of the cells in the node
	ptr_node->ptr_box_idx = NULL;
	ptr_node->cell_cap = 0;	 // Maximum capacity of the array of cells in the node
	ptr_node->cell_size = 0; // Number of existing cells in the node

	//** >> Struct of cells (Particles and cell mass)
	ptr_node->ptr_cell_struct = NULL;
	ptr_node->ptr_cell_struct_old = NULL;
	ptr_node->cell_struct_old_cap = 0;

	//** >> Total mass in the node **/
	ptr_node->local_mass = 0; 

	//** >> Grid points **/
	ptr_node->ptr_intr_grid_cell_idx_x = NULL; // X index of the interior grid point
	ptr_node->ptr_intr_grid_cell_idx_y = NULL; // Y index position of the cells in the node at level l
	ptr_node->ptr_intr_grid_cell_idx_z = NULL; // Z index position of the cells in the node at level l
	ptr_node->ptr_intr_grid_idx = NULL;	  // Indexes of the interior grid points of the block

	ptr_node->ptr_bder_grid_cell_idx_x = NULL;		 // X index of the interior grid point
	ptr_node->ptr_bder_grid_cell_idx_y = NULL;		 // Y index position of the cells in the node at level l
	ptr_node->ptr_bder_grid_cell_idx_z = NULL;		 // Z index position of the cells in the node at level l
	ptr_node->ptr_bder_grid_idx = NULL;

	ptr_node->ptr_SIMULATION_BOUNDARY_grid_cell_idx_x = NULL; // X index of the boundary of the simulation grid point
	ptr_node->ptr_SIMULATION_BOUNDARY_grid_cell_idx_y = NULL; // Y index position of the cells in the node at level l
	ptr_node->ptr_SIMULATION_BOUNDARY_grid_cell_idx_z = NULL; // Z index position of the cells in the node at level l
	ptr_node->ptr_SIMULATION_BOUNDARY_grid_idx = NULL;

	// Indexes of the interior grid points of the block
    ptr_node->grid_intr_cap = 0;   // Maximum cap of the grid interior points array of the block
    ptr_node->grid_bder_cap = 0; // Maximum cap of the grid border points array of the block
    ptr_node->grid_intr_size = 0; // Number of existing grid interior points in the block
    ptr_node->grid_bder_size = 0; // Number of existing grid border points in the block

	ptr_node->grid_SIMULATION_BOUNDARY_cap = 0;  // Maximum cap of the boundary simulation grid points array of the block
	ptr_node->grid_SIMULATION_BOUNDARY_size = 0; // Number of the boundary simulation grid points array of the block

	//** >> Potential, acceleration and density of the grid **/
	ptr_node->ptr_pot = NULL; // Array with the potential of the node. It is of the same size than the real box grid points
	ptr_node->ptr_pot_old = NULL;
	ptr_node->ptr_ax = NULL;	  // Same as potential but with the acceleration
	ptr_node->ptr_ay = NULL;	  // Same as potential but with the acceleration
	ptr_node->ptr_az = NULL;	  // Same as potential but with the acceleration
	ptr_node->ptr_d = NULL;	  // Array with the density grid of the node.
	ptr_node->grid_properties_cap = 0; // Maximum cap of the grid properties

	//** >> Tree structure **/
	ptr_node->pptr_chn = NULL; // Pointer to children pointers
	ptr_node->ptr_pt = NULL;   // Pointer to parent node
	ptr_node->chn_cap = 0;	   // Maximum capacity in the number of children nodes
	ptr_node->chn_size = 0;	   // Number of children of the node

	//** >> Auxililary arrays to go from old box to new box **/
	ptr_node->ptr_cell_ref = NULL;	// Index of the cell to be refined in the node cells array
	ptr_node->cell_ref_cap = 0;		// capacity of the refined cell array
	ptr_node->cell_ref_size = 0;	// Number of cells to be refined
	ptr_node->pptr_zones = NULL;	// Pointer to refined zones in the node. Each zone contain the cell index of ptr_cell_idx_x,Y,Z of it
	ptr_node->zones_cap = 0;		// capacity in the number of refined zones in the node
	ptr_node->zones_size = 0;		// Number of refined zones in the node
	ptr_node->ptr_zone_cap = NULL;	// capacity of each refined zone
	ptr_node->ptr_zone_size = NULL; // Number of cells in each refined zone

	ptr_node->ptr_aux_idx = NULL; // Auxiliary array, is used in the initial tree structure to save the index of the boxes elements and other
	ptr_node->ptr_aux_bool_boundary_simulation_contact_x = NULL;	// Only for periodic boundary conditions. It sais if the corresponding refinement zone touches the boundary of the simulation box at X axis
	ptr_node->ptr_aux_bool_boundary_simulation_contact_y = NULL;	// Only for periodic boundary conditions. It sais if the corresponding refinement zone touches the boundary of the simulation box at Y axis
	ptr_node->ptr_aux_bool_boundary_simulation_contact_z = NULL;	// Only for periodic boundary conditions. It sais if the corresponding refinement zone touches the boundary of the simulation box at Z axis
	ptr_node->ptr_aux_bool_boundary_anomalies_x = NULL;	// Only for periodic boundary conditions. It sais if the corresponding refinement zone crosses the simulation at X axis
	ptr_node->ptr_aux_bool_boundary_anomalies_y = NULL;
	ptr_node->ptr_aux_bool_boundary_anomalies_z = NULL;
	ptr_node->aux_idx_cap = 0;
	//Sub zones for periodic boundary conditions
	ptr_node->pptr_subzones = NULL;	// Pointer to refined subzones in the node
	ptr_node->subzones_cap = 0;		// capacity in the number of refined subzones in the node
	ptr_node->subzones_size = 0;		// Number of refined subzones in the node
	ptr_node->ptr_subzone_cap = NULL;	// capacity of each refined subzone
	ptr_node->ptr_subzone_size = NULL; // Number of cells in each refined subzone

	//** >> Links in Tree adaptation **/
	ptr_node->ptr_links_old_ord_old = NULL;
	ptr_node->ptr_links_new_ord_old = NULL;
	ptr_node->ptr_links_old_ord_new = NULL;
	ptr_node->ptr_links_new_ord_new = NULL;
	ptr_node->links_cap = 0;

	//** >> Boundary of the simulation box **/
	ptr_node->boundary_simulation_contact = false;
	ptr_node->boundary_simulation_contact_x = false; 
	ptr_node->boundary_simulation_contact_y = false;
	ptr_node->boundary_simulation_contact_z = false;

	//The following parameters are only used for the Periodic boundary conditions (boundary_type = 0), pbc = periodic boundary conditions
	ptr_node->pbc_anomalies_due_to_the_boundary = false;	//Anomalies happend when there are at least 2 neighboring cells separated by the simulation box boundary 
	ptr_node->pbc_crosses_the_boundary_simulation_box_x = false;	// when one node croses the simulation box at X direcction
	ptr_node->pbc_crosses_the_boundary_simulation_box_y = false;
	ptr_node->pbc_crosses_the_boundary_simulation_box_z = false;
	ptr_node->pbc_crosses_the_whole_simulation_box_x = false;
	ptr_node->pbc_crosses_the_whole_simulation_box_y = false;
	ptr_node->pbc_crosses_the_whole_simulation_box_z = false;
}

