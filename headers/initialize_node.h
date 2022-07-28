/*
 * initialize_node.h
 *
 * Header file of the initialize_node.c source file
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

#ifndef __INITIALIZENODE__
#define __INITIALIZENODE__

#include "common.h"

struct node
{
    //** >> Global node properties **/
    int ID; // Node ID
    int lv;  // Level of refinement

    //** >> Boxes **/
    int *ptr_box;   // Box contaning the cells status of the minimal box cells and more
    int *ptr_box_old;   // Auxiliary box contaning used to adatp the box to a new time-step
    int box_cap;        // Maximum capacity of the box
    int box_real_dim_x; // Real dimension X of the box
    int box_real_dim_y; // Real dimension Y of the box
    int box_real_dim_z; // Real dimension Z of the box
    int box_real_dim_x_old; // Auxiliary real dimension X of the box
    int box_real_dim_y_old; // Auxiliary real dimension X of the box
    //int box_real_dim_z_old; // Auxiliary real dimension X of the box
    int box_dim_x;      // Dimension X of the box (new and old)
    int box_dim_y;      // Dimension Y of the box (new and old)
    int box_dim_z;      // Dimension Z of the box (new and old)
    int box_ts_x;       // Index translation from real local index cell to box index at dimension X
    int box_ts_y;       // Index translation from real local index cell to box index at dimension Y
    int box_ts_z;       // Index translation from real local index cell to box index at dimension Z
    int box_ts_x_old;       // Auxiliary index translation from real local index cell to box index at dimension X
    int box_ts_y_old;       // Auxiliary ndex translation from real local index cell to box index at dimension Y
    int box_ts_z_old;       // Auxiliary ndex translation from real local index cell to box index at dimension Z
    int box_min_x;      // Already minimal box value index in the real local space at the dimension X
    int box_min_y;      // Already minimal box value index in the real local space at the dimension Y
    int box_min_z;      // Already minimal box value index in the real local space at the dimension Z
    int box_max_x;      // Already maximum box value index in the real local space at the dimension X
    int box_max_y;      // Already maximum box value index in the real local space at the dimension Y
    int box_max_z;      // Already maximum box value index in the real local space at the dimension Z
    bool box_check_fit; // Check if the new box will fit in the old one

    //** >> Cells in the node **/
    int *ptr_cell_idx_x; // X index position of the cells in the node at level l
    int *ptr_cell_idx_y; // Y index position of the cells in the node at level l
    int *ptr_cell_idx_z; // Z index position of the cells in the node at level l
    int *ptr_box_idx;
    int cell_cap;  // Maximum capacity of the array of cells in the node
    int cell_size; // Number of existing cells in the node

    //** >> Struct of cells (Particles and cell mass)
    struct cell_struct *ptr_cell_struct;
    struct cell_struct *ptr_cell_struct_old;
    int cell_struct_old_cap;

    //** >> Total mass in the node
    vtype local_mass;        // Total mass in the node

    //** >> Grid points **/
    int *ptr_intr_grid_cell_idx_x; // X index of the interior grid point 
    int *ptr_intr_grid_cell_idx_y; // Y index position of the cells in the node at level l
    int *ptr_intr_grid_cell_idx_z; // Z index position of the cells in the node at level l
    int *ptr_intr_grid_idx;   // Indexes of the interior grid points of the block

    int *ptr_bder_grid_cell_idx_x; // X index of the interior grid point
    int *ptr_bder_grid_cell_idx_y; // Y index position of the cells in the node at level l
    int *ptr_bder_grid_cell_idx_z; // Z index position of the cells in the node at level l
    int *ptr_bder_grid_idx;       // Indexes of the interior grid points of the block

    int *ptr_SIMULATION_BOUNDARY_grid_cell_idx_x; // X index of the boundary of the simulation grid point
    int *ptr_SIMULATION_BOUNDARY_grid_cell_idx_y; // Y index position of the cells in the node at level l
    int *ptr_SIMULATION_BOUNDARY_grid_cell_idx_z; // Z index position of the cells in the node at level l
    int *ptr_SIMULATION_BOUNDARY_grid_idx;

    int grid_intr_cap;  // Maximum cap of the grid interior points array of the block
    int grid_bder_cap;  // Maximum cap of the grid border points array of the block
    int grid_intr_size; // Number of existing grid interior points in the block
    int grid_bder_size; // Number of existing grid border points in the block

    int grid_SIMULATION_BOUNDARY_cap;  // Maximum cap of the boundary simulation grid points array of the block
    int grid_SIMULATION_BOUNDARY_size; // Number of the boundary simulation grid points array of the block

    //* Potential, Acceleration and density of the grid **/
    vtype *ptr_pot; // Array with the potential of the node. It is of the same size than the real box grid points
    vtype *ptr_pot_old;
    vtype *ptr_ax;  // Same as potential but with the acceleration
    vtype *ptr_ay;  // Same as potential but with the acceleration
    vtype *ptr_az;  // Same as potential but with the acceleration
    vtype *ptr_d;   // Array with the density grid of the node.
    int grid_properties_cap;   // Maximum cap of the grid properties

    //** >> Tree structure **/
    struct node **pptr_chn; // Pointer to children pointers
    struct node *ptr_pt;    // Pointer to parent node
    int chn_cap;            // Maximum capacity in the number of children nodes
    int chn_size;           // Number of children of the node

    //** >> Auxililary arrays to go from old box to new box **/
    int *ptr_cell_ref;  // Index of the cell to be refined in the node cells array
    int cell_ref_cap;   // capacity of the refined cell array
    int cell_ref_size;  // Number of cells to be refined
    int **pptr_zones;   // Pointer to refined zones in the node. Each zone contain the cell index of ptr_cell_idx_X,Y,Z of it
    int zones_cap;      // capacity in the number of refined zones in the node
    int zones_size;     // Number of refined zones in the node
    int *ptr_zone_cap;  // capacity of each refined zone
    int *ptr_zone_size; // Number of cells in each refined zone

    int *ptr_aux_idx; // Auxiliary array, is used in the initial tree structure to save the index of the boxes elements and other
    bool *ptr_aux_bool_boundary_simulation_contact_x; // Only for periodic boundary conditions. It sais if the corresponding refinement zone touches the boundary of the simulation box at X axis
    bool *ptr_aux_bool_boundary_simulation_contact_y; // Only for periodic boundary conditions. It sais if the corresponding refinement zone touches the boundary of the simulation box at Y axis
    bool *ptr_aux_bool_boundary_simulation_contact_z; // Only for periodic boundary conditions. It sais if the corresponding refinement zone touches the boundary of the simulation box at Z axis
    bool *ptr_aux_bool_boundary_anomalies_x; // Only for periodic boundary conditions. It sais if the corresponding refinement zone crosses the simulation at X axis
    bool *ptr_aux_bool_boundary_anomalies_y;
    bool *ptr_aux_bool_boundary_anomalies_z;
    int aux_idx_cap;

    // Sub zones for periodic boundary conditions
    int **pptr_subzones;    // Pointer to refined subzones in the node
    int subzones_cap;        // capacity in the number of refined subzones in the node
    int subzones_size;       // Number of refined subzones in the node
    int *ptr_subzone_cap;  // capacity of each refined subzone
    int *ptr_subzone_size; // Number of cells in each refined subzone

    //** >> Links in Tree adaptation **/
    int *ptr_links_old_ord_old;
    int *ptr_links_new_ord_old;
    int *ptr_links_old_ord_new;
    int *ptr_links_new_ord_new;
    int links_cap;

    //** >> Boundary of the simulation box **/
    bool boundary_simulation_contact;
    bool boundary_simulation_contact_x;
    bool boundary_simulation_contact_y;
    bool boundary_simulation_contact_z;

    bool pbc_anomalies_due_to_the_boundary;
    bool pbc_crosses_the_boundary_simulation_box_x;
    bool pbc_crosses_the_boundary_simulation_box_y;
    bool pbc_crosses_the_boundary_simulation_box_z;
    bool pbc_crosses_the_whole_simulation_box_x;
    bool pbc_crosses_the_whole_simulation_box_y;
    bool pbc_crosses_the_whole_simulation_box_z;
};



void initialize_node(struct node *ptr_head);

#endif

