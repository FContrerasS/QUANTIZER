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

    //** >> Cells in the node **/
    int *ptr_cell_idx_x; // X index position of the cells in the node at level l
    int *ptr_cell_idx_y; // Y index position of the cells in the node at level l
    int *ptr_cell_idx_z; // Z index position of the cells in the node at level l
    int cell_cap;        // Maximum capacity of the array of cells in the node
    int cell_size;       // Number of existing cells in the node

    //** >> Particles in the node **/
    int *ptr_ptcl; // The size is completely related with the level of refinement, the total number of particles, the number of cells in the node, and the refinement criteria
    int ptcl_cap;  // Maximum cap of the particles array in the node
    int ptcl_size; // Number of existing particles in the node

    //** >> Boxes **/
    int *ptr_box;   // Box contaning the cells status of the minimal box cells and more
    int *ptr_box_aux;   // Auxiliary box contaning used to adatp the box to a new time-step
    int box_cap;        // Maximum capacity of the box
    int box_real_dim_x; // Real dimension X of the box
    int box_real_dim_y; // Real dimension Y of the box
    int box_real_dim_z; // Real dimension Z of the box
    int box_dim_x;      // Dimension X of the box (new and old)
    int box_dim_y;      // Dimension Y of the box (new and old)
    int box_dim_z;      // Dimension Z of the box (new and old)
    int box_ts_x;       // Index translation from real local index cell to box index at dimension X
    int box_ts_y;       // Index translation from real local index cell to box index at dimension Y
    int box_ts_z;       // Index translation from real local index cell to box index at dimension Z
    int box_min_x;      // Already minimal box value index in the real local space at the dimension X
    int box_min_y;      // Already minimal box value index in the real local space at the dimension Y
    int box_min_z;      // Already minimal box value index in the real local space at the dimension Z
    int box_max_x;      // Already maximum box value index in the real local space at the dimension X
    int box_max_y;      // Already maximum box value index in the real local space at the dimension Y
    int box_max_z;      // Already maximum box value index in the real local space at the dimension Z

    //** >> Grid points **/
    int *ptr_grid_intr; // Indexes of the interior grid points of the block
    int *ptr_grid_bder; // Indexes of the border grid points of the block
    int grid_intr_cap;  // Maximum cap of the grid interior points array of the block
    int grid_bder_cap;  // Maximum cap of the grid border points array of the block
    int grid_intr_size; // Number of existing grid interior points in the block
    int grid_bder_size; // Number of existing grid border points in the block

    // Notes that the MIN and MAX values can change between diffents time-steps, but the Translation indexs always keep equal exept if there is a reallocation of the boxes

    //** >> Refinement Criterion **/
    vtype *ptr_box_mass; // Register of the mass in the cell to refinement criteria
    vtype *ptr_box_mass_aux; // Auxiliary mass box used to adapt the box mass
    vtype local_mass;    // Total mass in the node

    //* Potential, Acceleration and density of the grid **/
    vtype *ptr_pot; // Array with the potential of the node. It is of the same size than the real box grid points
    vtype *ptr_ax;  // Same as potential but with the acceleration
    vtype *ptr_ay;  // Same as potential but with the acceleration
    vtype *ptr_az;  // Same as potential but with the acceleration
    vtype *ptr_d;   // Array with the density grid of the node.

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
    int aux_idx_cap;
};

void initialize_node(struct node *ptr_head);

#endif

