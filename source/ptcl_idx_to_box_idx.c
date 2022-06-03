/*
 * ptcl_idx_to_box_idx.c
 *
 * Using the node and the particle index, the function return the corresponding
 * box index of the particle
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

#include "ptcl_idx_to_box_idx.h"

int ptcl_idx_to_box_idx(struct node *ptr_node, int ptcl_idx)
{

    int lv;      // Level of refinement

    vtype pos_x; // Particle position in the grid level
    vtype pos_y;
    vtype pos_z;

    int pos_x_floor; // floor of the particle position in the grid level
    int pos_y_floor;
    int pos_z_floor;

    int box_idx_x; // Box index in X direcction
    int box_idx_y; // Box index in Y direcction
    int box_idx_z; // Box index in Z direcction
    int box_idx;   // Box index

    lv = ptr_node->lv;

    //** >> Position of the particles in the grid level of the current node before updating**/
    pos_x = GL_ptcl_x[ptcl_idx] * (1 << lv);
    pos_y = GL_ptcl_y[ptcl_idx] * (1 << lv);
    pos_z = GL_ptcl_z[ptcl_idx] * (1 << lv);

    //** >> Floor of the particles positions in the grid level of the current node before updating**/
    pos_x_floor = (int)pos_x;
    pos_y_floor = (int)pos_y;
    pos_z_floor = (int)pos_z;

    //** >> Box index in the current node before updating **/
    box_idx_x = pos_x_floor - ptr_node->box_ts_x;
    box_idx_y = pos_y_floor - ptr_node->box_ts_y;
    box_idx_z = pos_z_floor - ptr_node->box_ts_z;
    box_idx = box_idx_x + box_idx_y * ptr_node->box_real_dim_x + box_idx_z * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;

    return box_idx;
}