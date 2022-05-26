/*
 * initial_potential_and_acceleration_head.c
 *
 * Compute the potential and acceleartion of the head node grid points in the 
 * boundary of the simulation
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

#include "initial_potential_and_acceleration_head.h"

void initial_potential_and_acceleration_head(struct node *ptr_head)
{
    int grid_side;  // Number of grid points per side
    int grid_limit; // Superior border in the internal grid box

    vtype cm[3]; // Center of mass
    vtype dist;  // Distance between the center of mass and the border

    vtype H; // Size of the grid side

    int box_grid_idx; // Grid box index

    vtype pos_x; // Particle position in the grid level
    vtype pos_y;
    vtype pos_z;

    vtype aux_coeff; // Auxiliar coefficient

    grid_side = box_side_lmin + 1;
    grid_limit = grid_side - bder_os_sim - 1;

    H = 1.0L / no_lmin_cell;

    //** >> Computing the center of mass **/
    cm[0] = 0; // X position
    cm[1] = 0; // Y position
    cm[2] = 0; // Z position

    for (int i = 0; i < GL_no_ptcl; i++)
    {
        // Finding the grid position of the particles [0,Ng-1]
        pos_x = GL_ptcl_x[i] * no_lmin_cell;
        pos_y = GL_ptcl_y[i] * no_lmin_cell;
        pos_z = GL_ptcl_z[i] * no_lmin_cell;

        cm[0] += GL_ptcl_mass[i] * pos_x;
        cm[1] += GL_ptcl_mass[i] * pos_y;
        cm[2] += GL_ptcl_mass[i] * pos_z;
    }

    //** Normalizing per the total mass **/
    cm[0] = cm[0] / total_mass;
    cm[1] = cm[1] / total_mass;
    cm[2] = cm[2] / total_mass;

    //** Moving accoring the real grid box dimensions **/
    cm[0] += bder_os_sim;
    cm[1] += bder_os_sim;
    cm[2] += bder_os_sim;

    for (int k = bder_os_sim; k < grid_limit + 1; k++)
    {
        for (int j = bder_os_sim; j < grid_limit + 1; j++)
        {
            for (int i = bder_os_sim; i < grid_limit + 1; i++)
            {
                if (i == bder_os_sim || i == grid_limit || j == bder_os_sim || j == grid_limit || k == bder_os_sim || k == grid_limit)
                {
                    // Computing distance between the center of mass and the border point
                    dist = (i - cm[0]) * (i - cm[0]) + (j - cm[1]) * (j - cm[1]) + (k - cm[2]) * (k - cm[2]);
                    dist = H * sqrt(dist);
                    box_grid_idx = i + j * grid_side + k * grid_side * grid_side;
                    aux_coeff = -_G_ * total_mass / (dist * dist * dist) * H;
                    ptr_head->ptr_pot[box_grid_idx] = -_G_ * total_mass / dist;
                    ptr_head->ptr_ax[box_grid_idx] = aux_coeff * myabs(i - cm[0]);
                    ptr_head->ptr_ay[box_grid_idx] = aux_coeff * myabs(j - cm[1]);
                    ptr_head->ptr_az[box_grid_idx] = aux_coeff * myabs(k - cm[2]);
                }
            }
        }
    }
}