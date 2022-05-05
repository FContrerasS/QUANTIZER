/*
 * grid_density.c
 *
 * Compute the density of the grid in each node
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

#include "grid_density.h"

static void computing_grid_density(struct node *ptr_node)
{
    int no_ptcl;    // Number of particles in the node
    int lv; // Level of refinement

    int box_grid_idx_x; // grid_idx component in the grid box
    int box_grid_idx_y;
    int box_grid_idx_z;
    int box_grid_idx; // Grid box grid_idx

    int ptcl_idx; // Particle grid_idx in the node

    vtype pos_x; // Particle position in the grid level
    vtype pos_y;
    vtype pos_z;

    int pos_x_floor; // floor of the particle position in the grid level
    int pos_y_floor;
    int pos_z_floor;

    vtype H;             // Size of the grid side
    vtype poisson_coeff; // Poisson coefficient

    vtype w_x; // Weight component of the CIC method. It corresponds to the distance between the particle and the grid point
    vtype w_y;
    vtype w_z;

    vtype w[8]; // Weight of the CIC method

    no_ptcl = ptr_node->ptcl_size;
    lv = ptr_node->lv;
    H = 1.0L / (1 << lv);
    poisson_coeff = 4 * _G_ * _PI_ / (H * H * H);

    for (int i = 0; i < no_ptcl; i++)
    {

        //** >> grid_idx of the particles in the node **/
        ptcl_idx = ptr_node->ptr_ptcl[i];

        //** >> Position of the particles in the grid level **/
        pos_x = GL_ptcl_x[ptcl_idx] * (1 << lv);
        pos_y = GL_ptcl_y[ptcl_idx] * (1 << lv);
        pos_z = GL_ptcl_z[ptcl_idx] * (1 << lv);

        // ** >> floor of the particles positions in the grid level **/
        pos_x_floor = (int) pos_x;
        pos_y_floor = (int) pos_y;
        pos_z_floor = (int) pos_z;

        //** >> Computing the weights of the nearest grid points of the particle **/
        // Each for cyle yields 2 options: X or 1-X, where X =  pos_x - pos_x_floor
        for (int kk = 0; kk < 2; kk++)
        {
            w_z = kk + (1 - 2 * kk) * (pos_z - pos_z_floor);
            for (int jj = 0; jj < 2; jj++)
            {
                w_y = jj + (1 - 2 * jj) * (pos_y - pos_y_floor);
                for (int ii = 0; ii < 2; ii++)
                {
                    w_x = ii + (1 - 2 * ii) * (pos_x - pos_x_floor);
                    w[kk * 4 + jj * 2 + ii] = (1 - w_x) * (1 - w_y) * (1 - w_z);
                }
            }
        }
        //** >> Particle density contributes to 8 enclosure grid points **/
        for (int kk = 0; kk < 2; kk++)
        {
            for (int jj = 0; jj < 2; jj++)
            {
                for (int ii = 0; ii < 2; ii++)
                {
                    box_grid_idx_x = (pos_x_floor + ii) - ptr_node->box_ts_x;
                    box_grid_idx_y = (pos_y_floor + jj) - ptr_node->box_ts_y;
                    box_grid_idx_z = (pos_z_floor + kk) - ptr_node->box_ts_z;
                    box_grid_idx = box_grid_idx_x + box_grid_idx_y * (ptr_node->box_real_dim_x + 1) + box_grid_idx_z * (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1);
                    ptr_node->ptr_d[box_grid_idx] += poisson_coeff * GL_ptcl_mass[ptcl_idx] * w[ii + 2 * jj + 4 * kk];
                }
            }
        }
    }
}

void grid_density()
{
    
    //** >> Density in the grid **/

    for (int lv = 0; lv < lmax - lmin + 1; lv++)
    {
        //number of parent of the level = GL_tentacles_size[lv];
        //** >> For cycle over parent nodes **/
        for (int i = 0; i < GL_tentacles_size[lv]; i++)
        {
            //ptr_node = GL_tentacles[lv][i];
            computing_grid_density(GL_tentacles[lv][i]);
        }
        if (GL_tentacles_size[lv+1] == 0)
        {
            lv = lmax - lmin + 1;
        }
    }
}
