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

static void computing_grid_density_head(struct node *ptr_head)
{
    int box_grid_idx_x; // grid_idx component in the grid box
    int box_grid_idx_y;
    int box_grid_idx_z;
    int box_grid_idx; // Grid box grid_idx
    int box_grid_idxNbr;   // Box index in the neigborhood

    vtype pos_x; // Particle position in the grid level
    vtype pos_y;
    vtype pos_z;

    int pos_x_floor; // floor of the particle position in the grid level
    int pos_y_floor;
    int pos_z_floor;

    vtype w_x_1; // Weight component of the CIC method. It corresponds to the distance between the particle and the grid point
    vtype w_y_1;
    vtype w_z_1;
    vtype w_x_2;
    vtype w_y_2;
    vtype w_z_2;

    vtype w[8]; // Weight of the CIC method

    int lv = ptr_head->lv;                              // Level of refinement
    vtype H = 1.0L / (1 << lv);                         // Size of the grid side
    vtype poisson_coeff = 4 * _G_ * _PI_ / (H * H * H); // Poisson coefficient

    int grid_box_real_dim_X = (ptr_head->box_real_dim_x + 1);
    int grid_box_real_dim_X_times_Y = (ptr_head->box_real_dim_x + 1) * (ptr_head->box_real_dim_y + 1);

    for (int i = 0; i < GL_no_ptcl; i++)
    {
        //** >> Position of the particles in the grid level **/
        pos_x = GL_ptcl_x[i] * (1 << lv);
        pos_y = GL_ptcl_y[i] * (1 << lv);
        pos_z = GL_ptcl_z[i] * (1 << lv);

        // ** >> floor of the particles positions in the grid level **/
        pos_x_floor = (int)pos_x;
        pos_y_floor = (int)pos_y;
        pos_z_floor = (int)pos_z;

        //** >> Computing the weights of the nearest grid points of the particle **/
        // Each for cyle yields 2 options: X or 1-X, where X =  pos_x - pos_x_floor

        w_x_1 = pos_x - pos_x_floor;
        w_y_1 = pos_y - pos_y_floor;
        w_z_1 = pos_z - pos_z_floor;
        w_x_2 = 1 - w_x_1;
        w_y_2 = 1 - w_y_1;
        w_z_2 = 1 - w_z_1;
        w[0] = w_x_2 * w_y_2 * w_z_2;
        w[1] = w_x_1 * w_y_2 * w_z_2;
        w[2] = w_x_2 * w_y_1 * w_z_2;
        w[3] = w_x_1 * w_y_1 * w_z_2;
        w[4] = w_x_2 * w_y_2 * w_z_1;
        w[5] = w_x_1 * w_y_2 * w_z_1;
        w[6] = w_x_2 * w_y_1 * w_z_1;
        w[7] = w_x_1 * w_y_1 * w_z_1;

        // for (int kk = 0; kk < 2; kk++)
        // {
        //     // w_z = kk + (1 - 2 * kk) * (pos_z - pos_z_floor);
        //     for (int jj = 0; jj < 2; jj++)
        //     {
        //         // w_y = jj + (1 - 2 * jj) * (pos_y - pos_y_floor);
        //         for (int ii = 0; ii < 2; ii++)
        //         {
        //             // w_x = ii + (1 - 2 * ii) * (pos_x - pos_x_floor);
        //             w[kk * 4 + jj * 2 + ii] = (1 - w_x) * (1 - w_y) * (1 - w_z);
        //         }
        //     }
        // }

        box_grid_idx_x = pos_x_floor - ptr_head->box_ts_x;
        box_grid_idx_y = pos_y_floor - ptr_head->box_ts_y;
        box_grid_idx_z = pos_z_floor - ptr_head->box_ts_z;
        box_grid_idx = box_grid_idx_x + box_grid_idx_y * grid_box_real_dim_X + box_grid_idx_z * grid_box_real_dim_X_times_Y;
        //** >> Particle density contributes to 8 enclosure grid points **/
        for (int kk = 0; kk < 2; kk++)
        {
            for (int jj = 0; jj < 2; jj++)
            {
                for (int ii = 0; ii < 2; ii++)
                {
                    // box_grid_idx_x = (pos_x_floor + ii) - ptr_head->box_ts_x;
                    // box_grid_idx_y = (pos_y_floor + jj) - ptr_head->box_ts_y;
                    // box_grid_idx_z = (pos_z_floor + kk) - ptr_head->box_ts_z;
                    //box_grid_idxNbr = box_grid_idx_x + box_grid_idx_y * (ptr_head->box_real_dim_x + 1) + box_grid_idx_z * (ptr_head->box_real_dim_x + 1) * (ptr_head->box_real_dim_y + 1);
                    box_grid_idxNbr = box_grid_idx + ii + jj * grid_box_real_dim_X + kk * grid_box_real_dim_X_times_Y;
                    ptr_head->ptr_d[box_grid_idxNbr] += poisson_coeff * GL_ptcl_mass[i] * w[ii + 2 * jj + 4 * kk];
                }
            }
        }
    }
}

static void computing_grid_density_branch(struct node *ptr_node)
{
    int box_grid_idx_x; // grid_idx component in the grid box
    int box_grid_idx_y;
    int box_grid_idx_z;
    int box_grid_idx; // Grid box grid_idx
    int box_grid_idxNbr; // Box index in the neigborhood

    int box_idx; // Box index of the node cell

    int ptcl_idx; // Particle grid_idx in the node

    vtype pos_x; // Particle position in the grid level
    vtype pos_y;
    vtype pos_z;

    int pos_x_floor; // floor of the particle position in the grid level
    int pos_y_floor;
    int pos_z_floor;

    vtype w_x_1; // Weight component of the CIC method. It corresponds to the distance between the particle and the grid point
    vtype w_y_1;
    vtype w_z_1;
    vtype w_x_2;
    vtype w_y_2;
    vtype w_z_2;

    vtype w[8]; // Weight of the CIC method

    int lv = ptr_node->lv;                              // Level of refinement
    vtype H = 1.0L / (1 << lv);                         // Size of the grid side
    vtype poisson_coeff = 4 * _G_ * _PI_ / (H * H * H); // Poisson coefficient

    int grid_box_real_dim_X = (ptr_node->box_real_dim_x + 1);
    int grid_box_real_dim_X_times_Y = (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1);

    for (int cell_idx = 0; cell_idx < ptr_node->cell_size; cell_idx++)
    {
        box_idx = ptr_node->ptr_box_idx[cell_idx];
        for (int j = 0; j < ptr_node->ptr_cell_struct[box_idx].ptcl_size; j++)
        {
            ptcl_idx = ptr_node->ptr_cell_struct[box_idx].ptr_ptcl[j];

            //** >> Position of the particles in the grid level **/
            pos_x = GL_ptcl_x[ptcl_idx] * (1 << lv);
            pos_y = GL_ptcl_y[ptcl_idx] * (1 << lv);
            pos_z = GL_ptcl_z[ptcl_idx] * (1 << lv);

            // ** >> floor of the particles positions in the grid level **/
            pos_x_floor = (int)pos_x;
            pos_y_floor = (int)pos_y;
            pos_z_floor = (int)pos_z;

            //** >> Computing the weights of the nearest grid points of the particle **/
            w_x_1 = pos_x - pos_x_floor;
            w_y_1 = pos_y - pos_y_floor;
            w_z_1 = pos_z - pos_z_floor;
            w_x_2 = 1 - w_x_1;
            w_y_2 = 1 - w_y_1;
            w_z_2 = 1 - w_z_1;
            w[0] = w_x_2 * w_y_2 * w_z_2;
            w[1] = w_x_1 * w_y_2 * w_z_2;
            w[2] = w_x_2 * w_y_1 * w_z_2;
            w[3] = w_x_1 * w_y_1 * w_z_2;
            w[4] = w_x_2 * w_y_2 * w_z_1;
            w[5] = w_x_1 * w_y_2 * w_z_1;
            w[6] = w_x_2 * w_y_1 * w_z_1;
            w[7] = w_x_1 * w_y_1 * w_z_1;

            //** >> Computing the weights of the nearest grid points of the particle **/
            // // Each for cyle yields 2 options: X or 1-X, where X =  pos_x - pos_x_floor
            // for (int kk = 0; kk < 2; kk++)
            // {
            //     w_z = kk + (1 - 2 * kk) * (pos_z - pos_z_floor);
            //     for (int jj = 0; jj < 2; jj++)
            //     {
            //         w_y = jj + (1 - 2 * jj) * (pos_y - pos_y_floor);
            //         for (int ii = 0; ii < 2; ii++)
            //         {
            //             w_x = ii + (1 - 2 * ii) * (pos_x - pos_x_floor);
            //             w[kk * 4 + jj * 2 + ii] = (1 - w_x) * (1 - w_y) * (1 - w_z);
            //         }
            //     }
            // }

            box_grid_idx_x = pos_x_floor - ptr_node->box_ts_x;
            box_grid_idx_y = pos_y_floor - ptr_node->box_ts_y;
            box_grid_idx_z = pos_z_floor - ptr_node->box_ts_z;
            box_grid_idx = box_grid_idx_x + box_grid_idx_y * grid_box_real_dim_X + box_grid_idx_z * grid_box_real_dim_X_times_Y;
            //** >> Particle density contributes to 8 enclosure grid points **/
            //** >> Particle density contributes to 8 enclosure grid points **/
            for (int kk = 0; kk < 2; kk++)
            {
                for (int jj = 0; jj < 2; jj++)
                {
                    for (int ii = 0; ii < 2; ii++)
                    {
                        // box_grid_idx_x = (pos_x_floor + ii) - ptr_node->box_ts_x;
                        // box_grid_idx_y = (pos_y_floor + jj) - ptr_node->box_ts_y;
                        // box_grid_idx_z = (pos_z_floor + kk) - ptr_node->box_ts_z;
                        // box_grid_idxNbr = box_grid_idx_x + box_grid_idx_y * (ptr_node->box_real_dim_x + 1) + box_grid_idx_z * (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1);
                        box_grid_idxNbr = box_grid_idx + ii + jj * grid_box_real_dim_X + kk * grid_box_real_dim_X_times_Y;
                        ptr_node->ptr_d[box_grid_idxNbr] += poisson_coeff * GL_ptcl_mass[ptcl_idx] * w[ii + 2 * jj + 4 * kk];
                    }
                }
            }
        }
    }
}

void grid_density()
{

    //** >> Density in the grid **/

    // Head
    computing_grid_density_head(GL_tentacles[0][0]);

    // Branches
    for (int lv = 1; lv < GL_tentacles_level_max + 1; lv++)
    {
        // number of parent of the level = GL_tentacles_size[lv];
        //** >> For cycle over parent nodes **/
        for (int i = 0; i < GL_tentacles_size[lv]; i++)
        {
            // ptr_node = GL_tentacles[lv][i];
            computing_grid_density_branch(GL_tentacles[lv][i]);
        }
    }
}