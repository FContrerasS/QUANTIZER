/*
 * particle_acceleration.c
 *
 * Compute the global particle acceleration
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

#include "particle_acceleration.h"

static void computing_particle_acceleration(const struct node *ptr_node)
{
    int lv;      // Level of refinement

    // int box_idx_x; // Box index in X direcction of the node cell
    // int box_idx_y; // Box index in Y direcction of the node cell
    // int box_idx_z; // Box index in Z direcction of the node cell
    int box_idx;   // Box index of the node cell

    int ptcl_idx; // Particle grid_idx in the node

    vtype pos_x; // Particle position in the grid level
    vtype pos_y;
    vtype pos_z;

    int pos_x_floor; // floor of the particle position in the grid level
    int pos_y_floor;
    int pos_z_floor;

    int box_grid_idx_x; // grid_idx component in the grid box
    int box_grid_idx_y;
    int box_grid_idx_z;
    int box_grid_idx; // Grid box grid_idx

    vtype w_x; // Weight component of the CIC method. It corresponds to the distance between the particle and the grid point
    vtype w_y;
    vtype w_z;

    vtype w[8]; // Weight of the CIC method

    lv = ptr_node->lv;


    if(lmin < lmax)
    {
        //** >> Case no more child, the node is a leaf **/
        if (ptr_node->chn_size == 0)
        {
            for (int cell_idx = 0; cell_idx < ptr_node->cell_size; cell_idx++)
            {
                box_idx = ptr_node->ptr_box_idx[cell_idx];
                // box_idx_x = ptr_node->ptr_cell_idx_x[cell_idx] - ptr_node->box_ts_x;
                // box_idx_y = ptr_node->ptr_cell_idx_y[cell_idx] - ptr_node->box_ts_y;
                // box_idx_z = ptr_node->ptr_cell_idx_z[cell_idx] - ptr_node->box_ts_z;
                // box_idx = box_idx_x + box_idx_y * ptr_node->box_real_dim_x + box_idx_z * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;
                for (int j = 0; j < ptr_node->ptr_cell_struct[box_idx].ptcl_size; j++)
                {
                    ptcl_idx = ptr_node->ptr_cell_struct[box_idx].ptr_ptcl[j];

                    //** >> Position of the particles in the grid level **/
                    pos_x = GL_ptcl_x[ptcl_idx] * (1 << lv);
                    pos_y = GL_ptcl_y[ptcl_idx] * (1 << lv);
                    pos_z = GL_ptcl_z[ptcl_idx] * (1 << lv);

                    //** >> Floor of the particles positions in the grid level **/
                    pos_x_floor = (int)pos_x;
                    pos_y_floor = (int)pos_y;
                    pos_z_floor = (int)pos_z;

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
                                GL_ptcl_ax[ptcl_idx] += ptr_node->ptr_ax[box_grid_idx] * w[ii + 2 * jj + 4 * kk];
                                GL_ptcl_ay[ptcl_idx] += ptr_node->ptr_ay[box_grid_idx] * w[ii + 2 * jj + 4 * kk];
                                GL_ptcl_az[ptcl_idx] += ptr_node->ptr_az[box_grid_idx] * w[ii + 2 * jj + 4 * kk];
                            }
                        }
                    }
                }
            }
        }
        //** >> Case there are more children, the node is a branch **/
        else
        {
            for (int cell_idx = 0; cell_idx < ptr_node->cell_size; cell_idx++)
            {
                box_idx = ptr_node->ptr_box_idx[cell_idx];
                // box_idx_x = ptr_node->ptr_cell_idx_x[cell_idx] - ptr_node->box_ts_x;
                // box_idx_y = ptr_node->ptr_cell_idx_y[cell_idx] - ptr_node->box_ts_y;
                // box_idx_z = ptr_node->ptr_cell_idx_z[cell_idx] - ptr_node->box_ts_z;
                // box_idx = box_idx_x + box_idx_y * ptr_node->box_real_dim_x + box_idx_z * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;

                if (ptr_node->ptr_box[box_idx] < 0)
                {
                    for (int j = 0; j < ptr_node->ptr_cell_struct[box_idx].ptcl_size; j++)
                    {
                        ptcl_idx = ptr_node->ptr_cell_struct[box_idx].ptr_ptcl[j];

                        //** >> Position of the particles in the grid level **/
                        pos_x = GL_ptcl_x[ptcl_idx] * (1 << lv);
                        pos_y = GL_ptcl_y[ptcl_idx] * (1 << lv);
                        pos_z = GL_ptcl_z[ptcl_idx] * (1 << lv);

                        //** >> Floor of the particles positions in the grid level **/
                        pos_x_floor = (int)pos_x;
                        pos_y_floor = (int)pos_y;
                        pos_z_floor = (int)pos_z;

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
                                    GL_ptcl_ax[ptcl_idx] += ptr_node->ptr_ax[box_grid_idx] * w[ii + 2 * jj + 4 * kk];
                                    GL_ptcl_ay[ptcl_idx] += ptr_node->ptr_ay[box_grid_idx] * w[ii + 2 * jj + 4 * kk];
                                    GL_ptcl_az[ptcl_idx] += ptr_node->ptr_az[box_grid_idx] * w[ii + 2 * jj + 4 * kk];
                                }
                            }
                        }
                    }
                }

            }
        }
    }
    else    //Case only 1 lv of refinement
    {
        for (int i = 0; i < GL_no_ptcl; i++)
        {
            //** >> Position of the particles in the grid level **/
            pos_x = GL_ptcl_x[i] * (1 << lv);
            pos_y = GL_ptcl_y[i] * (1 << lv);
            pos_z = GL_ptcl_z[i] * (1 << lv);

            //** >> Floor of the particles positions in the grid level **/
            pos_x_floor = (int)pos_x;
            pos_y_floor = (int)pos_y;
            pos_z_floor = (int)pos_z;

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
                        GL_ptcl_ax[i] += ptr_node->ptr_ax[box_grid_idx] * w[ii + 2 * jj + 4 * kk];
                        GL_ptcl_ay[i] += ptr_node->ptr_ay[box_grid_idx] * w[ii + 2 * jj + 4 * kk];
                        GL_ptcl_az[i] += ptr_node->ptr_az[box_grid_idx] * w[ii + 2 * jj + 4 * kk];
                    }
                }
            }
        }
    }

}

    int
    particle_acceleration()
{

    //** >> Particle acceleration **/

    // number of parents of the level = GL_tentacles_size[lv]

    for (int lv = GL_tentacles_level_max; lv > -1; lv--)
    {
        //** >> For cycle over parent nodes **/
        for (int i = 0; i < GL_tentacles_size[lv]; i++)
        {
            // ptr_node = GL_tentacles[lv][i];

            computing_particle_acceleration(GL_tentacles[lv][i]);
        }
    }

    return _SUCCESS_;
}