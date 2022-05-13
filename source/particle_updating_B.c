/*
 * particle_updating_B.c
 *
 * Upadte global velocity particles by a "Corrector step"
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

#include "particle_updating_B.h"

static int computing_particles_updating_B(struct node *ptr_node, vtype dt, bool status)
{

    int ptcl_idx; // Particle grid_idx in the node

    int box_idx_x; // Box index in X direcction of the node cell
    int box_idx_y; // Box index in Y direcction of the node cell
    int box_idx_z; // Box index in Z direcction of the node cell
    int box_idx;   // Box index of the node cell

    clock_t aux_clock;

    if(lmin < lmax)
    {
        if (ptr_node->chn_size == 0)
        {
            for (int cell_idx = 0; cell_idx < ptr_node->cell_size; cell_idx++)
            {
                box_idx_x = ptr_node->ptr_cell_idx_x[cell_idx] - ptr_node->box_ts_x;
                box_idx_y = ptr_node->ptr_cell_idx_y[cell_idx] - ptr_node->box_ts_y;
                box_idx_z = ptr_node->ptr_cell_idx_z[cell_idx] - ptr_node->box_ts_z;
                box_idx = box_idx_x + box_idx_y * ptr_node->box_real_dim_x + box_idx_z * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;

                for (int j = 0; j < ptr_node->ptr_cell_struct[box_idx].ptcl_size; j++)
                {
                    ptcl_idx = ptr_node->ptr_cell_struct[box_idx].ptr_ptcl[j];

                    //** >> Updating the new velocity of the particle **/
                    GL_ptcl_vx[ptcl_idx] += GL_ptcl_ax[ptcl_idx] * dt / 2;
                    GL_ptcl_vy[ptcl_idx] += GL_ptcl_ay[ptcl_idx] * dt / 2;
                    GL_ptcl_vz[ptcl_idx] += GL_ptcl_az[ptcl_idx] * dt / 2;

                    //** >> The status of the particle is changed from not updated to updated **/
                    GL_ptcl_updating_flag[ptcl_idx] = status;
                }
            }
        }
        else
        {
            for (int cell_idx = 0; cell_idx < ptr_node->cell_size; cell_idx++)
            {

                box_idx_x = ptr_node->ptr_cell_idx_x[cell_idx] - ptr_node->box_ts_x;
                box_idx_y = ptr_node->ptr_cell_idx_y[cell_idx] - ptr_node->box_ts_y;
                box_idx_z = ptr_node->ptr_cell_idx_z[cell_idx] - ptr_node->box_ts_z;
                box_idx = box_idx_x + box_idx_y * ptr_node->box_real_dim_x + box_idx_z * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;
                for (int j = 0; j < ptr_node->ptr_cell_struct[box_idx].ptcl_size; j++)
                {
                    ptcl_idx = ptr_node->ptr_cell_struct[box_idx].ptr_ptcl[j];

                    if (GL_ptcl_updating_flag[ptcl_idx] != status)
                    {
                        //** >> Updating the new velocity of the particle **/
                        GL_ptcl_vx[ptcl_idx] += GL_ptcl_ax[ptcl_idx] * dt / 2;
                        GL_ptcl_vy[ptcl_idx] += GL_ptcl_ay[ptcl_idx] * dt / 2;
                        GL_ptcl_vz[ptcl_idx] += GL_ptcl_az[ptcl_idx] * dt / 2;

                        //** >> The status of the particle is changed from not updated to updated **/
                        GL_ptcl_updating_flag[ptcl_idx] = status;
                    }
                }
            }
        }
    }
    else // Case only 1 level of refinement, i.e. lmin = lmax
    {
        for (int i = 0; i < GL_no_ptcl; i++)
        {

            //** >> Updating the new velocity of the particle **/
            GL_ptcl_vx[i] += GL_ptcl_ax[i] * dt / 2;
            GL_ptcl_vy[i] += GL_ptcl_ay[i] * dt / 2;
            GL_ptcl_vz[i] += GL_ptcl_az[i] * dt / 2;


        }
    }



    return _SUCCESS_;
}


int particle_updating_B(vtype dt)
{

    // Velocity are updated by a "Corrector step"

    //** >> Particle updating A **/
    bool status;    // Boolean value for the updating particles

    // number of parents of the level = GL_tentacles_size[lv]

    status = !GL_ptcl_updating_flag[0];

    for (int lv = GL_tentacles_level_max; lv > -1; lv--)
    {
        //** >> For cycle over parent nodes **/
        for (int i = 0; i < GL_tentacles_size[lv]; i++)
        {
            //ptr_node = GL_tentacles[lv][i];

            if (computing_particles_updating_B(GL_tentacles[lv][i], dt, status) == _FAILURE_)
            {
                printf("Error at function computing_particles_updating_A()\n");
                return _FAILURE_;
            }
        }
    }

    return _SUCCESS_;
}