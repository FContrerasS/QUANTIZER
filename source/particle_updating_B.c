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

    int no_ptcl; // Total number of particles in the node

    int ptcl_idx; // Particle grid_idx in the node

    no_ptcl = ptr_node->ptcl_size;

    if (ptr_node->chn_size == 0)    
    {
        for (int i = 0; i < no_ptcl; i++)
        {
            ptcl_idx = ptr_node->ptr_ptcl[i];
            //** >> Updating the new velocity of the particle **/
            GL_ptcl_vx[ptcl_idx] += GL_ptcl_ax[ptcl_idx] * dt / 2;
            GL_ptcl_vy[ptcl_idx] += GL_ptcl_ay[ptcl_idx] * dt / 2;
            GL_ptcl_vz[ptcl_idx] += GL_ptcl_az[ptcl_idx] * dt / 2;

            //** >> The status of the particle is changed from not updated to updated **/
            GL_ptcl_updating_flag[ptcl_idx] = status;
            
        }
    }
    else
    {
        for (int i = 0; i < no_ptcl; i++)
        {
            ptcl_idx = ptr_node->ptr_ptcl[i];
            //** >> Particle has not been updated yet
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

    return _SUCCESS_;
}


int particle_updating_B(vtype dt)
{

    // Velocity are updated by a "Corrector step"

    //** >> Particle updating A **/

    struct node *ptr_node;
    bool status;    // Boolean value for the updating particles

    int no_pts; // Number of parents in the cycle
    status = !GL_ptcl_updating_flag[0];

    for (int lv = GL_tentacles_level_max; lv > -1; lv--)
    {
        no_pts = GL_tentacles_size[lv];
        //** >> For cycle over parent nodes **/
        for (int i = 0; i < no_pts; i++)
        {
            ptr_node = GL_tentacles_old[lv][i];

            if (computing_particles_updating_B(ptr_node, dt, status) == _FAILURE_)
            {
                printf("Error at function computing_particles_updating_A()\n");
                return _FAILURE_;
            }
        }
    }

    return _SUCCESS_;
}