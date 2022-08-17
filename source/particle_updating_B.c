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

//** >> Local Functions
static void computing_particles_updating_B_HEAD_ONLY(struct node *ptr_node, vtype dt);
static void computing_particles_updating_B(struct node *ptr_node, vtype dt);

static void computing_particles_updating_B_HEAD_ONLY(struct node *ptr_node, vtype dt)
{

    for (int i = 0; i < GL_no_ptcl_final; i++)
    {
        //** >> Updating the new velocity of the particle **/
        GL_ptcl_vx[i] += GL_ptcl_ax[i] * dt * 0.5;
        GL_ptcl_vy[i] += GL_ptcl_ay[i] * dt * 0.5;
        GL_ptcl_vz[i] += GL_ptcl_az[i] * dt * 0.5;
    }
    
}

static void computing_particles_updating_B(struct node *ptr_node, vtype dt)
{

    int ptcl_idx; // Particle grid_idx in the node
    int box_idx;  // Box index of the node cell

    if (ptr_node->chn_size == 0)
    {
        for (int cell_idx = 0; cell_idx < ptr_node->cell_size; cell_idx++)
        {
            box_idx = ptr_node->ptr_box_idx[cell_idx];

            for (int j = 0; j < ptr_node->ptr_cell_struct[box_idx].ptcl_size; j++)
            {
                ptcl_idx = ptr_node->ptr_cell_struct[box_idx].ptr_ptcl[j];

                //** >> Updating the new velocity of the particle **/
                GL_ptcl_vx[ptcl_idx] += GL_ptcl_ax[ptcl_idx] * dt * 0.5;
                GL_ptcl_vy[ptcl_idx] += GL_ptcl_ay[ptcl_idx] * dt * 0.5;
                GL_ptcl_vz[ptcl_idx] += GL_ptcl_az[ptcl_idx] * dt * 0.5;
            }
        }
    }
    else
    {
        for (int cell_idx = 0; cell_idx < ptr_node->cell_size; cell_idx++)
        {
            box_idx = ptr_node->ptr_box_idx[cell_idx];

            if (ptr_node->ptr_box[box_idx] == -3)
            {
                for (int j = 0; j < ptr_node->ptr_cell_struct[box_idx].ptcl_size; j++)
                {
                    ptcl_idx = ptr_node->ptr_cell_struct[box_idx].ptr_ptcl[j];
                    
                    //** >> Updating the new velocity of the particle **/
                    GL_ptcl_vx[ptcl_idx] += GL_ptcl_ax[ptcl_idx] * dt * 0.5;
                    GL_ptcl_vy[ptcl_idx] += GL_ptcl_ay[ptcl_idx] * dt * 0.5;
                    GL_ptcl_vz[ptcl_idx] += GL_ptcl_az[ptcl_idx] * dt * 0.5;
                }
            }
        }
    }
}

void particle_updating_B(vtype dt)
{

    // Velocity are updated by a "Corrector step"

    //** >> Particle updating A **/

    // number of parents of the level = GL_tentacles_size[lv]

    if(lmin < lmax)
    {

        for (int lv = GL_tentacles_level_max; lv > -1; lv--)
        {
            //** >> For cycle over parent nodes **/
            for (int i = 0; i < GL_tentacles_size[lv]; i++)
            {
                // ptr_node = GL_tentacles[lv][i];
                computing_particles_updating_B(GL_tentacles[lv][i], dt);
            }
        }
    }
    else
    {
        computing_particles_updating_B_HEAD_ONLY(GL_tentacles[0][0], dt);
    }

}