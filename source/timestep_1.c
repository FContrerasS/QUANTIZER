/*
 * timestep_1.c
 *
 * Compute of the time-step considering velocity
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

#include "timestep_1.h"

static vtype timestep_computation_1(const struct node *ptr_node, bool status)
{

    int no_ptcl; // Total number of particles in the node
    int lv;      // Level of refinement

    vtype aux_v_pow2; // Auxiliary square velocity

    vtype H;

    int ptcl_idx; // Particle grid_idx in the node

    vtype mydt; // The time-step of the node
    vtype myvmax;   // Maximum velocity of all particles

    no_ptcl = ptr_node->ptcl_size;
    lv = ptr_node->lv;
    H = 1.0L / (1 << lv);

    myvmax = 0; // Minium velocity designated by the user

    //** >> Case no more child, the node is a leaf **/
    if (ptr_node->chn_size == 0)
    {
        for (int i = 0; i < no_ptcl; i++)
        {
            //** >> grid_idx of the particles in the node **/
            ptcl_idx = ptr_node->ptr_ptcl[i];

            //** >> Velocity at x-axis
            aux_v_pow2 = GL_ptcl_vx[ptcl_idx] * GL_ptcl_vx[ptcl_idx];
            if (myvmax < aux_v_pow2)
            {
                myvmax = aux_v_pow2;
            }
            //** >> Velocity at y-axis
            aux_v_pow2 = GL_ptcl_vy[ptcl_idx] * GL_ptcl_vy[ptcl_idx];
            if (myvmax < aux_v_pow2)
            {
                myvmax = aux_v_pow2;
            }
            //** >> Velocity at z-axis
            aux_v_pow2 = GL_ptcl_vz[ptcl_idx] * GL_ptcl_vz[ptcl_idx];
            if (myvmax < aux_v_pow2)
            {
                myvmax = aux_v_pow2;
            }

            //** >> The status of the particle is changed from not updated to updated **/
            GL_ptcl_updating_flag[ptcl_idx] = status;
        }
    }
    //** >> Case there are more children, the node is a branch **/
    else
    {
        for (int i = 0; i < no_ptcl; i++)
        {

            //** >> grid_idx of the particles in the node **/
            ptcl_idx = ptr_node->ptr_ptcl[i];

            if (GL_ptcl_updating_flag[ptcl_idx] != status)
            {

                //** >> Velocity at x-axis
                aux_v_pow2 = GL_ptcl_vx[ptcl_idx] * GL_ptcl_vx[ptcl_idx];
                if (myvmax < aux_v_pow2)
                {
                    myvmax = aux_v_pow2;
                }
                //** >> Velocity at y-axis
                aux_v_pow2 = GL_ptcl_vy[ptcl_idx] * GL_ptcl_vy[ptcl_idx];
                if (myvmax < aux_v_pow2)
                {
                    myvmax = aux_v_pow2;
                }
                //** >> Velocity at z-axis
                aux_v_pow2 = GL_ptcl_vz[ptcl_idx] * GL_ptcl_vz[ptcl_idx];
                if (myvmax < aux_v_pow2)
                {
                    myvmax = aux_v_pow2;
                }

                //** >> The status of the particle is changed from not updated to updated **/
                GL_ptcl_updating_flag[ptcl_idx] = status;
            }
        }
    }

    myvmax = sqrt(myvmax);

    mydt = myvmax < _CFL_ * H / _MAX_dt_ ? _MAX_dt_ : _CFL_ * H / myvmax;

    return mydt;
}

int timestep_1(vtype *ptr_dt)
{

    //** >> Time-step computing **/

    struct node *ptr_node;
    bool status; // Boolean value for the updating particles

    int no_pts; // Number of parents in the cycle
    status = !GL_ptcl_updating_flag[0];

    vtype dt_min; // minimum dt
    vtype aux_dt; // Auxiliary time

    dt_min = _MAX_dt_;

    for (int lv = GL_tentacles_level_max; lv > -1; lv--)
    {
        no_pts = GL_tentacles_size[lv];

        //** >> For cycle over parent nodes **/
        for (int i = 0; i < no_pts; i++)
        {
            ptr_node = GL_tentacles_old[lv][i];

            //** >> Computing the time step of the node
            aux_dt = timestep_computation_1(ptr_node, status);

            if (dt_min > aux_dt)
            {
                dt_min = aux_dt;
            }
        }
    }

    *ptr_dt = dt_min;

    return _SUCCESS_;
}