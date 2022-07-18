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

//** >> Local Functions
static vtype timestep_computation_1_HEAD_ONLY(const struct node *ptr_node);
static vtype timestep_computation_1(const struct node *ptr_node);

static vtype timestep_computation_1_HEAD_ONLY(const struct node *ptr_node)
{



    vtype aux_v_pow2; // Auxiliary square velocity
    vtype mydt;   // The time-step of the node
    vtype myvmax; // Maximum velocity of all particles

    vtype H = 1.0L / (1 << ptr_node->lv);

    myvmax = 0; // Minium velocity designated by the user

    for (int i = 0; i < GL_no_ptcl_final; i++)
    {
        //** >> Velocity at x-axis
        aux_v_pow2 = myabs(GL_ptcl_vx[i]);
        if (myvmax < aux_v_pow2)
        {
            myvmax = aux_v_pow2;
        }
        //** >> Velocity at y-axis
        aux_v_pow2 = myabs(GL_ptcl_vy[i]);
        if (myvmax < aux_v_pow2)
        {
            myvmax = aux_v_pow2;
        }
        //** >> Velocity at z-axis
        aux_v_pow2 = myabs(GL_ptcl_vz[i]);
        if (myvmax < aux_v_pow2)
        {
            myvmax = aux_v_pow2;
        }
    }

    mydt = myvmax < (_CFL_ * H / _MAX_dt_) ? _MAX_dt_ : (_CFL_ * H / myvmax);

    return mydt;
}

static vtype timestep_computation_1(const struct node *ptr_node)
{
    int lv;      // Level of refinement

    vtype aux_v_pow2; // Auxiliary square velocity

    vtype H;

    int box_idx;   // Box index of the node cell

    int ptcl_idx; // Particle grid_idx in the node

    vtype mydt; // The time-step of the node
    vtype myvmax;   // Maximum velocity of all particles

    lv = ptr_node->lv;
    H = 1.0L / (1 << lv);

    myvmax = 0; // Minium velocity designated by the user

    //** >> Case no more child, the node is a leaf **/
    if (ptr_node->chn_size == 0)
    {
        for (int cell_idx = 0; cell_idx < ptr_node->cell_size; cell_idx++)
        {
            box_idx = ptr_node->ptr_box_idx[cell_idx];

            for (int j = 0; j < ptr_node->ptr_cell_struct[box_idx].ptcl_size; j++)
            {
                ptcl_idx = ptr_node->ptr_cell_struct[box_idx].ptr_ptcl[j];

                //** >> Velocity at x-axis
                aux_v_pow2 = myabs(GL_ptcl_vx[ptcl_idx]);
                if (myvmax < aux_v_pow2)
                {
                    myvmax = aux_v_pow2;
                }
                //** >> Velocity at y-axis
                aux_v_pow2 = myabs(GL_ptcl_vy[ptcl_idx]);
                if (myvmax < aux_v_pow2)
                {
                    myvmax = aux_v_pow2;
                }
                //** >> Velocity at z-axis
                aux_v_pow2 = myabs(GL_ptcl_vz[ptcl_idx]);
                if (myvmax < aux_v_pow2)
                {
                    myvmax = aux_v_pow2;
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

            if (ptr_node->ptr_box[box_idx] < 0)
            {
                for (int j = 0; j < ptr_node->ptr_cell_struct[box_idx].ptcl_size; j++)
                {
                    ptcl_idx = ptr_node->ptr_cell_struct[box_idx].ptr_ptcl[j];

                    //** >> Velocity at x-axis
                    aux_v_pow2 = myabs(GL_ptcl_vx[ptcl_idx]);
                    if (myvmax < aux_v_pow2)
                    {
                        myvmax = aux_v_pow2;
                    }
                    //** >> Velocity at y-axis
                    aux_v_pow2 = myabs(GL_ptcl_vy[ptcl_idx]) ;
                    if (myvmax < aux_v_pow2)
                    {
                        myvmax = aux_v_pow2;
                    }
                    //** >> Velocity at z-axis
                    aux_v_pow2 = myabs(GL_ptcl_vz[ptcl_idx]);
                    if (myvmax < aux_v_pow2)
                    {
                        myvmax = aux_v_pow2;
                    }
                }
            }
        }
    }

    mydt = myvmax < (_CFL_ * H / _MAX_dt_) ? _MAX_dt_ : (_CFL_ * H / myvmax);

    return mydt;
}

int timestep_1(vtype *ptr_dt)
{

    //** >> Time-step computing **/

    vtype dt_min; // minimum dt
    vtype aux_dt; // Auxiliary time

    dt_min = _MAX_dt_;


    if(lmin < lmax)
    {
        for (int lv = GL_tentacles_level_max; lv > -1; lv--)
        {
            // NUmber of parents of the level = GL_tentacles_size[lv];

            //** >> For cycle over parent nodes **/
            for (int i = 0; i < GL_tentacles_size[lv]; i++)
            {
                // ptr_node = GL_tentacles[lv][i];

                //** >> Computing the time step of the node
                aux_dt = timestep_computation_1(GL_tentacles[lv][i]);

                if (dt_min > aux_dt)
                {
                    dt_min = aux_dt;
                }
            }
        }
    }
    else
    {
        aux_dt = timestep_computation_1_HEAD_ONLY(GL_tentacles[0][0]);
        if (dt_min > aux_dt)
        {
            dt_min = aux_dt;
        }
    }

    *ptr_dt = dt_min;

    return _SUCCESS_;
}