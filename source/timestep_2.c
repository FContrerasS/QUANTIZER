/*
 * timestep_2.c
 *
 * Compute of the time-step considering velocity and acceleration
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

#include "timestep_2.h"

static vtype timestep_computation_2(const struct node *ptr_node, bool status)
{

    int no_ptcl; // Total number of particles in the node
    int lv;      // Level of refinement

    int ptcl_idx; // Particle grid_idx in the node

    vtype aux_dt;
    vtype H;

    vtype mydt; // The time-step of the node

    no_ptcl = ptr_node->ptcl_size;
    lv = ptr_node->lv;
    H = 1.0L / (1 << lv);

    mydt = _MAX_dt_;

    //** >> Case no more child, the node is a leaf **/
    if (ptr_node->chn_size == 0)
    {
        for (int i = 0; i < no_ptcl; i++)
        {
            //** >> grid_idx of the particles in the node **/
            ptcl_idx = ptr_node->ptr_ptcl[i];

            //** >> x-axis
            if (myabs(GL_ptcl_ax[ptcl_idx]) <= 1.0e-12)
            {
                if (myabs(GL_ptcl_vx[ptcl_idx]) > 1.0e-12)
                {
                    aux_dt = _CFL_ * H / myabs(GL_ptcl_vx[ptcl_idx]);
                    mydt = mydt < aux_dt ? mydt : aux_dt;
                }
            }
            else
            {
                aux_dt = (-myabs(GL_ptcl_vx[ptcl_idx]) + sqrt(GL_ptcl_vx[ptcl_idx] * GL_ptcl_vx[ptcl_idx] + 2 * myabs(GL_ptcl_ax[ptcl_idx]) * H * _CFL_)) / myabs(GL_ptcl_ax[ptcl_idx]);
                mydt = mydt < aux_dt ? mydt : aux_dt;
            }

            //** >> y-axis
            if (myabs(GL_ptcl_ay[ptcl_idx]) <= 1.0e-12)
            {
                if (myabs(GL_ptcl_vy[ptcl_idx]) > 1.0e-12)
                {
                    aux_dt = _CFL_ * H / myabs(GL_ptcl_vy[ptcl_idx]);
                    mydt = mydt < aux_dt ? mydt : aux_dt;
                }
            }
            else
            {
                aux_dt = (-myabs(GL_ptcl_vy[ptcl_idx]) + sqrt(GL_ptcl_vy[ptcl_idx] * GL_ptcl_vy[ptcl_idx] + 2 * myabs(GL_ptcl_ay[ptcl_idx]) * H * _CFL_)) / myabs(GL_ptcl_ay[ptcl_idx]);
                mydt = mydt < aux_dt ? mydt : aux_dt;
            }

            //** >> z-axis
            if (myabs(GL_ptcl_az[ptcl_idx]) <= 1.0e-12)
            {
                if (myabs(GL_ptcl_vz[ptcl_idx]) > 1.0e-12)
                {
                    aux_dt = _CFL_ * H / myabs(GL_ptcl_vz[ptcl_idx]);
                    mydt = mydt < aux_dt ? mydt : aux_dt;
                }
            }
            else
            {
                aux_dt = (-myabs(GL_ptcl_vz[ptcl_idx]) + sqrt(GL_ptcl_vz[ptcl_idx] * GL_ptcl_vz[ptcl_idx] + 2 * myabs(GL_ptcl_az[ptcl_idx]) * H * _CFL_)) / myabs(GL_ptcl_az[ptcl_idx]);
                mydt = mydt < aux_dt ? mydt : aux_dt;
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

            if(GL_ptcl_updating_flag[ptcl_idx] !=status)
            {
                //** >> grid_idx of the particles in the node **/
                ptcl_idx = ptr_node->ptr_ptcl[i];

                //** >> x-axis
                if (myabs(GL_ptcl_ax[ptcl_idx]) <= 1.0e-12)
                {
                    if (myabs(GL_ptcl_vx[ptcl_idx]) > 1.0e-12)
                    {
                        aux_dt = _CFL_ * H / myabs(GL_ptcl_vx[ptcl_idx]);
                        mydt = mydt < aux_dt ? mydt : aux_dt;
                    }
                }
                else
                {
                    aux_dt = (-myabs(GL_ptcl_vx[ptcl_idx]) + sqrt(GL_ptcl_vx[ptcl_idx] * GL_ptcl_vx[ptcl_idx] + 2 * myabs(GL_ptcl_ax[ptcl_idx]) * H * _CFL_)) / myabs(GL_ptcl_ax[ptcl_idx]);
                    mydt = mydt < aux_dt ? mydt : aux_dt;
                }

                //** >> y-axis
                if (myabs(GL_ptcl_ay[ptcl_idx]) <= 1.0e-12)
                {
                    if (myabs(GL_ptcl_vy[ptcl_idx]) > 1.0e-12)
                    {
                        aux_dt = _CFL_ * H / myabs(GL_ptcl_vy[ptcl_idx]);
                        mydt = mydt < aux_dt ? mydt : aux_dt;
                    }
                }
                else
                {
                    aux_dt = (-myabs(GL_ptcl_vy[ptcl_idx]) + sqrt(GL_ptcl_vy[ptcl_idx] * GL_ptcl_vy[ptcl_idx] + 2 * myabs(GL_ptcl_ay[ptcl_idx]) * H * _CFL_)) / myabs(GL_ptcl_ay[ptcl_idx]);
                    mydt = mydt < aux_dt ? mydt : aux_dt;
                }

                //** >> z-axis
                if (myabs(GL_ptcl_az[ptcl_idx]) <= 1.0e-12)
                {
                    if (myabs(GL_ptcl_vz[ptcl_idx]) > 1.0e-12)
                    {
                        aux_dt = _CFL_ * H / myabs(GL_ptcl_vz[ptcl_idx]);
                        mydt = mydt < aux_dt ? mydt : aux_dt;
                    }
                }
                else
                {
                    aux_dt = (-myabs(GL_ptcl_vz[ptcl_idx]) + sqrt(GL_ptcl_vz[ptcl_idx] * GL_ptcl_vz[ptcl_idx] + 2 * myabs(GL_ptcl_az[ptcl_idx]) * H * _CFL_)) / myabs(GL_ptcl_az[ptcl_idx]);
                    mydt = mydt < aux_dt ? mydt : aux_dt;
                }

                //** >> The status of the particle is changed from not updated to updated **/
                GL_ptcl_updating_flag[ptcl_idx] = status;
            }
        }
    }

    return mydt;
}

int timestep_2(vtype *ptr_dt)
{

    //** >> Time-step computing **/

    struct node *ptr_node = NULL;
    bool status; // Boolean value for the updating particles

    int no_pts;   // Number of parents in the cycle
    status = !GL_ptcl_updating_flag[0];

    vtype dt_min;   // minimum dt
    vtype aux_dt; // Auxiliary time

    dt_min = _MAX_dt_;

    for (int lv = GL_tentacles_level_max; lv > -1; lv--)
    {
        no_pts = GL_tentacles_size[lv];

        //** >> For cycle over parent nodes **/
        for (int i = 0; i < no_pts; i++)
        {
            ptr_node = GL_tentacles[lv][i];

            //** >> Computing the time step of the node
            aux_dt = timestep_computation_2(ptr_node,status);
            
            if (dt_min > aux_dt)
            {
                dt_min = aux_dt;
            }
        }
    }

    *ptr_dt = dt_min;

    ptr_node = NULL;

    return _SUCCESS_;
}