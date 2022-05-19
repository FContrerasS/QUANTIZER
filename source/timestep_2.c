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

static vtype timestep_computation_2(const struct node *ptr_node)
{
    int lv;      // Level of refinement

    int box_idx;   // Box index of the node cell
    int ptcl_idx; // Particle grid_idx in the node

    vtype aux_dt;
    vtype H;

    vtype mydt; // The time-step of the node

    lv = ptr_node->lv;
    H = 1.0L / (1 << lv);

    mydt = _MAX_dt_;

    if(lmin < lmax)
    {

        //** >> Case no more child, the node is a leaf **/
        if (ptr_node->chn_size == 0)
        {
            for (int cell_idx = 0; cell_idx < ptr_node->cell_size; cell_idx++)
            {
                box_idx = ptr_node->ptr_box_idx[cell_idx];

                for (int j = 0; j < ptr_node->ptr_cell_struct[box_idx].ptcl_size; j++)
                {
                    ptcl_idx = ptr_node->ptr_cell_struct[box_idx].ptr_ptcl[j];

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
                }
            }
        }
        //** >> Case there are more children, the node is a branch **/
        else
        {
            for (int cell_idx = 0; cell_idx < ptr_node->cell_size; cell_idx++)
            {
                box_idx = ptr_node->ptr_box_idx[cell_idx];
                
                if(ptr_node->ptr_box[box_idx] < 0)
                {
                    for (int j = 0; j < ptr_node->ptr_cell_struct[box_idx].ptcl_size; j++)
                    {
                        ptcl_idx = ptr_node->ptr_cell_struct[box_idx].ptr_ptcl[j];

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
                    }
                }
            }
        }
    }
    else // Case only 1 level of refinement, i.e. lmin = lmax
    {
        for (int i = 0; i < GL_no_ptcl; i++)
        {

            //** >> x-axis
            if (myabs(GL_ptcl_ax[i]) <= 1.0e-12)
            {
                if (myabs(GL_ptcl_vx[i]) > 1.0e-12)
                {
                    aux_dt = _CFL_ * H / myabs(GL_ptcl_vx[i]);
                    mydt = mydt < aux_dt ? mydt : aux_dt;
                }
            }
            else
            {
                aux_dt = (-myabs(GL_ptcl_vx[i]) + sqrt(GL_ptcl_vx[i] * GL_ptcl_vx[i] + 2 * myabs(GL_ptcl_ax[i]) * H * _CFL_)) / myabs(GL_ptcl_ax[i]);
                mydt = mydt < aux_dt ? mydt : aux_dt;
            }

            //** >> y-axis
            if (myabs(GL_ptcl_ay[i]) <= 1.0e-12)
            {
                if (myabs(GL_ptcl_vy[i]) > 1.0e-12)
                {
                    aux_dt = _CFL_ * H / myabs(GL_ptcl_vy[i]);
                    mydt = mydt < aux_dt ? mydt : aux_dt;
                }
            }
            else
            {
                aux_dt = (-myabs(GL_ptcl_vy[i]) + sqrt(GL_ptcl_vy[i] * GL_ptcl_vy[i] + 2 * myabs(GL_ptcl_ay[i]) * H * _CFL_)) / myabs(GL_ptcl_ay[i]);
                mydt = mydt < aux_dt ? mydt : aux_dt;
            }

            //** >> z-axis
            if (myabs(GL_ptcl_az[i]) <= 1.0e-12)
            {
                if (myabs(GL_ptcl_vz[i]) > 1.0e-12)
                {
                    aux_dt = _CFL_ * H / myabs(GL_ptcl_vz[i]);
                    mydt = mydt < aux_dt ? mydt : aux_dt;
                }
            }
            else
            {
                aux_dt = (-myabs(GL_ptcl_vz[i]) + sqrt(GL_ptcl_vz[i] * GL_ptcl_vz[i] + 2 * myabs(GL_ptcl_az[i]) * H * _CFL_)) / myabs(GL_ptcl_az[i]);
                mydt = mydt < aux_dt ? mydt : aux_dt;
            }
        }
    }


    return mydt;
}

int timestep_2(vtype *ptr_dt)
{

    //** >> Time-step computing **/

    vtype dt_min; // minimum dt
    vtype aux_dt; // Auxiliary time

    dt_min = _MAX_dt_;

    for (int lv = GL_tentacles_level_max; lv > -1; lv--)
    {
        // Number of parent of the level = GL_tentacles_size[lv];

        //** >> For cycle over parent nodes **/
        for (int i = 0; i < GL_tentacles_size[lv]; i++)
        {
            // ptr_node = GL_tentacles[lv][i];

            //** >> Computing the time step of the node
            aux_dt = timestep_computation_2(GL_tentacles[lv][i]);

            if (dt_min > aux_dt)
            {
                dt_min = aux_dt;
            }
        }
    }

    *ptr_dt = dt_min;

    return _SUCCESS_;
}