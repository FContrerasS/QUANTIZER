/*
 * grid_acceleration.c
 *
 * Compute acceleration in the grid of each node
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

#include "grid_acceleration.h"

static void computing_grid_acceleration(struct node *ptr_node)
{
    vtype H;
    vtype one_over_2H;
    vtype aux_acc; // Auxiliary value for the acceleration

    int lv;

    int cell_idx_x; // The cell index
    int cell_idx_y; // The cell index
    int cell_idx_z; // The cell index

    int box_idx_x; // Box grid index in X direcction
    int box_idx_y; // Box grid index in Y direcction
    int box_idx_z; // Box grid index in Z direcction
    int box_grid_idx;   // Box grid index

    int pt_box_idx_1_x; // Parent box index 1 at X direcction
    int pt_box_idx_1_y; // Parent box index 1 at Y direcction
    int pt_box_idx_1_z; // Parent box index 1 at Z direcction
    int pt_box_idx_2_x; // Parent box index 2 at X direcction
    int pt_box_idx_2_y; // Parent box index 2 at Y direcction
    int pt_box_idx_2_z; // Parent box index 2 at Z direcction
    int pt_box_grid_idx_1; // Parent box grid index 1
    int pt_box_grid_idx_2; // Parent box grid index 2


    lv = ptr_node->l;

    H = 1.0L / (1 << lv);
    one_over_2H = 1.0L/(2*H);

    //** >> Case head node **/
    if (lv == lmin)
    {
        int box_grid_idxNbr_1; // Box index on the right (+x)
        int box_grid_idxNbr_2; // Box index on the left (-x)
        int box_grid_idxNbr_3; // Box index behind (+y)
        int box_grid_idxNbr_4; // Box index in front  (-y)
        int box_grid_idxNbr_5; // Box index up (+z)
        int box_grid_idxNbr_6; // Box index down (-z)
        //** >> Cycle over interior grid points **/
        for (int i = 0; i < ptr_node->grid_intr_size; i++)
        {
            //** >> Box grid index **/
            box_grid_idx = ptr_node->ptr_grid_intr[i];

            //** >> Neighboring box grid indexes **/
            box_grid_idxNbr_1 = box_grid_idx + 1;
            box_grid_idxNbr_2 = box_grid_idx - 1;
            box_grid_idxNbr_3 = box_grid_idx + (ptr_node->box_real_dim_x + 1);
            box_grid_idxNbr_4 = box_grid_idx - (ptr_node->box_real_dim_x + 1);
            box_grid_idxNbr_5 = box_grid_idx + (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1);
            box_grid_idxNbr_6 = box_grid_idx - (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1);

            //** >> Filling the acceleration **/
            ptr_node->ptr_ax[box_grid_idx] = one_over_2H * (ptr_node->ptr_pot[box_grid_idxNbr_2] - ptr_node->ptr_pot[box_grid_idxNbr_1]);
            ptr_node->ptr_ay[box_grid_idx] = one_over_2H * (ptr_node->ptr_pot[box_grid_idxNbr_4] - ptr_node->ptr_pot[box_grid_idxNbr_3]);
            ptr_node->ptr_az[box_grid_idx] = one_over_2H * (ptr_node->ptr_pot[box_grid_idxNbr_6] - ptr_node->ptr_pot[box_grid_idxNbr_5]);
  
        }
    }
    //** >> Case branch nodes **/
    else
    {
        //** >> Border grid points **/
        for (int i = 0; i< ptr_node->grid_bder_size; i++)
        {
            box_grid_idx = ptr_node->ptr_grid_bder[i];
            box_idx_z = box_grid_idx / ((ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1));
            box_idx_y = (box_grid_idx - box_idx_z * (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1)) / (ptr_node->box_real_dim_x + 1);
            box_idx_x = box_grid_idx - box_idx_z * (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1) - box_idx_y * (ptr_node->box_real_dim_x + 1);
            cell_idx_x = box_idx_x + ptr_node->box_ts_x;
            cell_idx_y = box_idx_y + ptr_node->box_ts_y;
            cell_idx_z = box_idx_z + ptr_node->box_ts_z;

            
            pt_box_idx_1_x = (cell_idx_x >> 1) - ptr_node->ptr_pt->box_ts_x;
            pt_box_idx_2_x = ((cell_idx_x + 1) >> 1) - ptr_node->ptr_pt->box_ts_x;
            pt_box_idx_1_y = (cell_idx_y >> 1) - ptr_node->ptr_pt->box_ts_y;
            pt_box_idx_2_y = ((cell_idx_y + 1) >> 1) - ptr_node->ptr_pt->box_ts_y;
            pt_box_idx_1_z = (cell_idx_z >> 1) - ptr_node->ptr_pt->box_ts_z;
            pt_box_idx_2_z = ((cell_idx_z + 1) >> 1) - ptr_node->ptr_pt->box_ts_z;
            pt_box_grid_idx_1 = pt_box_idx_1_x + pt_box_idx_1_y * (ptr_node->ptr_pt->box_real_dim_x + 1) + pt_box_idx_1_z * (ptr_node->ptr_pt->box_real_dim_x + 1) * (ptr_node->ptr_pt->box_real_dim_y + 1);

            //** x-axis acceleration **/
            pt_box_grid_idx_2 = pt_box_idx_2_x + pt_box_idx_1_y * (ptr_node->ptr_pt->box_real_dim_x + 1) + pt_box_idx_1_z * (ptr_node->ptr_pt->box_real_dim_x + 1) * (ptr_node->ptr_pt->box_real_dim_y + 1);
            aux_acc = (ptr_node->ptr_pt->ptr_ax[pt_box_grid_idx_1] + ptr_node->ptr_pt->ptr_ax[pt_box_grid_idx_2]) * 0.5;
            ptr_node->ptr_ax[box_grid_idx] = aux_acc;

            //** y-axis acceleration **/
            pt_box_grid_idx_2 = pt_box_idx_1_x + pt_box_idx_2_y * (ptr_node->ptr_pt->box_real_dim_x + 1) + pt_box_idx_1_z * (ptr_node->ptr_pt->box_real_dim_x + 1) * (ptr_node->ptr_pt->box_real_dim_y + 1);
            aux_acc = (ptr_node->ptr_pt->ptr_ay[pt_box_grid_idx_1] + ptr_node->ptr_pt->ptr_ay[pt_box_grid_idx_2]) * 0.5;
            ptr_node->ptr_ay[box_grid_idx] = aux_acc;

            //** z-axis acceleration **/
            pt_box_grid_idx_2 = pt_box_idx_1_x + pt_box_idx_1_y * (ptr_node->ptr_pt->box_real_dim_x + 1) + pt_box_idx_2_z * (ptr_node->ptr_pt->box_real_dim_x + 1) * (ptr_node->ptr_pt->box_real_dim_y + 1);
            aux_acc = (ptr_node->ptr_pt->ptr_az[pt_box_grid_idx_1] + ptr_node->ptr_pt->ptr_az[pt_box_grid_idx_2]) * 0.5;
            ptr_node->ptr_az[box_grid_idx] = aux_acc;
        }

        //** >> Interior grid points **/
        for (int i = 0; i < ptr_node->grid_intr_size; i++)
        {
            box_grid_idx = ptr_node->ptr_grid_intr[i];
            ptr_node->ptr_ax[box_grid_idx] = one_over_2H * (ptr_node->ptr_pot[box_grid_idx - 1] - ptr_node->ptr_pot[box_grid_idx + 1]);
            ptr_node->ptr_ay[box_grid_idx] = one_over_2H * (ptr_node->ptr_pot[box_grid_idx - (ptr_node->box_real_dim_x + 1)] - ptr_node->ptr_pot[box_grid_idx + (ptr_node->box_real_dim_x + 1)]);
            ptr_node->ptr_az[box_grid_idx] = one_over_2H * (ptr_node->ptr_pot[box_grid_idx - (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1)] - ptr_node->ptr_pot[box_grid_idx + (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1)]);
        }
    }
}

    int
    grid_acceleration()
{
    //** >> Acceleration in the grid **/

    struct node *ptr_node;

    int no_pts;   // Number of parents in the cycle

    for (int lv = 0; lv < lmax - lmin + 1; lv++)
    {
        no_pts = GL_tentacles_size[lv];

        //** >> For cycle over parent nodes **/
        for (int i = 0; i < no_pts; i++)
        {
            ptr_node = GL_tentacles_old[lv][i];
            computing_grid_acceleration(ptr_node);
        }
        if (GL_tentacles_size[lv + 1] == 0)
        {
            lv = lmax - lmin + 1;
        }
    }


    return _SUCCESS_;
}