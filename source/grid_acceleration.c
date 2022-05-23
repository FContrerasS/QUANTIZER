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

static void computing_grid_acceleration_head(struct node *ptr_head)
{
    int lv = ptr_head->lv;

    int box_grid_idx; // Box grid index

    vtype H = 1.0L / (1 << lv);
    vtype one_over_2H = 1.0L / (2 * H);

    //** >> Case head node **/

    int box_grid_idxNbr_1; // Box index on the right (+x)
    int box_grid_idxNbr_2; // Box index on the left (-x)
    int box_grid_idxNbr_3; // Box index behind (+y)
    int box_grid_idxNbr_4; // Box index in front  (-y)
    int box_grid_idxNbr_5; // Box index up (+z)
    int box_grid_idxNbr_6; // Box index down (-z)

    int grid_box_real_dim_X_head = (ptr_head->box_real_dim_x + 1);
    int grid_box_real_dim_X_times_Y_head = (ptr_head->box_real_dim_x + 1) * (ptr_head->box_real_dim_y + 1);

    //** >> Cycle over interior grid points **/
    for (int i = 0; i < ptr_head->grid_intr_size; i++)
    {
        //** >> Box grid index **/
        box_grid_idx = ptr_head->ptr_intr_grid_idx[i];

        //** >> Neighboring box grid indexes **/
        box_grid_idxNbr_1 = box_grid_idx + 1;
        box_grid_idxNbr_2 = box_grid_idx - 1;
        box_grid_idxNbr_3 = box_grid_idx + grid_box_real_dim_X_head;
        box_grid_idxNbr_4 = box_grid_idx - grid_box_real_dim_X_head;
        box_grid_idxNbr_5 = box_grid_idx + grid_box_real_dim_X_times_Y_head;
        box_grid_idxNbr_6 = box_grid_idx - grid_box_real_dim_X_times_Y_head;

        //** >> Filling the acceleration **/
        ptr_head->ptr_ax[box_grid_idx] = one_over_2H * (ptr_head->ptr_pot[box_grid_idxNbr_2] - ptr_head->ptr_pot[box_grid_idxNbr_1]);
        ptr_head->ptr_ay[box_grid_idx] = one_over_2H * (ptr_head->ptr_pot[box_grid_idxNbr_4] - ptr_head->ptr_pot[box_grid_idxNbr_3]);
        ptr_head->ptr_az[box_grid_idx] = one_over_2H * (ptr_head->ptr_pot[box_grid_idxNbr_6] - ptr_head->ptr_pot[box_grid_idxNbr_5]);
    }
}

static void computing_grid_acceleration_branch(struct node *ptr_node)
{
    vtype aux_acc; // Auxiliary value for the acceleration

    // int cell_idx_x; // The cell index
    // int cell_idx_y; // The cell index
    // int cell_idx_z; // The cell index

    // int box_idx_x;    // Box grid index in X direcction
    // int box_idx_y;    // Box grid index in Y direcction
    // int box_idx_z;    // Box grid index in Z direcction
    int box_grid_idx; // Box grid index

    int box_idx_1_x_node;    // Parent box index 1 at X direcction
    int box_idx_1_y_node;    // Parent box index 1 at Y direcction
    int box_idx_1_z_node;    // Parent box index 1 at Z direcction
    int box_idx_2_x_node;    // Parent box index 2 at X direcction
    int box_idx_2_y_node;    // Parent box index 2 at Y direcction
    int box_idx_2_z_node;    // Parent box index 2 at Z direcction
    int box_grid_idx_1_node; // Parent box grid index 1
    int box_grid_idx_2_node; // Parent box grid index 2

    int lv = ptr_node->lv;

    vtype H = 1.0L / (1 << lv);
    vtype one_over_2H = 1.0L / (2 * H);

    int grid_box_real_dim_X_node = (ptr_node->box_real_dim_x + 1);
    int grid_box_real_dim_X_times_Y_node = (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1);

    int grid_box_real_dim_X_pt = (ptr_node->ptr_pt->box_real_dim_x + 1);
    int grid_box_real_dim_X_times_Y_pt = (ptr_node->ptr_pt->box_real_dim_x + 1) * (ptr_node->ptr_pt->box_real_dim_y + 1);

    //** >> Border grid points **/
    for (int i = 0; i < ptr_node->grid_bder_size; i++)
    {
        box_grid_idx = ptr_node->ptr_bder_grid_idx[i];
        // box_idx_z = box_grid_idx / grid_box_real_dim_X_times_Y_node;
        // box_idx_y = (box_grid_idx - box_idx_z * grid_box_real_dim_X_times_Y_node) / grid_box_real_dim_X_node;
        // box_idx_x = box_grid_idx - box_idx_z * grid_box_real_dim_X_times_Y_node - box_idx_y * grid_box_real_dim_X_node;
        // cell_idx_x = box_idx_x + ptr_node->box_ts_x;
        // cell_idx_y = box_idx_y + ptr_node->box_ts_y;
        // cell_idx_z = box_idx_z + ptr_node->box_ts_z;

        // cell_idx_x = ptr_node->ptr_bder_grid_cell_idx_x[i];
        // cell_idx_y = ptr_node->ptr_bder_grid_cell_idx_y[i];
        // cell_idx_z = ptr_node->ptr_bder_grid_cell_idx_z[i];

        box_idx_1_x_node = (ptr_node->ptr_bder_grid_cell_idx_x[i] >> 1) - ptr_node->ptr_pt->box_ts_x;
        box_idx_2_x_node = ((ptr_node->ptr_bder_grid_cell_idx_x[i] + 1) >> 1) - ptr_node->ptr_pt->box_ts_x;
        box_idx_1_y_node = (ptr_node->ptr_bder_grid_cell_idx_y[i] >> 1) - ptr_node->ptr_pt->box_ts_y;
        box_idx_2_y_node = ((ptr_node->ptr_bder_grid_cell_idx_y[i] + 1) >> 1) - ptr_node->ptr_pt->box_ts_y;
        box_idx_1_z_node = (ptr_node->ptr_bder_grid_cell_idx_z[i] >> 1) - ptr_node->ptr_pt->box_ts_z;
        box_idx_2_z_node = ((ptr_node->ptr_bder_grid_cell_idx_z[i] + 1) >> 1) - ptr_node->ptr_pt->box_ts_z;

        box_grid_idx_1_node = box_idx_1_x_node + box_idx_1_y_node * grid_box_real_dim_X_pt + box_idx_1_z_node * grid_box_real_dim_X_times_Y_pt;

        //** x-axis acceleration **/
        box_grid_idx_2_node = box_idx_2_x_node + box_idx_1_y_node * grid_box_real_dim_X_pt + box_idx_1_z_node * grid_box_real_dim_X_times_Y_pt;
        aux_acc = (ptr_node->ptr_pt->ptr_ax[box_grid_idx_1_node] + ptr_node->ptr_pt->ptr_ax[box_grid_idx_2_node]) * 0.5;
        ptr_node->ptr_ax[box_grid_idx] = aux_acc;

        //** y-axis acceleration **/
        box_grid_idx_2_node = box_idx_1_x_node + box_idx_2_y_node * grid_box_real_dim_X_pt + box_idx_1_z_node * grid_box_real_dim_X_times_Y_pt;
        aux_acc = (ptr_node->ptr_pt->ptr_ay[box_grid_idx_1_node] + ptr_node->ptr_pt->ptr_ay[box_grid_idx_2_node]) * 0.5;
        ptr_node->ptr_ay[box_grid_idx] = aux_acc;

        //** z-axis acceleration **/
        box_grid_idx_2_node = box_idx_1_x_node + box_idx_1_y_node * grid_box_real_dim_X_pt + box_idx_2_z_node * grid_box_real_dim_X_times_Y_pt;
        aux_acc = (ptr_node->ptr_pt->ptr_az[box_grid_idx_1_node] + ptr_node->ptr_pt->ptr_az[box_grid_idx_2_node]) * 0.5;
        ptr_node->ptr_az[box_grid_idx] = aux_acc;
    }

    //** >> Interior grid points **/
    for (int i = 0; i < ptr_node->grid_intr_size; i++)
    {
        box_grid_idx = ptr_node->ptr_intr_grid_idx[i];
        ptr_node->ptr_ax[box_grid_idx] = one_over_2H * (ptr_node->ptr_pot[box_grid_idx - 1] - ptr_node->ptr_pot[box_grid_idx + 1]);
        ptr_node->ptr_ay[box_grid_idx] = one_over_2H * (ptr_node->ptr_pot[box_grid_idx - grid_box_real_dim_X_node] - ptr_node->ptr_pot[box_grid_idx + grid_box_real_dim_X_node]);
        ptr_node->ptr_az[box_grid_idx] = one_over_2H * (ptr_node->ptr_pot[box_grid_idx - grid_box_real_dim_X_times_Y_node] - ptr_node->ptr_pot[box_grid_idx + grid_box_real_dim_X_times_Y_node]);
    }
}

int grid_acceleration()
{
    //** >> Acceleration in the grid **/

    // Head
    computing_grid_acceleration_head(GL_tentacles[0][0]);

    // Branches
    for (int lv = 1; lv < GL_tentacles_level_max + 1; lv++)
    {
        // Number of parent of the level = GL_tentacles_size[lv];

        //** >> For cycle over parent nodes **/
        for (int i = 0; i < GL_tentacles_size[lv]; i++)
        {
            // ptr_node = GL_tentacles[lv][i];
            computing_grid_acceleration_branch(GL_tentacles[lv][i]);
        }
    }

    return _SUCCESS_;
}