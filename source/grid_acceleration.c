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

//** >> Local Functions
static void computing_grid_acceleration_head(struct node *ptr_head);
static void computing_grid_acceleration_head2(struct node *ptr_head);
static void computing_simulation_boundary_grid_point_acceleration(struct node *ptr_node);
static void computing_grid_acceleration_branch(struct node *ptr_node);
static void computing_grid_acceleration_branch2(struct node *ptr_node);

static void computing_grid_acceleration_head(struct node *ptr_head)
{
    int lv = ptr_head->lv;

    int box_grid_idx; // Box grid index

    vtype H = 1.0L / (1 << lv);
    vtype one_over_2H = 1.0L / (2.0 * H);

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
        box_grid_idx = ptr_head->ptr_intr_box_grid_idx[i];

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

static void computing_grid_acceleration_head2(struct node *ptr_head)
{
    int lv = ptr_head->lv;

    int box_grid_idx; // Box grid index

    vtype H = 1.0L / (1 << lv);
    vtype one_over_2H = 1.0L / (2.0 * H);
    vtype one_over_12H = 1.0L / (12.0 * H);
    vtype two_over_3H = 2.0L / (3.0 * H);

    //** >> Case head node **/

    int box_grid_idxNbr_1; // Box index on the right (+x)
    int box_grid_idxNbr_2; // Box index on the left (-x)
    int box_grid_idxNbr_3; // Box index behind (+y)
    int box_grid_idxNbr_4; // Box index in front  (-y)
    int box_grid_idxNbr_5; // Box index up (+z)
    int box_grid_idxNbr_6; // Box index down (-z)

    int box_grid_idxNbr_7; // Box index down (+2x)
    int box_grid_idxNbr_8; // Box index down (-2x)
    int box_grid_idxNbr_9; // Box index down (+2y)
    int box_grid_idxNbr_10; // Box index down (-2y)
    int box_grid_idxNbr_11; // Box index down (+2z)
    int box_grid_idxNbr_12; // Box index down (-2z)

    int box_idx_x;
    int box_idx_y;
    int box_idx_z;
    int box_idx;

    int grid_box_real_dim_X_head = (ptr_head->box_real_dim_x + 1);
    int grid_box_real_dim_X_times_Y_head = (ptr_head->box_real_dim_x + 1) * (ptr_head->box_real_dim_y + 1);

    int box_real_dim_X_head = ptr_head->box_real_dim_x;
    int box_real_dim_X_times_Y_head = ptr_head->box_real_dim_x * ptr_head->box_real_dim_y;

    //** >> Cycle over interior grid points **/
    for (int i = 0; i < ptr_head->grid_intr_size; i++)
    {
        //** >> Box grid index **/
        box_grid_idx = ptr_head->ptr_intr_box_grid_idx[i];

        box_idx_x = ptr_head->ptr_intr_grid_cell_idx_x[i] - ptr_head->box_ts_x;
        box_idx_y = ptr_head->ptr_intr_grid_cell_idx_y[i] - ptr_head->box_ts_y;
        box_idx_z = ptr_head->ptr_intr_grid_cell_idx_z[i] - ptr_head->box_ts_z;

        box_idx = box_idx_x + box_idx_y * ptr_head->box_real_dim_x + box_idx_z * ptr_head->box_real_dim_x * ptr_head->box_real_dim_y;

        //** >> Neighboring box grid indexes **/
        box_grid_idxNbr_1 = box_grid_idx + 1;
        box_grid_idxNbr_2 = box_grid_idx - 1;
        box_grid_idxNbr_3 = box_grid_idx + grid_box_real_dim_X_head;
        box_grid_idxNbr_4 = box_grid_idx - grid_box_real_dim_X_head;
        box_grid_idxNbr_5 = box_grid_idx + grid_box_real_dim_X_times_Y_head;
        box_grid_idxNbr_6 = box_grid_idx - grid_box_real_dim_X_times_Y_head;

        //** >> Filling the acceleration **/
        //x
        if (ptr_head->ptr_box[box_idx - 2] > -4 && ptr_head->ptr_box[box_idx + 1] > -4)
        {
            box_grid_idxNbr_7 = box_grid_idx + 2;
            box_grid_idxNbr_8 = box_grid_idx - 2;
            ptr_head->ptr_ax[box_grid_idx] = (ptr_head->ptr_pot[box_grid_idxNbr_2] - ptr_head->ptr_pot[box_grid_idxNbr_1]) * two_over_3H - (ptr_head->ptr_pot[box_grid_idxNbr_8] - ptr_head->ptr_pot[box_grid_idxNbr_7]) * one_over_12H;
        }
        else if(boundary_type == 0)
        {
            box_grid_idxNbr_7 = box_grid_idx + 2;
            box_grid_idxNbr_8 = box_grid_idx - 2;
            box_grid_idxNbr_7 -= ptr_head->ptr_box[box_idx - 2] == -6 ? 0 : (1 << lv);
            box_grid_idxNbr_8 += ptr_head->ptr_box[box_idx - 2] == -6 ? (1 << lv) : 0;
            ptr_head->ptr_ax[box_grid_idx] = (ptr_head->ptr_pot[box_grid_idxNbr_2] - ptr_head->ptr_pot[box_grid_idxNbr_1]) * two_over_3H - (ptr_head->ptr_pot[box_grid_idxNbr_8] - ptr_head->ptr_pot[box_grid_idxNbr_7]) * one_over_12H;
        }
        else
        {
            ptr_head->ptr_ax[box_grid_idx] = one_over_2H * (ptr_head->ptr_pot[box_grid_idxNbr_2] - ptr_head->ptr_pot[box_grid_idxNbr_1]);
        }

        //y
        if (ptr_head->ptr_box[box_idx - 2 * box_real_dim_X_head] > -4 && ptr_head->ptr_box[box_idx + box_real_dim_X_head] > -4)
        {
            box_grid_idxNbr_9  = box_grid_idx + 2 * grid_box_real_dim_X_head;
            box_grid_idxNbr_10 = box_grid_idx - 2 * grid_box_real_dim_X_head;
            ptr_head->ptr_ay[box_grid_idx] = (ptr_head->ptr_pot[box_grid_idxNbr_4] - ptr_head->ptr_pot[box_grid_idxNbr_3]) * two_over_3H - (ptr_head->ptr_pot[box_grid_idxNbr_10] - ptr_head->ptr_pot[box_grid_idxNbr_9]) * one_over_12H;
        }
        else if (boundary_type == 0)
        {
            box_grid_idxNbr_9  = box_grid_idx + 2 * grid_box_real_dim_X_head;
            box_grid_idxNbr_10 = box_grid_idx - 2 * grid_box_real_dim_X_head;
            box_grid_idxNbr_9  -= ptr_head->ptr_box[box_idx - 2 * box_real_dim_X_head] == -6 ? 0 : (1 << lv) * grid_box_real_dim_X_head;
            box_grid_idxNbr_10 += ptr_head->ptr_box[box_idx - 2 * box_real_dim_X_head] == -6 ? (1 << lv) * grid_box_real_dim_X_head : 0;
            ptr_head->ptr_ay[box_grid_idx] = (ptr_head->ptr_pot[box_grid_idxNbr_4] - ptr_head->ptr_pot[box_grid_idxNbr_3]) * two_over_3H - (ptr_head->ptr_pot[box_grid_idxNbr_10] - ptr_head->ptr_pot[box_grid_idxNbr_9]) * one_over_12H;
        }
        else
        {
            ptr_head->ptr_ay[box_grid_idx] = one_over_2H * (ptr_head->ptr_pot[box_grid_idxNbr_4] - ptr_head->ptr_pot[box_grid_idxNbr_3]);
        }

        //z
        if (ptr_head->ptr_box[box_idx - 2 * box_real_dim_X_times_Y_head] > -4 && ptr_head->ptr_box[box_idx + box_real_dim_X_times_Y_head] > -4)
        {
            box_grid_idxNbr_11 = box_grid_idx + 2 * grid_box_real_dim_X_times_Y_head;
            box_grid_idxNbr_12 = box_grid_idx - 2 * grid_box_real_dim_X_times_Y_head;
            ptr_head->ptr_az[box_grid_idx] = (ptr_head->ptr_pot[box_grid_idxNbr_6] - ptr_head->ptr_pot[box_grid_idxNbr_5]) * two_over_3H - (ptr_head->ptr_pot[box_grid_idxNbr_12] - ptr_head->ptr_pot[box_grid_idxNbr_11]) * one_over_12H;
        }
        else if (boundary_type == 0)
        {
            box_grid_idxNbr_11 = box_grid_idx + 2 * grid_box_real_dim_X_times_Y_head;
            box_grid_idxNbr_12 = box_grid_idx - 2 * grid_box_real_dim_X_times_Y_head;
            box_grid_idxNbr_11 -= ptr_head->ptr_box[box_idx - 2 * box_real_dim_X_times_Y_head] == -6 ? 0 : (1 << lv) * grid_box_real_dim_X_times_Y_head;
            box_grid_idxNbr_12 += ptr_head->ptr_box[box_idx - 2 * box_real_dim_X_times_Y_head] == -6 ? (1 << lv) * grid_box_real_dim_X_times_Y_head : 0;
            ptr_head->ptr_az[box_grid_idx] = (ptr_head->ptr_pot[box_grid_idxNbr_6] - ptr_head->ptr_pot[box_grid_idxNbr_5]) * two_over_3H - (ptr_head->ptr_pot[box_grid_idxNbr_12] - ptr_head->ptr_pot[box_grid_idxNbr_11]) * one_over_12H;
        }
        else
        {
            ptr_head->ptr_az[box_grid_idx] = one_over_2H * (ptr_head->ptr_pot[box_grid_idxNbr_6] - ptr_head->ptr_pot[box_grid_idxNbr_5]);
        }
    }
}

static void computing_simulation_boundary_grid_point_acceleration(struct node *ptr_node)
{
    int lv = ptr_node->lv;

    int box_grid_idx; // Box grid index

    int aux_i;
    int aux_j;
    int aux_k;

    vtype H = 1.0L / (1 << lv);

    vtype dist;

    vtype aux_coeff_1 = -_G_ * GL_total_mass_initial;
    vtype aux_coeff_2;

    //** >> Simulation Boundary grid points **/
    for (int i = 0; i < ptr_node->grid_SIMULATION_BOUNDARY_size; i++)
    {
        box_grid_idx = ptr_node->ptr_SIMULATION_BOUNDARY_box_grid_idx[i];

        if (ptr_node->ptr_SIMULATION_BOUNDARY_grid_cell_idx_x[i] == 0 ||
            ptr_node->ptr_SIMULATION_BOUNDARY_grid_cell_idx_x[i] == (1 << lv))
        {
            aux_i = GL_cm[0] < 0.5 ? 0 : (1 << lv) * H;
        }
        else
        {
            aux_i = ptr_node->ptr_SIMULATION_BOUNDARY_grid_cell_idx_x[i] * H;
        }

        if (ptr_node->ptr_SIMULATION_BOUNDARY_grid_cell_idx_y[i] == 0 ||
            ptr_node->ptr_SIMULATION_BOUNDARY_grid_cell_idx_y[i] == (1 << lv))
        {
            aux_j = GL_cm[1] < 0.5 ? 0 : (1 << lv) * H;
        }
        else
        {
            aux_j = ptr_node->ptr_SIMULATION_BOUNDARY_grid_cell_idx_y[i];
        }

        if (ptr_node->ptr_SIMULATION_BOUNDARY_grid_cell_idx_z[i] == 0 ||
            ptr_node->ptr_SIMULATION_BOUNDARY_grid_cell_idx_z[i] == (1 << lv))
        {
            aux_k = GL_cm[2] < 0.5 ? 0 : (1 << lv) * H;
        }
        else
        {
            aux_k = ptr_node->ptr_SIMULATION_BOUNDARY_grid_cell_idx_z[i] * H;
        }

        dist = (aux_i - GL_cm[0]) * (aux_i - GL_cm[0]) + (aux_j - GL_cm[1]) * (aux_j - GL_cm[1]) + (aux_k - GL_cm[2]) * (aux_k - GL_cm[2]);
        dist = sqrt(dist);

        aux_coeff_2 = aux_coeff_1 / (dist * dist * dist);

        ptr_node->ptr_ax[box_grid_idx] = aux_coeff_2 * (aux_i - GL_cm[0]);
        ptr_node->ptr_ay[box_grid_idx] = aux_coeff_2 * (aux_j - GL_cm[1]);
        ptr_node->ptr_az[box_grid_idx] = aux_coeff_2 * (aux_k - GL_cm[2]);
    }
}

static void computing_grid_acceleration_branch(struct node *ptr_node)
{
    //vtype aux_acc; // Auxiliary value for the acceleration

    struct node *ptr_pt = ptr_node->ptr_pt;
    

    int box_grid_idx; // Box grid index

    int box_idx_1_x_pt;    // Parent box index 1 at X direcction
    int box_idx_1_y_pt;    // Parent box index 1 at Y direcction
    int box_idx_1_z_pt;    // Parent box index 1 at Z direcction
    int box_idx_2_x_pt;    // Parent box index 2 at X direcction
    int box_idx_2_y_pt;    // Parent box index 2 at Y direcction
    int box_idx_2_z_pt;    // Parent box index 2 at Z direcction
    int box_grid_idx_1_pt; // Parent box grid index 1
    int box_grid_idx_2_pt; // Parent box grid index 2
    int box_grid_idx_3_pt; // Parent box grid index 3
    int box_grid_idx_4_pt; // Parent box grid index 4
    int box_grid_idx_5_pt; // Parent box grid index 5
    int box_grid_idx_6_pt; // Parent box grid index 6
    int box_grid_idx_7_pt; // Parent box grid index 7
    int box_grid_idx_8_pt; // Parent box grid index 8

    int lv = ptr_node->lv;

    vtype H = 1.0L / (1 << lv);
    vtype one_over_2H = 1.0L / (2.0 * H);

    int aux_idx_1_x;
    int aux_idx_1_y;
    int aux_idx_1_z;
    int aux_idx_2_x;
    int aux_idx_2_y;
    int aux_idx_2_z;


    int grid_box_real_dim_X_node = (ptr_node->box_real_dim_x + 1);
    int grid_box_real_dim_X_times_Y_node = (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1);

    int grid_box_real_dim_X_pt = (ptr_pt->box_real_dim_x + 1);
    int grid_box_real_dim_X_times_Y_pt = (ptr_pt->box_real_dim_x + 1) * (ptr_pt->box_real_dim_y + 1);

    //** >> Simulation Boundary grid points **/
    computing_simulation_boundary_grid_point_acceleration(ptr_node);

    //** >> Border grid points **/
    for (int i = 0; i < ptr_node->grid_bder_size; i++)
    {
        //Notes here that a middle grid child grid point will never be a border grid point because the child cells are
        // composed by 8 grid cells where that grid point is in the middle.
        box_grid_idx = ptr_node->ptr_bder_box_grid_idx[i];

        aux_idx_1_x = ptr_node->ptr_bder_grid_cell_idx_x[i] >> 1;
        aux_idx_1_y = ptr_node->ptr_bder_grid_cell_idx_y[i] >> 1;
        aux_idx_1_z = ptr_node->ptr_bder_grid_cell_idx_z[i] >> 1;
        aux_idx_2_x = (ptr_node->ptr_bder_grid_cell_idx_x[i] + 1) >> 1;
        aux_idx_2_y = (ptr_node->ptr_bder_grid_cell_idx_y[i] + 1) >> 1;
        aux_idx_2_z = (ptr_node->ptr_bder_grid_cell_idx_z[i] + 1) >> 1;

        box_idx_1_x_pt = aux_idx_1_x - ptr_pt->box_ts_x;
        box_idx_1_y_pt = aux_idx_1_y - ptr_pt->box_ts_y;
        box_idx_1_z_pt = aux_idx_1_z - ptr_pt->box_ts_z;
        box_idx_2_x_pt = aux_idx_2_x - ptr_pt->box_ts_x;
        box_idx_2_y_pt = aux_idx_2_y - ptr_pt->box_ts_y;
        box_idx_2_z_pt = aux_idx_2_z - ptr_pt->box_ts_z;

        if (ptr_pt->pbc_crosses_the_boundary_simulation_box == true)
        {
            if (aux_idx_1_x > ptr_pt->box_max_x)
            {
                box_idx_1_x_pt -= (1 << (lv - 1));
                box_idx_2_x_pt -= (1 << (lv - 1));
            }

            if (aux_idx_1_y > ptr_pt->box_max_y)
            {
                box_idx_1_y_pt -= (1 << (lv - 1));
                box_idx_2_y_pt -= (1 << (lv - 1));
            }

            if (aux_idx_1_z > ptr_pt->box_max_z)
            {
                box_idx_1_z_pt -= (1 << (lv - 1));
                box_idx_2_z_pt -= (1 << (lv - 1));
            }
        }

        // box_idx_1_x_node = (ptr_node->ptr_bder_grid_cell_idx_x[i] >> 1) - ptr_pt->box_ts_x;
        // box_idx_2_x_node = ((ptr_node->ptr_bder_grid_cell_idx_x[i] + 1) >> 1) - ptr_pt->box_ts_x;
        // box_idx_1_y_node = (ptr_node->ptr_bder_grid_cell_idx_y[i] >> 1) - ptr_pt->box_ts_y;
        // box_idx_2_y_node = ((ptr_node->ptr_bder_grid_cell_idx_y[i] + 1) >> 1) - ptr_pt->box_ts_y;
        // box_idx_1_z_node = (ptr_node->ptr_bder_grid_cell_idx_z[i] >> 1) - ptr_pt->box_ts_z;
        // box_idx_2_z_node = ((ptr_node->ptr_bder_grid_cell_idx_z[i] + 1) >> 1) - ptr_pt->box_ts_z;

        box_grid_idx_1_pt = box_idx_1_x_pt + box_idx_1_y_pt * grid_box_real_dim_X_pt + box_idx_1_z_pt * grid_box_real_dim_X_times_Y_pt;
        box_grid_idx_2_pt = box_idx_2_x_pt + box_idx_1_y_pt * grid_box_real_dim_X_pt + box_idx_1_z_pt * grid_box_real_dim_X_times_Y_pt;
        box_grid_idx_3_pt = box_idx_1_x_pt + box_idx_2_y_pt * grid_box_real_dim_X_pt + box_idx_1_z_pt * grid_box_real_dim_X_times_Y_pt;
        box_grid_idx_4_pt = box_idx_2_x_pt + box_idx_2_y_pt * grid_box_real_dim_X_pt + box_idx_1_z_pt * grid_box_real_dim_X_times_Y_pt;
        box_grid_idx_5_pt = box_idx_1_x_pt + box_idx_1_y_pt * grid_box_real_dim_X_pt + box_idx_2_z_pt * grid_box_real_dim_X_times_Y_pt;
        box_grid_idx_6_pt = box_idx_2_x_pt + box_idx_1_y_pt * grid_box_real_dim_X_pt + box_idx_2_z_pt * grid_box_real_dim_X_times_Y_pt;
        box_grid_idx_7_pt = box_idx_1_x_pt + box_idx_2_y_pt * grid_box_real_dim_X_pt + box_idx_2_z_pt * grid_box_real_dim_X_times_Y_pt;
        box_grid_idx_8_pt = box_idx_2_x_pt + box_idx_2_y_pt * grid_box_real_dim_X_pt + box_idx_2_z_pt * grid_box_real_dim_X_times_Y_pt;

        //** x-axis acceleration **/
        // box_grid_idx_2_pt = box_idx_2_x_pt + box_idx_1_y_pt * grid_box_real_dim_X_pt + box_idx_1_z_pt * grid_box_real_dim_X_times_Y_pt;
        // aux_acc = (ptr_pt->ptr_ax[box_grid_idx_1_pt] + ptr_pt->ptr_ax[box_grid_idx_2_pt]) * 0.5;
        // ptr_node->ptr_ax[box_grid_idx] = aux_acc;

        // //** y-axis acceleration **/
        // box_grid_idx_2_pt = box_idx_1_x_pt + box_idx_2_y_pt * grid_box_real_dim_X_pt + box_idx_1_z_pt * grid_box_real_dim_X_times_Y_pt;
        // aux_acc = (ptr_pt->ptr_ay[box_grid_idx_1_pt] + ptr_pt->ptr_ay[box_grid_idx_2_pt]) * 0.5;
        // ptr_node->ptr_ay[box_grid_idx] = aux_acc;

        // //** z-axis acceleration **/
        // box_grid_idx_2_pt = box_idx_1_x_pt + box_idx_1_y_pt * grid_box_real_dim_X_pt + box_idx_2_z_pt * grid_box_real_dim_X_times_Y_pt;
        // aux_acc = (ptr_pt->ptr_az[box_grid_idx_1_pt] + ptr_pt->ptr_az[box_grid_idx_2_pt]) * 0.5;
        // ptr_node->ptr_az[box_grid_idx] = aux_acc;

        //** x-axis acceleration **/
        ptr_node->ptr_ax[box_grid_idx] = 0.125 * (ptr_pt->ptr_ax[box_grid_idx_1_pt] + ptr_pt->ptr_ax[box_grid_idx_2_pt] +
                                                  ptr_pt->ptr_ax[box_grid_idx_3_pt] + ptr_pt->ptr_ax[box_grid_idx_4_pt] +
                                                  ptr_pt->ptr_ax[box_grid_idx_5_pt] + ptr_pt->ptr_ax[box_grid_idx_6_pt] +
                                                  ptr_pt->ptr_ax[box_grid_idx_7_pt] + ptr_pt->ptr_ax[box_grid_idx_8_pt]);

        //** y-axis acceleration **/
        ptr_node->ptr_ay[box_grid_idx] = 0.125 * (ptr_pt->ptr_ay[box_grid_idx_1_pt] + ptr_pt->ptr_ay[box_grid_idx_2_pt] +
                                                  ptr_pt->ptr_ay[box_grid_idx_3_pt] + ptr_pt->ptr_ay[box_grid_idx_4_pt] +
                                                  ptr_pt->ptr_ay[box_grid_idx_5_pt] + ptr_pt->ptr_ay[box_grid_idx_6_pt] +
                                                  ptr_pt->ptr_ay[box_grid_idx_7_pt] + ptr_pt->ptr_ay[box_grid_idx_8_pt]);

        //** z-axis acceleration **/
        ptr_node->ptr_az[box_grid_idx] = 0.125 * (ptr_pt->ptr_az[box_grid_idx_1_pt] + ptr_pt->ptr_az[box_grid_idx_2_pt] +
                                                  ptr_pt->ptr_az[box_grid_idx_3_pt] + ptr_pt->ptr_az[box_grid_idx_4_pt] +
                                                  ptr_pt->ptr_az[box_grid_idx_5_pt] + ptr_pt->ptr_az[box_grid_idx_6_pt] +
                                                  ptr_pt->ptr_az[box_grid_idx_7_pt] + ptr_pt->ptr_az[box_grid_idx_8_pt]);
        
    }

    //** >> Interior grid points **/
    for (int i = 0; i < ptr_node->grid_intr_size; i++)
    {
        box_grid_idx = ptr_node->ptr_intr_box_grid_idx[i];
        ptr_node->ptr_ax[box_grid_idx] = one_over_2H * (ptr_node->ptr_pot[box_grid_idx - 1] - ptr_node->ptr_pot[box_grid_idx + 1]);
        ptr_node->ptr_ay[box_grid_idx] = one_over_2H * (ptr_node->ptr_pot[box_grid_idx - grid_box_real_dim_X_node] - ptr_node->ptr_pot[box_grid_idx + grid_box_real_dim_X_node]);
        ptr_node->ptr_az[box_grid_idx] = one_over_2H * (ptr_node->ptr_pot[box_grid_idx - grid_box_real_dim_X_times_Y_node] - ptr_node->ptr_pot[box_grid_idx + grid_box_real_dim_X_times_Y_node]);
    }
}

static void computing_grid_acceleration_branch2(struct node *ptr_node)
{
    //vtype aux_acc; // Auxiliary value for the acceleration

    struct node *ptr_pt = ptr_node->ptr_pt;

    int box_grid_idx; // Box grid index

    int box_idx_x_node;
    int box_idx_y_node;
    int box_idx_z_node;
    int box_idx_node;

    int box_grid_idxNbr_1; // Box index on the right (+x)
    int box_grid_idxNbr_2; // Box index on the left (-x)
    int box_grid_idxNbr_3; // Box index behind (+y)
    int box_grid_idxNbr_4; // Box index in front  (-y)
    int box_grid_idxNbr_5; // Box index up (+z)
    int box_grid_idxNbr_6; // Box index down (-z)

    int box_grid_idxNbr_7;  // Box index down (+2x)
    int box_grid_idxNbr_8;  // Box index down (-2x)
    int box_grid_idxNbr_9;  // Box index down (+2y)
    int box_grid_idxNbr_10; // Box index down (-2y)
    int box_grid_idxNbr_11; // Box index down (+2z)
    int box_grid_idxNbr_12; // Box index down (-2z)

    int box_idx_1_x_pt;    // Parent box index 1 at X direcction
    int box_idx_1_y_pt;    // Parent box index 1 at Y direcction
    int box_idx_1_z_pt;    // Parent box index 1 at Z direcction
    int box_idx_2_x_pt;    // Parent box index 2 at X direcction
    int box_idx_2_y_pt;    // Parent box index 2 at Y direcction
    int box_idx_2_z_pt;    // Parent box index 2 at Z direcction
    int box_grid_idx_1_pt; // Parent box grid index 1
    int box_grid_idx_2_pt; // Parent box grid index 2
    int box_grid_idx_3_pt; // Parent box grid index 3
    int box_grid_idx_4_pt; // Parent box grid index 4
    int box_grid_idx_5_pt; // Parent box grid index 5
    int box_grid_idx_6_pt; // Parent box grid index 6
    int box_grid_idx_7_pt; // Parent box grid index 7
    int box_grid_idx_8_pt; // Parent box grid index 8

    int aux_idx_1_x;
    int aux_idx_1_y;
    int aux_idx_1_z;
    int aux_idx_2_x;
    int aux_idx_2_y;
    int aux_idx_2_z;

    int aux_idx_x;
    int aux_idx_y;
    int aux_idx_z;

    int lv = ptr_node->lv;

    vtype H = 1.0L / (1 << lv);
    vtype one_over_2H = 1.0L / (2.0 * H);
    vtype one_over_12H = 1.0L / (12.0 * H);
    vtype two_over_3H = 2.0L / (3.0 * H);

    int grid_box_real_dim_X_node = (ptr_node->box_real_dim_x + 1);
    int grid_box_real_dim_X_times_Y_node = (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1);

    int grid_box_real_dim_X_pt = (ptr_pt->box_real_dim_x + 1);
    int grid_box_real_dim_X_times_Y_pt = (ptr_pt->box_real_dim_x + 1) * (ptr_pt->box_real_dim_y + 1);

    int box_real_dim_X_node = ptr_node->box_real_dim_x;
    int box_real_dim_X_times_Y_node = ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;

    bool check;

    //** >> Simulation Boundary grid points **/
    computing_simulation_boundary_grid_point_acceleration(ptr_node);

    //** >> Border grid points **/
    for (int i = 0; i < ptr_node->grid_bder_size; i++)
    {
        box_grid_idx = ptr_node->ptr_bder_box_grid_idx[i];

        aux_idx_1_x = ptr_node->ptr_bder_grid_cell_idx_x[i] >> 1;
        aux_idx_1_y = ptr_node->ptr_bder_grid_cell_idx_y[i] >> 1;
        aux_idx_1_z = ptr_node->ptr_bder_grid_cell_idx_z[i] >> 1;
        aux_idx_2_x = (ptr_node->ptr_bder_grid_cell_idx_x[i] + 1) >> 1;
        aux_idx_2_y = (ptr_node->ptr_bder_grid_cell_idx_y[i] + 1) >> 1;
        aux_idx_2_z = (ptr_node->ptr_bder_grid_cell_idx_z[i] + 1) >> 1;

        box_idx_1_x_pt = aux_idx_1_x - ptr_pt->box_ts_x;
        box_idx_1_y_pt = aux_idx_1_y - ptr_pt->box_ts_y;
        box_idx_1_z_pt = aux_idx_1_z - ptr_pt->box_ts_z;
        box_idx_2_x_pt = aux_idx_2_x - ptr_pt->box_ts_x;
        box_idx_2_y_pt = aux_idx_2_y - ptr_pt->box_ts_y;
        box_idx_2_z_pt = aux_idx_2_z - ptr_pt->box_ts_z;

        if (ptr_pt->pbc_crosses_the_boundary_simulation_box == true)
        {
            if (aux_idx_1_x > ptr_node->box_max_x)
            {
                box_idx_1_x_pt -= (1 << (lv - 1));
                box_idx_2_x_pt -= (1 << (lv - 1));
            }

            if (aux_idx_1_y > ptr_node->box_max_y)
            {
                box_idx_1_y_pt -= (1 << (lv - 1));
                box_idx_2_y_pt -= (1 << (lv - 1));
            }

            if (aux_idx_1_z > ptr_node->box_max_z)
            {
                box_idx_1_z_pt -= (1 << (lv - 1));
                box_idx_2_z_pt -= (1 << (lv - 1));
            }
        }

        box_grid_idx_1_pt = box_idx_1_x_pt + box_idx_1_y_pt * grid_box_real_dim_X_pt + box_idx_1_z_pt * grid_box_real_dim_X_times_Y_pt;
        box_grid_idx_2_pt = box_idx_2_x_pt + box_idx_1_y_pt * grid_box_real_dim_X_pt + box_idx_1_z_pt * grid_box_real_dim_X_times_Y_pt;
        box_grid_idx_3_pt = box_idx_1_x_pt + box_idx_2_y_pt * grid_box_real_dim_X_pt + box_idx_1_z_pt * grid_box_real_dim_X_times_Y_pt;
        box_grid_idx_4_pt = box_idx_2_x_pt + box_idx_2_y_pt * grid_box_real_dim_X_pt + box_idx_1_z_pt * grid_box_real_dim_X_times_Y_pt;
        box_grid_idx_5_pt = box_idx_1_x_pt + box_idx_1_y_pt * grid_box_real_dim_X_pt + box_idx_2_z_pt * grid_box_real_dim_X_times_Y_pt;
        box_grid_idx_6_pt = box_idx_2_x_pt + box_idx_1_y_pt * grid_box_real_dim_X_pt + box_idx_2_z_pt * grid_box_real_dim_X_times_Y_pt;
        box_grid_idx_7_pt = box_idx_1_x_pt + box_idx_2_y_pt * grid_box_real_dim_X_pt + box_idx_2_z_pt * grid_box_real_dim_X_times_Y_pt;
        box_grid_idx_8_pt = box_idx_2_x_pt + box_idx_2_y_pt * grid_box_real_dim_X_pt + box_idx_2_z_pt * grid_box_real_dim_X_times_Y_pt;

        // box_grid_idx_1_pt = box_idx_1_x_pt + box_idx_1_y_pt * grid_box_real_dim_X_pt + box_idx_1_z_pt * grid_box_real_dim_X_times_Y_pt;

        // //** x-axis acceleration **/
        // box_grid_idx_2_pt = box_idx_2_x_pt + box_idx_1_y_pt * grid_box_real_dim_X_pt + box_idx_1_z_pt * grid_box_real_dim_X_times_Y_pt;
        // aux_acc = (ptr_pt->ptr_ax[box_grid_idx_1_pt] + ptr_pt->ptr_ax[box_grid_idx_2_pt]) * 0.5;
        // ptr_node->ptr_ax[box_grid_idx] = aux_acc;

        // //** y-axis acceleration **/
        // box_grid_idx_2_pt = box_idx_1_x_pt + box_idx_2_y_pt * grid_box_real_dim_X_pt + box_idx_1_z_pt * grid_box_real_dim_X_times_Y_pt;
        // aux_acc = (ptr_pt->ptr_ay[box_grid_idx_1_pt] + ptr_pt->ptr_ay[box_grid_idx_2_pt]) * 0.5;
        // ptr_node->ptr_ay[box_grid_idx] = aux_acc;

        // //** z-axis acceleration **/
        // box_grid_idx_2_pt = box_idx_1_x_pt + box_idx_1_y_pt * grid_box_real_dim_X_pt + box_idx_2_z_pt * grid_box_real_dim_X_times_Y_pt;
        // aux_acc = (ptr_pt->ptr_az[box_grid_idx_1_pt] + ptr_pt->ptr_az[box_grid_idx_2_pt]) * 0.5;
        // ptr_node->ptr_az[box_grid_idx] = aux_acc;

        //** x-axis acceleration **/
        ptr_node->ptr_ax[box_grid_idx] = 0.125 * (ptr_pt->ptr_ax[box_grid_idx_1_pt] + ptr_pt->ptr_ax[box_grid_idx_2_pt] +
                                                  ptr_pt->ptr_ax[box_grid_idx_3_pt] + ptr_pt->ptr_ax[box_grid_idx_4_pt] +
                                                  ptr_pt->ptr_ax[box_grid_idx_5_pt] + ptr_pt->ptr_ax[box_grid_idx_6_pt] +
                                                  ptr_pt->ptr_ax[box_grid_idx_7_pt] + ptr_pt->ptr_ax[box_grid_idx_8_pt]);

        //** y-axis acceleration **/
        ptr_node->ptr_ay[box_grid_idx] = 0.125 * (ptr_pt->ptr_ay[box_grid_idx_1_pt] + ptr_pt->ptr_ay[box_grid_idx_2_pt] +
                                                  ptr_pt->ptr_ay[box_grid_idx_3_pt] + ptr_pt->ptr_ay[box_grid_idx_4_pt] +
                                                  ptr_pt->ptr_ay[box_grid_idx_5_pt] + ptr_pt->ptr_ay[box_grid_idx_6_pt] +
                                                  ptr_pt->ptr_ay[box_grid_idx_7_pt] + ptr_pt->ptr_ay[box_grid_idx_8_pt]);

        //** z-axis acceleration **/
        ptr_node->ptr_az[box_grid_idx] = 0.125 * (ptr_pt->ptr_az[box_grid_idx_1_pt] + ptr_pt->ptr_az[box_grid_idx_2_pt] +
                                                  ptr_pt->ptr_az[box_grid_idx_3_pt] + ptr_pt->ptr_az[box_grid_idx_4_pt] +
                                                  ptr_pt->ptr_az[box_grid_idx_5_pt] + ptr_pt->ptr_az[box_grid_idx_6_pt] +
                                                  ptr_pt->ptr_az[box_grid_idx_7_pt] + ptr_pt->ptr_az[box_grid_idx_8_pt]);
    }

    //** >> Interior grid points **/
    for (int i = 0; i < ptr_node->grid_intr_size; i++)
    {
        //** >> Box grid index **/
        box_grid_idx = ptr_node->ptr_intr_box_grid_idx[i];

        // box_idx_x_node = ptr_node->ptr_intr_grid_cell_idx_x[i] - ptr_node->box_ts_x;
        // box_idx_y_node = ptr_node->ptr_intr_grid_cell_idx_y[i] - ptr_node->box_ts_y;
        // box_idx_z_node = ptr_node->ptr_intr_grid_cell_idx_z[i] - ptr_node->box_ts_z;

        aux_idx_x = ptr_node->ptr_intr_grid_cell_idx_x[i];
        aux_idx_y = ptr_node->ptr_intr_grid_cell_idx_y[i];
        aux_idx_z = ptr_node->ptr_intr_grid_cell_idx_z[i];

        box_idx_x_node = aux_idx_x - ptr_node->box_ts_x;
        box_idx_y_node = aux_idx_y - ptr_node->box_ts_y;
        box_idx_z_node = aux_idx_z - ptr_node->box_ts_z;

        if (ptr_node->pbc_crosses_the_boundary_simulation_box == true)
        {
            if (aux_idx_x > ptr_node->box_max_x)
            {
                box_idx_x_node -= (1 << lv);
            }

            if (aux_idx_y > ptr_node->box_max_y)
            {
                box_idx_y_node -= (1 << lv);
            }

            if (aux_idx_z > ptr_node->box_max_z)
            {
                box_idx_z_node -= (1 << lv);
            }
        }

        box_idx_node = box_idx_x_node + box_idx_y_node * box_real_dim_X_node + box_idx_z_node * box_real_dim_X_times_Y_node;

        //** >> Neighboring box grid indexes **/
        box_grid_idxNbr_1 = box_grid_idx + 1;
        box_grid_idxNbr_2 = box_grid_idx - 1;
        box_grid_idxNbr_3 = box_grid_idx + grid_box_real_dim_X_node;
        box_grid_idxNbr_4 = box_grid_idx - grid_box_real_dim_X_node;
        box_grid_idxNbr_5 = box_grid_idx + grid_box_real_dim_X_times_Y_node;
        box_grid_idxNbr_6 = box_grid_idx - grid_box_real_dim_X_times_Y_node;

        //** >> Filling the acceleration **/
        // x
        check = false;
        if (ptr_node->ptr_box[box_idx_node - 2] > -4 && ptr_node->ptr_box[box_idx_node + 1] > -4)
        {
            box_grid_idxNbr_7 = box_grid_idx + 2;
            box_grid_idxNbr_8 = box_grid_idx - 2;
            check = true;
        }
        else if (ptr_node->pbc_crosses_the_whole_simulation_box_x == true)
        {
            if (ptr_node->ptr_box[box_idx_node - 2] == -6 && ptr_node->ptr_box[box_idx_node + 1] > -4)
            {
                if (ptr_node->ptr_box[box_idx_node - 2 + (1 << lv)] > -4)
                {
                    box_grid_idxNbr_7 = box_grid_idx + 2;
                    box_grid_idxNbr_8 = box_grid_idx - 2 + (1 << lv);
                    check = true;
                }
            }
            else if (ptr_node->ptr_box[box_idx_node + 1] == -6 && ptr_node->ptr_box[box_idx_node - 2] > -4)
            {
                box_grid_idxNbr_7 = box_grid_idx + 2 - (1 << lv);
                box_grid_idxNbr_8 = box_grid_idx - 2;
                check = true;
            } 
        }

        if (check == true)
        {
            ptr_node->ptr_ax[box_grid_idx] = (ptr_node->ptr_pot[box_grid_idxNbr_2] - ptr_node->ptr_pot[box_grid_idxNbr_1]) * two_over_3H - (ptr_node->ptr_pot[box_grid_idxNbr_8] - ptr_node->ptr_pot[box_grid_idxNbr_7]) * one_over_12H;
        }
        else
        {
            ptr_node->ptr_ax[box_grid_idx] = one_over_2H * (ptr_node->ptr_pot[box_grid_idxNbr_2] - ptr_node->ptr_pot[box_grid_idxNbr_1]);
        }

        // y
        check = false;
        if (ptr_node->ptr_box[box_idx_node - 2 * box_real_dim_X_node] > -4 && ptr_node->ptr_box[box_idx_node + box_real_dim_X_node] > -4)
        {
            box_grid_idxNbr_9  = box_grid_idx + 2 * grid_box_real_dim_X_node;
            box_grid_idxNbr_10 = box_grid_idx - 2 * grid_box_real_dim_X_node;
        }
        else if (ptr_node->pbc_crosses_the_whole_simulation_box_y == true)
        {
            if (ptr_node->ptr_box[box_idx_node - 2 * box_real_dim_X_node] == -6 && ptr_node->ptr_box[box_idx_node + box_real_dim_X_node] > -4)
            {
                if (ptr_node->ptr_box[box_idx_node + (-2 + (1 << lv)) * box_real_dim_X_node] > -4)
                {
                    box_grid_idxNbr_9  = box_grid_idx + 2 * grid_box_real_dim_X_node;
                    box_grid_idxNbr_10 = box_grid_idx + (-2 + (1 << lv)) * grid_box_real_dim_X_node;
                    check = true;
                }
            }
            else if (ptr_node->ptr_box[box_idx_node + box_real_dim_X_node] == -6 && ptr_node->ptr_box[box_idx_node - 2 * box_real_dim_X_node] > -4)
            {
                if (ptr_node->ptr_box[box_idx_node + (1 - (1 << lv)) * box_real_dim_X_node] > -4)
                {
                    box_grid_idxNbr_9  = box_grid_idx + (2 - (1 << lv)) * grid_box_real_dim_X_node;
                    box_grid_idxNbr_10 = box_grid_idx - 2 * grid_box_real_dim_X_node;
                    check = true;
                }

            }
        }

        if (check == true)
        {
            ptr_node->ptr_ay[box_grid_idx] = (ptr_node->ptr_pot[box_grid_idxNbr_4] - ptr_node->ptr_pot[box_grid_idxNbr_3]) * two_over_3H - (ptr_node->ptr_pot[box_grid_idxNbr_10] - ptr_node->ptr_pot[box_grid_idxNbr_9]) * one_over_12H;
        }
        else
        {
            ptr_node->ptr_ay[box_grid_idx] = one_over_2H * (ptr_node->ptr_pot[box_grid_idxNbr_4] - ptr_node->ptr_pot[box_grid_idxNbr_3]);
        }


        // z
        check = false;
        if (ptr_node->ptr_box[box_idx_node - 2 * box_real_dim_X_times_Y_node] > -4 && ptr_node->ptr_box[box_idx_node + box_real_dim_X_times_Y_node] > -4)
        {
            box_grid_idxNbr_11 = box_grid_idx + 2 * grid_box_real_dim_X_times_Y_node;
            box_grid_idxNbr_12 = box_grid_idx - 2 * grid_box_real_dim_X_times_Y_node;
        }
        else if (ptr_node->pbc_crosses_the_whole_simulation_box_y == true)
        {
            if (ptr_node->ptr_box[box_idx_node - 2 * grid_box_real_dim_X_times_Y_node] == -6 && ptr_node->ptr_box[box_idx_node + grid_box_real_dim_X_times_Y_node] > -4)
            {
                if (ptr_node->ptr_box[box_idx_node + (-2 + (1 << lv)) * grid_box_real_dim_X_times_Y_node] > -4)
                {
                    box_grid_idxNbr_11 = box_grid_idx + 2 * grid_box_real_dim_X_times_Y_node;
                    box_grid_idxNbr_12 = box_grid_idx + (-2 + (1 << lv)) * grid_box_real_dim_X_times_Y_node;
                    check = true;
                }
            }
            else if (ptr_node->ptr_box[box_idx_node + grid_box_real_dim_X_times_Y_node] == -6 && ptr_node->ptr_box[box_idx_node - 2 * grid_box_real_dim_X_times_Y_node] > -4)
            {
                if (ptr_node->ptr_box[box_idx_node + (1 - (1 << lv)) * grid_box_real_dim_X_times_Y_node] > -4)
                {
                    box_grid_idxNbr_11 = box_grid_idx + (2 - (1 << lv)) * grid_box_real_dim_X_times_Y_node;
                    box_grid_idxNbr_12 = box_grid_idx - 2 * grid_box_real_dim_X_times_Y_node;
                    check = true;
                }
            }
        }
        
        if (check == true)
        {
            ptr_node->ptr_az[box_grid_idx] = (ptr_node->ptr_pot[box_grid_idxNbr_6] - ptr_node->ptr_pot[box_grid_idxNbr_5]) * two_over_3H - (ptr_node->ptr_pot[box_grid_idxNbr_12] - ptr_node->ptr_pot[box_grid_idxNbr_11]) * one_over_12H;
        }
        else
        {
            ptr_node->ptr_az[box_grid_idx] = one_over_2H * (ptr_node->ptr_pot[box_grid_idxNbr_6] - ptr_node->ptr_pot[box_grid_idxNbr_5]);
        }
    }
}

int grid_acceleration(void)
{
    //** >> Acceleration in the grid **/

    //** >> HEAD
    if (force_stencil == 0)
    {
        //** >> HEAD
        computing_grid_acceleration_head(GL_tentacles[0][0]);
        
        //** >>  BRANCHES
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
    }
    else if (force_stencil == 1)
    {
        //** >> HEAD
        computing_grid_acceleration_head2(GL_tentacles[0][0]);

        //** >>  BRANCHES
        for (int lv = 1; lv < GL_tentacles_level_max + 1; lv++)
        {
            // Number of parent of the level = GL_tentacles_size[lv];

            //** >> For cycle over parent nodes **/
            for (int i = 0; i < GL_tentacles_size[lv]; i++)
            {
                // ptr_node = GL_tentacles[lv][i];
                computing_grid_acceleration_branch2(GL_tentacles[lv][i]);
            }
        }
    }
    else
    {
        printf("Error, the force stencil needs to be equal to 0 or 1\n");
        return _FAILURE_;
    }

    return _SUCCESS_;
}