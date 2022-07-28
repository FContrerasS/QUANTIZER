/*
 * tree_construction.c
 *
 * Create the initial tree
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

#include "tree_construction.h"

//** >> Local Functions
static int fill_cell_ref_PERIODIC_BOUNDARY(struct node *ptr_node);
static int fill_cell_ref(struct node *ptr_node);
static int fill_zones_ref_PERIODIC_BOUNDARY(struct node *ptr_node);
static int fill_zones_ref(struct node *ptr_node);
static int fill_subzones_ref_PERIODIC_BOUNDARY(struct node *ptr_node, int zone_idx);
static int fill_child_nodes_PERIODIC_BOUNDARY(struct node *ptr_node);
static int fill_child_nodes(struct node *ptr_node);
static int fill_tentacles(const struct node *ptr_node);

static int fill_cell_ref_PERIODIC_BOUNDARY(struct node *ptr_node)
{
    //** >> Adding cells which satisfy the refinement criteria to the array ptr_cell_ref and chaning the box the status of refinement -1 **/
    int box_idx_node;    // Box index
    int box_idxNbr_node; // Box index in the neigborhood

    int cell_ref_idx = 0; // Index of the position in the cell refined array

    int box_real_dim_X_node = ptr_node->box_real_dim_x;
    int box_real_dim_X_times_Y_node = ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;

    int ii_aux, jj_aux, kk_aux;

    //** >> Changing the box status from EXIST (-3) to REFINEMENT REQUIRED (-1) **/
    for (int i = 0; i < ptr_node->cell_size; i++)
    {
        box_idx_node = ptr_node->ptr_box_idx[i];

        // Refinement criterion in the box_mass in no border box points
        if (ptr_node->ptr_cell_struct[box_idx_node].cell_mass >= ref_criterion_mass || ptr_node->ptr_cell_struct[box_idx_node].ptcl_size >= ref_criterion_ptcl)
        {
            if (ptr_node->ptr_box[box_idx_node] == -3) // Cell has not been added yet
            {
                cell_ref_idx++;
                //** >> Chaning the cell box status from EXIST (-3) to REFINEMENT REQUIRED (-1) **/
                ptr_node->ptr_box[box_idx_node] = -1;
            }

            //** >> Changing the neighboring cell status **/
            for (int kk = -n_exp; kk < n_exp + 1; kk++)
            {
                for (int jj = -n_exp; jj < n_exp + 1; jj++)
                {
                    for (int ii = -n_exp; ii < n_exp + 1; ii++)
                    {
                        box_idxNbr_node = box_idx_node + ii + jj * box_real_dim_X_node + kk * box_real_dim_X_times_Y_node;

                        if (ptr_node->ptr_box[box_idxNbr_node] == -5)
                        {
                            if (ptr_node->ptr_box[box_idx_node + ii] == -5)
                            {
                                ii_aux = ii < 0 ? ii + ptr_node->box_dim_x : ii - ptr_node->box_dim_x;
                            }
                            else
                            {
                                ii_aux = ii;
                            }

                            if (ptr_node->ptr_box[box_idx_node + jj * box_real_dim_X_node] == -5)
                            {
                                jj_aux = jj < 0 ? jj + ptr_node->box_dim_y : jj - ptr_node->box_dim_y;
                            }
                            else
                            {
                                jj_aux = jj;
                            }

                            if (ptr_node->ptr_box[box_idx_node + kk * box_real_dim_X_times_Y_node] == -5)
                            {
                                kk_aux = kk < 0 ? kk + ptr_node->box_dim_z : kk - ptr_node->box_dim_z;
                            }
                            else
                            {
                                kk_aux = kk;
                            }

                            box_idxNbr_node = box_idx_node + ii_aux + jj_aux * box_real_dim_X_node + kk_aux * box_real_dim_X_times_Y_node;
                        }

                        //** >> Asking if the neighboring cell has not been changed yet **/
                        // The border (-2) of the simulation can be added
                        if (ptr_node->ptr_box[box_idxNbr_node] == -3) // Cell has not been added yet
                        {
                            cell_ref_idx++;
                            //** >> Chaning the cell box status from EXIST (-3) to REFINEMENT REQUIRED (-1) **/
                            ptr_node->ptr_box[box_idxNbr_node] = -1;
                        }
                    }
                }
            }
        }
    }

    //** >> Adding the infomation about size of the ptr_cell_ref array **/
    ptr_node->cell_ref_size = cell_ref_idx;
    ptr_node->cell_ref_cap = cell_ref_idx;
    ptr_node->ptr_cell_ref = (int *)malloc(ptr_node->cell_ref_cap * sizeof(int));

    //** >> Adding cells refined to the array ptr_cell_ref **/
    cell_ref_idx = 0;
    for (int i = 0; i < ptr_node->cell_size; i++)
    {
        box_idx_node = ptr_node->ptr_box_idx[i];

        if (ptr_node->ptr_box[box_idx_node] == -1) // Cell require refinement
        {
            ptr_node->ptr_cell_ref[cell_ref_idx] = i;
            cell_ref_idx++;
        }
    }

    return _SUCCESS_;
} // end function fill_cell_ref

static int fill_cell_ref(struct node *ptr_node)
{
    //** >> Adding cells which satisfy the refinement criteria to the array ptr_cell_ref and chaning the box the status of refinement -1 **/
    int box_idx_node;    // Box index
    int box_idxNbr_node; // Box index in the neigborhood

    int cell_ref_idx = 0; // Index of the position in the cell refined array

    int box_real_dim_X_node = ptr_node->box_real_dim_x;
    int box_real_dim_X_times_Y_node = ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;

    //** >> Changing the box status from EXIST (-3) to REFINEMENT REQUIRED (-1) **/
    for (int i = 0; i < ptr_node->cell_size; i++)
    {
        box_idx_node = ptr_node->ptr_box_idx[i];

        // Refinement criterion in the box_mass in no border box points
        if (ptr_node->ptr_cell_struct[box_idx_node].cell_mass >= ref_criterion_mass || ptr_node->ptr_cell_struct[box_idx_node].ptcl_size >= ref_criterion_ptcl) 
        {
            if (ptr_node->ptr_box[box_idx_node] == -3) // Cell has not been added yet
            {
                cell_ref_idx++;
                //** >> Chaning the cell box status from EXIST (-3) to REFINEMENT REQUIRED (-1) **/
                ptr_node->ptr_box[box_idx_node] = -1;
            }

            //** >> Changing the neighboring cell status **/
            for (int kk = -n_exp; kk < n_exp + 1; kk++)
            {
                for (int jj = -n_exp; jj < n_exp + 1; jj++)
                {
                    for (int ii = -n_exp; ii < n_exp + 1; ii++)
                    {
                        box_idxNbr_node = box_idx_node + ii + jj * box_real_dim_X_node + kk * box_real_dim_X_times_Y_node;
                        //** >> Asking if the neighboring cell has not been changed yet **/
                        // The border (-2) of the simulation can be added
                        if (ptr_node->ptr_box[box_idxNbr_node] == -3) // Cell has not been added yet
                        {
                            cell_ref_idx++;
                            //** >> Chaning the cell box status from EXIST (-3) to REFINEMENT REQUIRED (-1) **/
                            ptr_node->ptr_box[box_idxNbr_node] = -1;
                        }
                    }
                }
            }
        }
    }

    //** >> Adding the infomation about size of the ptr_cell_ref array **/
    ptr_node->cell_ref_size = cell_ref_idx;
    ptr_node->cell_ref_cap = cell_ref_idx;
    ptr_node->ptr_cell_ref = (int *)malloc(ptr_node->cell_ref_cap * sizeof(int));

    //** >> Adding cells refined to the array ptr_cell_ref **/
    cell_ref_idx = 0;
    for (int i = 0; i < ptr_node->cell_size; i++)
    {
        box_idx_node = ptr_node->ptr_box_idx[i];

        if (ptr_node->ptr_box[box_idx_node] == -1) // Cell require refinement
        {
            ptr_node->ptr_cell_ref[cell_ref_idx] = i;
            cell_ref_idx++;
        }
    }

    return _SUCCESS_;
} // end function fill_cell_ref

static int fill_zones_ref_PERIODIC_BOUNDARY(struct node *ptr_node)
{
    //** >> Filling the auxiliary diferent refinement zones **/
    int cntr_cell_add_all_zones = 0; // Counter Number of cells added to any refinement zone
    int cntr_cell_add;               // Counter number of cells added per inspection
    int cntr_insp;                   // Numer of inspected cells in the zone

    int zone_size;                          // Number of cells in the zone
    int zone_idx_max = ptr_node->zones_cap; // Maximum id of the zone. It is equal to to the capacity in the number of zones
    int zone_idx = 0;                       // Index of the zone. Initialized at zone 0

    int box_idx_node; // Box index

    int box_idxNbr_x_plus;  // Box index in the neigborhood on the right
    int box_idxNbr_x_minus; // Box index in the neigborhood on the left
    int box_idxNbr_y_plus;  // Box index in the neigborhood behind
    int box_idxNbr_y_minus; // Box index in the neigborhood in front
    int box_idxNbr_z_plus;  // Box index in the neigborhood up
    int box_idxNbr_z_minus; // Box index in the neigborhood down

    int cell_idx;         // The cell index is simply i of the for loop
    int cell_ref_idx = 0; // The index in the cell refined array ptr_cell_ref

    int box_real_dim_X_node = ptr_node->box_real_dim_x;
    int box_real_dim_X_times_Y_node = ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;

    // The auxiliary array ptr_aux_idx will be used to store the box indexes of the zone

    //** >> Space checking of auxiliary array ptr_aux_idx **/
    // The size will always be bigger or equal than the number of refined cells

    if (space_check(&(ptr_node->aux_idx_cap), ptr_node->cell_ref_cap, 1.0f, "p7i1b1b1b1b1b1b1", &(ptr_node->ptr_aux_idx), &(ptr_node->ptr_aux_bool_boundary_simulation_contact_x), &(ptr_node->ptr_aux_bool_boundary_simulation_contact_y), &(ptr_node->ptr_aux_bool_boundary_simulation_contact_z), &(ptr_node->ptr_aux_bool_boundary_anomalies_x), &(ptr_node->ptr_aux_bool_boundary_anomalies_y), &(ptr_node->ptr_aux_bool_boundary_anomalies_z)) == _FAILURE_)
    {
        printf("Error, in space_check function\n");
        return _FAILURE_;
    }

    //** >>  Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0) **/
    while (ptr_node->cell_ref_size > cntr_cell_add_all_zones) // The loop while works as long as the number of cell addeed is less than the total refined cells
    {
        // Notes that the initiality we inspect the elements in the refined cell array until an element has been found that is not found in any of the current refinement zones
        cell_idx = ptr_node->ptr_cell_ref[cell_ref_idx]; // Index of the cells array in the node
        box_idx_node = ptr_node->ptr_box_idx[cell_idx];

        if (ptr_node->ptr_box[box_idx_node] == -1) // A cell without zone has been founded
        {
            zone_size = 0; // Initial number of element in the zone

            //** >> Including the first element of the box to the auxiliary array ptr_aux_idx **/
            ptr_node->ptr_aux_idx[0] = box_idx_node;

            //** >>  Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0) **/
            ptr_node->ptr_box[box_idx_node] = zone_idx;

            zone_size++;               // +1 to the number of cells in the zone
            cntr_cell_add_all_zones++; // +1 to the number of cells added in total

            //** >> Building of the zone of ID = zone_idx **/
            cntr_insp = 0;                                                                     // Counter for the number of elements inspectioned in the current zone array
            while (zone_size > cntr_insp && ptr_node->cell_ref_size > cntr_cell_add_all_zones) // The cycle ends when all the elements of the zone have been inspected or when the number of cells refined is equal to the total number of cells added.
            {
                // Note that the number of elements in the zone increases if the neighbors of the inspected cell must be added to the zone

                cntr_cell_add = 0; // Counter number of cells added per cell inspected
                box_idx_node = ptr_node->ptr_aux_idx[cntr_insp];

                box_idxNbr_x_plus = box_idx_node + 1;
                box_idxNbr_x_minus = box_idx_node - 1;
                box_idxNbr_y_plus = box_idx_node + box_real_dim_X_node;
                box_idxNbr_y_minus = box_idx_node - box_real_dim_X_node;
                box_idxNbr_z_plus = box_idx_node + box_real_dim_X_times_Y_node;
                box_idxNbr_z_minus = box_idx_node - box_real_dim_X_times_Y_node;

                if (ptr_node->ptr_box[box_idxNbr_x_plus] == -5)
                {
                    ptr_node->ptr_aux_bool_boundary_simulation_contact_x[zone_idx] = true;
                    box_idxNbr_x_plus -= ptr_node->box_dim_x;
                    if (ptr_node->ptr_box[box_idxNbr_x_plus] == -1 || ptr_node->ptr_box[box_idxNbr_x_plus] == zone_idx)
                    {
                        ptr_node->ptr_aux_bool_boundary_anomalies_x[zone_idx] = true;
                    }
                }
                else if (ptr_node->ptr_box[box_idxNbr_x_minus] == -5)
                {
                    ptr_node->ptr_aux_bool_boundary_simulation_contact_x[zone_idx] = true;
                    box_idxNbr_x_minus += ptr_node->box_dim_x;
                    if (ptr_node->ptr_box[box_idxNbr_x_minus] == -1 || ptr_node->ptr_box[box_idxNbr_x_minus] == zone_idx)
                    {
                        ptr_node->ptr_aux_bool_boundary_anomalies_x[zone_idx] = true;
                    }
                }

                if (ptr_node->ptr_box[box_idxNbr_y_plus] == -5)
                {
                    ptr_node->ptr_aux_bool_boundary_simulation_contact_y[zone_idx] = true;
                    box_idxNbr_y_plus -= ptr_node->box_dim_y * box_real_dim_X_node;
                    if (ptr_node->ptr_box[box_idxNbr_y_plus] == -1 || ptr_node->ptr_box[box_idxNbr_y_plus] == zone_idx)
                    {
                        ptr_node->ptr_aux_bool_boundary_anomalies_y[zone_idx] = true;
                    }
                }
                else if (ptr_node->ptr_box[box_idxNbr_y_minus] == -5)
                {
                    ptr_node->ptr_aux_bool_boundary_simulation_contact_y[zone_idx] = true;
                    box_idxNbr_y_minus += ptr_node->box_dim_y * box_real_dim_X_node;
                    if (ptr_node->ptr_box[box_idxNbr_y_minus] == -1 || ptr_node->ptr_box[box_idxNbr_y_minus] == zone_idx)
                    {
                        ptr_node->ptr_aux_bool_boundary_anomalies_y[zone_idx] = true;
                    }
                }

                if (ptr_node->ptr_box[box_idxNbr_z_plus] == -5)
                {
                    ptr_node->ptr_aux_bool_boundary_simulation_contact_z[zone_idx] = true;
                    box_idxNbr_z_plus -= ptr_node->box_dim_z * box_real_dim_X_times_Y_node;
                    if (ptr_node->ptr_box[box_idxNbr_z_plus] == -1 || ptr_node->ptr_box[box_idxNbr_z_plus] == zone_idx)
                    {
                        ptr_node->ptr_aux_bool_boundary_anomalies_z[zone_idx] = true;
                    }
                }
                else if (ptr_node->ptr_box[box_idxNbr_z_minus] == -5)
                {
                    ptr_node->ptr_aux_bool_boundary_simulation_contact_z[zone_idx] = true;
                    box_idxNbr_z_minus += ptr_node->box_dim_z * box_real_dim_X_times_Y_node;
                    if (ptr_node->ptr_box[box_idxNbr_z_minus] == -1 || ptr_node->ptr_box[box_idxNbr_z_minus] == zone_idx)
                    {
                        ptr_node->ptr_aux_bool_boundary_anomalies_z[zone_idx] = true;
                    }
                }

                //** Checking the nearest 6 neighbors of face
                // First neighbor
                if (ptr_node->ptr_box[box_idxNbr_x_plus] == -1)
                {
                    ptr_node->ptr_aux_idx[zone_size + cntr_cell_add] = box_idxNbr_x_plus; // Including the neighboring element of the box to the auxiliary array
                    ptr_node->ptr_box[box_idxNbr_x_plus] = zone_idx;                      // Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0)
                    cntr_cell_add++;                                                      // +1 to the number of cells added in the current inspection
                }
                // Second neighbor
                if (ptr_node->ptr_box[box_idxNbr_x_minus] == -1)
                {
                    ptr_node->ptr_aux_idx[zone_size + cntr_cell_add] = box_idxNbr_x_minus; // Including the neighboring element of the box to the auxiliary array
                    ptr_node->ptr_box[box_idxNbr_x_minus] = zone_idx;                      // Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0)
                    cntr_cell_add++;                                                       // +1 to the number of cells added in the current inspection
                }
                // Third neighbor
                if (ptr_node->ptr_box[box_idxNbr_y_plus] == -1)
                {
                    ptr_node->ptr_aux_idx[zone_size + cntr_cell_add] = box_idxNbr_y_plus; // Including the neighboring element of the box to the auxiliary array
                    ptr_node->ptr_box[box_idxNbr_y_plus] = zone_idx;                      // Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0)
                    cntr_cell_add++;                                                      // +1 to the number of cells added in the current inspection
                }
                // Fourth neighbor
                if (ptr_node->ptr_box[box_idxNbr_y_minus] == -1)
                {
                    ptr_node->ptr_aux_idx[zone_size + cntr_cell_add] = box_idxNbr_y_minus; // Including the neighboring element of the box to the auxiliary array
                    ptr_node->ptr_box[box_idxNbr_y_minus] = zone_idx;                      // Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0)
                    cntr_cell_add++;                                                       // +1 to the number of cells added in the current inspection
                }
                // Fifth neighbor
                if (ptr_node->ptr_box[box_idxNbr_z_plus] == -1)
                {
                    ptr_node->ptr_aux_idx[zone_size + cntr_cell_add] = box_idxNbr_z_plus; // Including the neighboring element of the box to the auxiliary array
                    ptr_node->ptr_box[box_idxNbr_z_plus] = zone_idx;                      // Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0)
                    cntr_cell_add++;                                                      // +1 to the number of cells added in the current inspection
                }
                // Sixth neighbor
                if (ptr_node->ptr_box[box_idxNbr_z_minus] == -1)
                {
                    ptr_node->ptr_aux_idx[zone_size + cntr_cell_add] = box_idxNbr_z_minus; // Including the neighboring element of the box to the auxiliary array
                    ptr_node->ptr_box[box_idxNbr_z_minus] = zone_idx;                      // Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0)
                    cntr_cell_add++;                                                       // +1 to the number of cells added in the current inspection
                }
                zone_size += cntr_cell_add;               // Increasing the number of cells in the zone
                cntr_cell_add_all_zones += cntr_cell_add; // Increasing the number of cells added to the zones
                cntr_insp++;                              // Increasing the number of inspections
            }                                             // End while cycle, now the box contains the the information about all cells of the zone "zone_idx"

            //** >> Space checking of refinement zones arrays, refinement capacity array and refinement size array **/
            if (zone_idx_max < zone_idx + 1)
            {
                if (space_check(&(zone_idx_max), zone_idx + 1, 2.0f, "p3i2i1i1", &(ptr_node->pptr_zones), &(ptr_node->ptr_zone_cap), &(ptr_node->ptr_zone_size)) == _FAILURE_)
                {
                    printf("Error, in space_check function\n");
                    return _FAILURE_;
                }
            }
            ptr_node->ptr_zone_size[zone_idx] = zone_size;
            zone_idx++; // Increasing the zone number
        }               // Zone defined in the box
        cell_ref_idx++;
    } // At this point the box contains the information of all refinement zones

    ptr_node->zones_cap = zone_idx_max; // Maximum amount of zones
    ptr_node->zones_size = zone_idx;    // Total amount of zones

    //** >> Space checking in each zone **/
    for (int i = 0; i < ptr_node->zones_size; i++)
    {
        if (space_check(&(ptr_node->ptr_zone_cap[i]), ptr_node->ptr_zone_size[i], 2.0f, "p1i1", &(ptr_node->pptr_zones[i])) == _FAILURE_)
        {
            printf("Error, in space_check function\n");
            return _FAILURE_;
        }
    }

    //** >> Initializing the ptr_aux_idx array to be used as a counter of elemnents in each zone**/
    // Notes that the number of refined cells is always bigger or equal than the number of refined zones, so that the capacity of ptr_aux_idx is always enough to counter the number of refined zones
    for (int i = 0; i < ptr_node->zones_size; i++)
    {
        ptr_node->ptr_aux_idx[i] = 0;
    }

    //** >> Adding the cells to the zone array pptr_zones **/
    for (int i = 0; i < ptr_node->cell_ref_size; i++)
    {
        cell_idx = ptr_node->ptr_cell_ref[i];
        box_idx_node = ptr_node->ptr_box_idx[cell_idx];
        zone_idx = ptr_node->ptr_box[box_idx_node];
        cntr_cell_add = ptr_node->ptr_aux_idx[zone_idx];          // Counter the element in the zone "zone_idx"
        ptr_node->pptr_zones[zone_idx][cntr_cell_add] = cell_idx; // Adding the index of the cell array in the block to the zone
        ptr_node->ptr_aux_idx[zone_idx] += 1;                     // Counter the number of elements added in the zone "zone_idx"
    }

    return _SUCCESS_;
} // End function fill_zones_ref

static int fill_zones_ref(struct node *ptr_node)
{

    //** >> Filling the auxiliary diferent refinement zones **/
    int cntr_cell_add_all_zones = 0; // Counter Number of cells added to any refinement zone
    int cntr_cell_add;               // Counter number of cells added per inspection
    int cntr_insp;                   // Numer of inspected cells in the zone

    int zone_size;                          // Number of cells in the zone
    int zone_idx_max = ptr_node->zones_cap; // Maximum id of the zone. It is equal to to the capacity in the number of zones
    int zone_idx = 0;                       // Index of the zone. Initialized at zone 0

    int box_idx_node; // Box index

    int box_idxNbr_x_plus;  // Box index in the neigborhood on the right
    int box_idxNbr_x_minus; // Box index in the neigborhood on the left
    int box_idxNbr_y_plus;  // Box index in the neigborhood behind
    int box_idxNbr_y_minus; // Box index in the neigborhood in front
    int box_idxNbr_z_plus;  // Box index in the neigborhood up
    int box_idxNbr_z_minus; // Box index in the neigborhood down

    int cell_idx;         // The cell index is simply i of the for loop
    int cell_ref_idx = 0; // The index in the cell refined array ptr_cell_ref

    int box_real_dim_X_node = ptr_node->box_real_dim_x;
    int box_real_dim_X_times_Y_node = ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;

    // The auxiliary array ptr_aux_idx will be used to store the box indexes of the zone

    //** >> Space checking of auxiliary array ptr_aux_idx **/
    // The size will always be bigger or equal than the number of refined cells

    if (space_check(&(ptr_node->aux_idx_cap), ptr_node->cell_ref_cap, 1.0f, "p1i1", &(ptr_node->ptr_aux_idx)) == _FAILURE_)
    {
        printf("Error, in space_check function\n");
        return _FAILURE_;
    }

    //** >>  Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0) **/
    while (ptr_node->cell_ref_size > cntr_cell_add_all_zones) // The loop while works as long as the number of cell addeed is less than the total refined cells
    {
        // Notes that the initiality we inspect the elements in the refined cell array until an element has been found that is not found in any of the current refinement zones
        cell_idx = ptr_node->ptr_cell_ref[cell_ref_idx]; // Index of the cells array in the node
        box_idx_node = ptr_node->ptr_box_idx[cell_idx];

        if (ptr_node->ptr_box[box_idx_node] == -1) // A cell without zone has been founded
        {
            zone_size = 0; // Initial number of element in the zone

            //** >> Including the first element of the box to the auxiliary array ptr_aux_idx **/
            ptr_node->ptr_aux_idx[0] = box_idx_node;

            //** >>  Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0) **/
            ptr_node->ptr_box[box_idx_node] = zone_idx;

            zone_size++;               // +1 to the number of cells in the zone
            cntr_cell_add_all_zones++; // +1 to the number of cells added in total

            //** >> Building of the zone of ID = zone_idx **/
            cntr_insp = 0;                                                                     // Counter for the number of elements inspectioned in the current zone array
            while (zone_size > cntr_insp && ptr_node->cell_ref_size > cntr_cell_add_all_zones) // The cycle ends when all the elements of the zone have been inspected or when the number of cells refined is equal to the total number of cells added.
            {
                // Note that the number of elements in the zone increases if the neighbors of the inspected cell must be added to the zone

                cntr_cell_add = 0; // Counter number of cells added per cell inspected
                box_idx_node = ptr_node->ptr_aux_idx[cntr_insp];

                box_idxNbr_x_plus = box_idx_node + 1;
                box_idxNbr_x_minus = box_idx_node - 1;
                box_idxNbr_y_plus = box_idx_node + box_real_dim_X_node;
                box_idxNbr_y_minus = box_idx_node - box_real_dim_X_node;
                box_idxNbr_z_plus = box_idx_node + box_real_dim_X_times_Y_node;
                box_idxNbr_z_minus = box_idx_node - box_real_dim_X_times_Y_node;

                //** Checking the nearest 6 neighbors of face
                // First neighbor
                if (ptr_node->ptr_box[box_idxNbr_x_plus] == -1)
                {
                    ptr_node->ptr_aux_idx[zone_size + cntr_cell_add] = box_idxNbr_x_plus; // Including the neighboring element of the box to the auxiliary array
                    ptr_node->ptr_box[box_idxNbr_x_plus] = zone_idx;                      // Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0)
                    cntr_cell_add++;                                                      // +1 to the number of cells added in the current inspection
                }
                // Second neighbor
                if (ptr_node->ptr_box[box_idxNbr_x_minus] == -1)
                {
                    ptr_node->ptr_aux_idx[zone_size + cntr_cell_add] = box_idxNbr_x_minus; // Including the neighboring element of the box to the auxiliary array
                    ptr_node->ptr_box[box_idxNbr_x_minus] = zone_idx;                      // Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0)
                    cntr_cell_add++;                                                       // +1 to the number of cells added in the current inspection
                }
                // Third neighbor
                if (ptr_node->ptr_box[box_idxNbr_y_plus] == -1)
                {
                    ptr_node->ptr_aux_idx[zone_size + cntr_cell_add] = box_idxNbr_y_plus; // Including the neighboring element of the box to the auxiliary array
                    ptr_node->ptr_box[box_idxNbr_y_plus] = zone_idx;                      // Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0)
                    cntr_cell_add++;                                                      // +1 to the number of cells added in the current inspection
                }
                // Fourth neighbor
                if (ptr_node->ptr_box[box_idxNbr_y_minus] == -1)
                {
                    ptr_node->ptr_aux_idx[zone_size + cntr_cell_add] = box_idxNbr_y_minus; // Including the neighboring element of the box to the auxiliary array
                    ptr_node->ptr_box[box_idxNbr_y_minus] = zone_idx;                      // Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0)
                    cntr_cell_add++;                                                       // +1 to the number of cells added in the current inspection
                }
                // Fifth neighbor
                if (ptr_node->ptr_box[box_idxNbr_z_plus] == -1)
                {
                    ptr_node->ptr_aux_idx[zone_size + cntr_cell_add] = box_idxNbr_z_plus; // Including the neighboring element of the box to the auxiliary array
                    ptr_node->ptr_box[box_idxNbr_z_plus] = zone_idx;                      // Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0)
                    cntr_cell_add++;                                                      // +1 to the number of cells added in the current inspection
                }
                // Sixth neighbor
                if (ptr_node->ptr_box[box_idxNbr_z_minus] == -1)
                {
                    ptr_node->ptr_aux_idx[zone_size + cntr_cell_add] = box_idxNbr_z_minus; // Including the neighboring element of the box to the auxiliary array
                    ptr_node->ptr_box[box_idxNbr_z_minus] = zone_idx;                      // Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0)
                    cntr_cell_add++;                                                       // +1 to the number of cells added in the current inspection
                }
                zone_size += cntr_cell_add;               // Increasing the number of cells in the zone
                cntr_cell_add_all_zones += cntr_cell_add; // Increasing the number of cells added to the zones
                cntr_insp++;                              // Increasing the number of inspections
            }                                             // End while cycle, now the box contains the the information about all cells of the zone "zone_idx"

            //** >> Space checking of refinement zones arrays, refinement capacity array and refinement size array **/
            if (zone_idx_max < zone_idx + 1)
            {
                if (space_check(&(zone_idx_max), zone_idx + 1, 2.0f, "p3i2i1i1", &(ptr_node->pptr_zones), &(ptr_node->ptr_zone_cap), &(ptr_node->ptr_zone_size)) == _FAILURE_)
                {
                    printf("Error, in space_check function\n");
                    return _FAILURE_;
                }
            }
            ptr_node->ptr_zone_size[zone_idx] = zone_size;
            zone_idx++; // Increasing the zone number
        }               // Zone defined in the box
        cell_ref_idx++;
    } // At this point the box contains the information of all refinement zones

    ptr_node->zones_cap = zone_idx_max; // Maximum amount of zones
    ptr_node->zones_size = zone_idx;    // Total amount of zones

    //** >> Space checking in each zone **/
    for (int i = 0; i < ptr_node->zones_size; i++)
    {
        if (space_check(&(ptr_node->ptr_zone_cap[i]), ptr_node->ptr_zone_size[i], 2.0f, "p1i1", &(ptr_node->pptr_zones[i])) == _FAILURE_)
        {
            printf("Error, in space_check function\n");
            return _FAILURE_;
        }
    }

    //** >> Initializing the ptr_aux_idx array to be used as a counter of elemnents in each zone**/
    // Notes that the number of refined cells is always bigger or equal than the number of refined zones, so that the capacity of ptr_aux_idx is always enough to counter the number of refined zones
    for (int i = 0; i < ptr_node->zones_size; i++)
    {
        ptr_node->ptr_aux_idx[i] = 0;
    }

    //** >> Adding the cells to the zone array pptr_zones **/
    for (int i = 0; i < ptr_node->cell_ref_size; i++)
    {
        cell_idx = ptr_node->ptr_cell_ref[i];
        box_idx_node = ptr_node->ptr_box_idx[cell_idx];
        zone_idx = ptr_node->ptr_box[box_idx_node];
        cntr_cell_add = ptr_node->ptr_aux_idx[zone_idx];          // Counter the element in the zone "zone_idx"
        ptr_node->pptr_zones[zone_idx][cntr_cell_add] = cell_idx; // Adding the index of the cell array in the block to the zone
        ptr_node->ptr_aux_idx[zone_idx] += 1;                     // Counter the number of elements added in the zone "zone_idx"
    }

    return _SUCCESS_;
} // End function fill_zones_ref

static int fill_subzones_ref_PERIODIC_BOUNDARY(struct node *ptr_node, int zone_idx)
{
    //** >> Filling the auxiliary diferent refinement subzones **/
    int cntr_cell_add_all_subzones = 0; // Counter Number of cells added to any refinement subzone
    int cntr_cell_add;               // Counter number of cells added per inspection
    int cntr_insp;                   // Numer of inspected cells in the subzone

    int subzone_size;                          // Number of cells in the subzone
    int subzone_idx_max = ptr_node->subzones_cap; // Maximum id of the subzone. It is equal to to the capacity in the number of subzones
    int subzone_idx = 0;                       // Index of the subzone. Initialized at subzone 0

    int box_idx_node; // Box index

    int box_idxNbr_x_plus;  // Box index in the neigborhood on the right
    int box_idxNbr_x_minus; // Box index in the neigborhood on the left
    int box_idxNbr_y_plus;  // Box index in the neigborhood behind
    int box_idxNbr_y_minus; // Box index in the neigborhood in front
    int box_idxNbr_z_plus;  // Box index in the neigborhood up
    int box_idxNbr_z_minus; // Box index in the neigborhood down

    int cell_idx;         // The cell index is simply i of the for loop
    int cell_ref_idx = 0; // The index in the cell refined array pptr_zones[zone_idx]

    int box_real_dim_X_node = ptr_node->box_real_dim_x;
    int box_real_dim_X_times_Y_node = ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;

    // The auxiliary array ptr_aux_idx will be used to store the box indexes of the subzone

    //** >>  Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0) **/
    while (ptr_node->ptr_zone_size[zone_idx] > cntr_cell_add_all_subzones) // The loop while works as long as the number of cell addeed is less than the total refined cells
    {
        // Notes that the initiality we inspect the elements in the refined cell array until an element has been found that is not found in any of the current refinement subzones
        cell_idx = ptr_node->pptr_zones[zone_idx][cell_ref_idx]; // Index of the cells array in the node
        box_idx_node = ptr_node->ptr_box_idx[cell_idx];

        if (ptr_node->ptr_box[box_idx_node] == -1) // A cell without subzone has been founded
        {
            subzone_size = 0; // Initial number of element in the subzone

            //** >> Including the first element of the box to the auxiliary array ptr_aux_idx **/
            ptr_node->ptr_aux_idx[0] = box_idx_node;

            //** >>  Changing the box status from zone ID (>= 0) to the REFINEMENT REQUIRED (-1) in order not to repeat elements **/
            ptr_node->ptr_box[box_idx_node] = -1;

            subzone_size++;               // +1 to the number of cells in the zone
            cntr_cell_add_all_subzones++; // +1 to the number of cells added in total

            //** >> Building of the subzone of ID = subzone_idx **/
            cntr_insp = 0;                                                                     // Counter for the number of elements inspectioned in the current subzone array
            while (subzone_size > cntr_insp && ptr_node->cell_ref_size > cntr_cell_add_all_subzones) // The cycle ends when all the elements of the subzone have been inspected or when the number of cells in the zone is equal to the total number of cells added in all the subzones.
            {
                // Note that the number of elements in the subzone increases if the neighbors of the inspected cell must be added to the subzone

                cntr_cell_add = 0; // Counter number of cells added per cell inspected
                box_idx_node = ptr_node->ptr_aux_idx[cntr_insp];

                box_idxNbr_x_plus = box_idx_node + 1;
                box_idxNbr_x_minus = box_idx_node - 1;
                box_idxNbr_y_plus = box_idx_node + box_real_dim_X_node;
                box_idxNbr_y_minus = box_idx_node - box_real_dim_X_node;
                box_idxNbr_z_plus = box_idx_node + box_real_dim_X_times_Y_node;
                box_idxNbr_z_minus = box_idx_node - box_real_dim_X_times_Y_node;

                //** Checking the nearest 6 neighbors of face
                // First neighbor
                if (ptr_node->ptr_box[box_idxNbr_x_plus] == zone_idx)
                {
                    ptr_node->ptr_aux_idx[subzone_size + cntr_cell_add] = box_idxNbr_x_plus; // Including the neighboring element of the box to the auxiliary array
                    ptr_node->ptr_box[box_idxNbr_x_plus] = -1;                               // Changing the box status from refinement zone ID (>= 0) to the REFINEMENT REQUIRED (-1)
                    cntr_cell_add++;                                                      // +1 to the number of cells added in the current inspection
                }
                // Second neighbor
                if (ptr_node->ptr_box[box_idxNbr_x_minus] == zone_idx)
                {
                    ptr_node->ptr_aux_idx[subzone_size + cntr_cell_add] = box_idxNbr_x_minus; // Including the neighboring element of the box to the auxiliary array
                    ptr_node->ptr_box[box_idxNbr_x_minus] = -1;                               // Changing the box status from refinement zone ID (>= 0) to the REFINEMENT REQUIRED (-1)
                    cntr_cell_add++;                                                       // +1 to the number of cells added in the current inspection
                }
                // Third neighbor
                if (ptr_node->ptr_box[box_idxNbr_y_plus] == zone_idx)
                {
                    ptr_node->ptr_aux_idx[subzone_size + cntr_cell_add] = box_idxNbr_y_plus; // Including the neighboring element of the box to the auxiliary array
                    ptr_node->ptr_box[box_idxNbr_y_plus] = -1;                               // Changing the box status from refinement zone ID (>= 0) to the REFINEMENT REQUIRED (-1)
                    cntr_cell_add++;                                                      // +1 to the number of cells added in the current inspection
                }
                // Fourth neighbor
                if (ptr_node->ptr_box[box_idxNbr_y_minus] == zone_idx)
                {
                    ptr_node->ptr_aux_idx[subzone_size + cntr_cell_add] = box_idxNbr_y_minus; // Including the neighboring element of the box to the auxiliary array
                    ptr_node->ptr_box[box_idxNbr_y_minus] = -1;                               // Changing the box status from refinement zone ID (>= 0) to the REFINEMENT REQUIRED (-1)
                    cntr_cell_add++;                                                       // +1 to the number of cells added in the current inspection
                }
                // Fifth neighbor
                if (ptr_node->ptr_box[box_idxNbr_z_plus] == zone_idx)
                {
                    ptr_node->ptr_aux_idx[subzone_size + cntr_cell_add] = box_idxNbr_z_plus; // Including the neighboring element of the box to the auxiliary array
                    ptr_node->ptr_box[box_idxNbr_z_plus] = -1;                               // Changing the box status from refinement zone ID (>= 0) to the REFINEMENT REQUIRED (-1)
                    cntr_cell_add++;                                                      // +1 to the number of cells added in the current inspection
                }
                // Sixth neighbor
                if (ptr_node->ptr_box[box_idxNbr_z_minus] == zone_idx)
                {
                    ptr_node->ptr_aux_idx[subzone_size + cntr_cell_add] = box_idxNbr_z_minus; // Including the neighboring element of the box to the auxiliary array
                    ptr_node->ptr_box[box_idxNbr_z_minus] = -1;                               // Changing the box status from refinement zone ID (>= 0) to the REFINEMENT REQUIRED (-1)
                    cntr_cell_add++;                                                       // +1 to the number of cells added in the current inspection
                }
                subzone_size += cntr_cell_add;            // Increasing the number of cells in the subzone
                cntr_cell_add_all_subzones += cntr_cell_add; // Increasing the number of cells added to the subzones
                cntr_insp++;                              // Increasing the number of inspections
            }                                             // End while cycle, now the box contains the the information about all cells of the subzone "subzone_idx"

            //** >> Space checking of refinement zones arrays, refinement capacity array and refinement size array **/
            if (subzone_idx_max < subzone_idx + 1)
            {
                if (space_check(&(subzone_idx_max), subzone_idx + 1, 2.0f, "p3i2i1i1", &(ptr_node->pptr_subzones), &(ptr_node->ptr_subzone_cap), &(ptr_node->ptr_subzone_size)) == _FAILURE_)
                {
                    printf("Error, in space_check function\n");
                    return _FAILURE_;
                }
            }
            ptr_node->ptr_subzone_size[subzone_idx] = subzone_size;
            subzone_idx++; // Increasing the subzone number
        }               // subZone defined in the box
        cell_ref_idx++;
    } // At this point the box contains the information of all refinement subzones

    ptr_node->subzones_cap = subzone_idx_max; // Maximum amount of subzones
    ptr_node->subzones_size = subzone_idx;    // Total amount of subzones

    //** >> Space checking in each subzone **/
    for (int i = 0; i < ptr_node->subzones_size; i++)
    {
        if (space_check(&(ptr_node->ptr_subzone_cap[i]), ptr_node->ptr_subzone_size[i], 2.0f, "p1i1", &(ptr_node->pptr_subzones[i])) == _FAILURE_)
        {
            printf("Error, in space_check function\n");
            return _FAILURE_;
        }
    }

    //** >> Initializing the ptr_aux_idx array to be used as a counter of elemnents in each subzone**/
    // Notes that the number of refined cells is always bigger or equal than the number of refined subzones, so that the capacity of ptr_aux_idx is always enough to counter the number of refined subzones
    for (int i = 0; i < ptr_node->subzones_size; i++)
    {
        ptr_node->ptr_aux_idx[i] = 0;
    }

    //** >> Adding the cells to the zone array pptr_zones **/
    for (int i = 0; i < ptr_node->ptr_zone_size[zone_idx]; i++)
    {
        cell_idx = ptr_node->pptr_zones[zone_idx][i];
        box_idx_node = ptr_node->ptr_box_idx[cell_idx];
        subzone_idx = ptr_node->ptr_box[box_idx_node];
        cntr_cell_add = ptr_node->ptr_aux_idx[subzone_idx];          // Counter the element in the subzone "subzone_idx"
        ptr_node->pptr_subzones[subzone_idx][cntr_cell_add] = cell_idx; // Adding the index of the cell array in the block to the subzone
        ptr_node->ptr_aux_idx[subzone_idx] += 1;                     // Counter the number of elements added in the subzone "subzone_idx"
        
        ptr_node->ptr_box[box_idx_node] = zone_idx; // Returning the value zone_idx at the box in the zone partitioned 
    }

    return _SUCCESS_;
}

static int fill_child_nodes_PERIODIC_BOUNDARY(struct node *ptr_node)
{

    struct node *ptr_ch; // Pointer to the child node

    int cell_idx;  // Cell index
    int aux_idx_x; // Auxiliary indexes
    int aux_idx_y;
    int aux_idx_z;

    int box_idx_x_ch;  // Box index in X direcction of child node
    int box_idx_y_ch;  // Box index in Y direcction of child node
    int box_idx_z_ch;  // Box index in Z direcction of child node
    int box_idx_ch;    // Box index
    int box_idxNbr_ch; // Box index of the neighboring child cell

    int box_idx_node; // Box index

    int box_idxNbr_i0_j0_k0_ch;    // Box index at x=0,y=0,z=0
    int box_idxNbr_im1_j0_k0_ch;   // Box index of the neighbor at x=-1,y=0,z=0
    int box_idxNbr_i0_jm1_k0_ch;   // Box index of the neighbor at x=0,y=-1,z=0
    int box_idxNbr_im1_jm1_k0_ch;  // Box index of the neighbor at x=-1,y=-1,z=0
    int box_idxNbr_i0_j0_km1_ch;   // Box index of the neighbor at x=0,y=0,z=-1
    int box_idxNbr_im1_j0_km1_ch;  // Box index of the neighbor at x=-1,y=0,z=-1
    int box_idxNbr_i0_jm1_km1_ch;  // Box index of the neighbor at x=0,y=-1,z=-1
    int box_idxNbr_im1_jm1_km1_ch; // Box index of the neighbor at x=-1,y=-1,z=-1

    int box_idx_ptcl_ch;

    int box_grid_idx_ch; // Box grid index

    bool is_bder_grid_point; // Ask if the grid point is interior
    bool grid_point_exist;
    bool grid_point_simulation_boundary;

    int ptcl_idx;

    int pos_x; // Distance between the real box and the minimal box when the last is localized in the middle of the real one
    int pos_y;
    int pos_z;

    int cap;  // Capacity of arrays in child nodes
    int size; // Size or number of elements in some array in child nodes

    int lv = ptr_node->lv;

    int box_real_dim_X_ch;
    int box_real_dim_X_times_Y_ch;

    int grid_box_real_dim_X_ch;
    int grid_box_real_dim_X_times_Y_ch;

    int size_bder_grid_points_ch;
    int size_intr_grid_points_ch;

    //** >> Creating child nodes in parent node **/
    ptr_node->chn_cap = ptr_node->zones_cap; // Same amount than refinement zones
    ptr_node->chn_size = ptr_node->zones_size;

    if (ptr_node->chn_cap > 0)
    {
        ptr_node->pptr_chn = (struct node **)malloc(ptr_node->chn_cap * sizeof(struct node *)); // Allocating children pointers
    }

    for (int i = 0; i < ptr_node->chn_cap; i++)
    {
        ptr_node->pptr_chn[i] = NULL;
    }

    //** >> Filling child nodes **/
    for (int zone_idx = 0; zone_idx < ptr_node->zones_size; zone_idx++)
    {
        //** >> Create and Initialize child node **/
        ptr_ch = (struct node *)malloc(sizeof(struct node));
        initialize_node(ptr_ch);

        //** >> Global properties **/
        ptr_ch->ID = zone_idx;
        ptr_ch->lv = lv + 1;

        //** >> Boxes **/
        // MIN and MAX cell indexes values of the node.

        ptr_ch->boundary_simulation_contact_x = ptr_node->ptr_aux_bool_boundary_simulation_contact_x[zone_idx];
        ptr_ch->boundary_simulation_contact_y = ptr_node->ptr_aux_bool_boundary_simulation_contact_y[zone_idx];
        ptr_ch->boundary_simulation_contact_z = ptr_node->ptr_aux_bool_boundary_simulation_contact_z[zone_idx];

        ptr_ch->boundary_simulation_contact = ptr_ch->boundary_simulation_contact_x ||
                                              ptr_ch->boundary_simulation_contact_y ||
                                              ptr_ch->boundary_simulation_contact_z;

        ptr_ch->pbc_crosses_the_boundary_simulation_box_x = ptr_node->ptr_aux_bool_boundary_anomalies_x[zone_idx]; // The node has crossed the simulation box at X direction
        ptr_ch->pbc_crosses_the_boundary_simulation_box_y = ptr_node->ptr_aux_bool_boundary_anomalies_y[zone_idx];
        ptr_ch->pbc_crosses_the_boundary_simulation_box_z = ptr_node->ptr_aux_bool_boundary_anomalies_z[zone_idx];
        ptr_ch->pbc_anomalies_due_to_the_boundary = ptr_ch->pbc_crosses_the_boundary_simulation_box_x ||
                                                ptr_ch->pbc_crosses_the_boundary_simulation_box_y ||
                                                ptr_ch->pbc_crosses_the_boundary_simulation_box_z;

        //** >> Filling the different Sub zones of refinement **/
        if (ptr_ch->pbc_anomalies_due_to_the_boundary == true)
        {
            if (fill_subzones_ref_PERIODIC_BOUNDARY(ptr_node, zone_idx) == _FAILURE_)
            {
                printf("Error at function fill_subzones_ref_PERIODIC_BOUNDARY()\n");
                return _FAILURE_;
            }
        }

        if (ptr_ch->pbc_crosses_the_boundary_simulation_box_x == false)
        {
            for (int j = 0; j < ptr_node->ptr_zone_size[zone_idx]; j++)
            {
                cell_idx = ptr_node->pptr_zones[zone_idx][j]; // Cell index in the parent zone
                if (ptr_ch->box_min_x > ptr_node->ptr_cell_idx_x[cell_idx])
                {
                    ptr_ch->box_min_x = ptr_node->ptr_cell_idx_x[cell_idx];
                }
                if (ptr_ch->box_max_x < ptr_node->ptr_cell_idx_x[cell_idx])
                {
                    ptr_ch->box_max_x = ptr_node->ptr_cell_idx_x[cell_idx];
                }
            }
        }
        else
        {
            for (int subzone_idx = 0; subzone_idx < ptr_node->subzones_size; subzone_idx++)
            {
                ptr_node->ptr_aux_idx[0] = INT_MAX; // Minimum value
                ptr_node->ptr_aux_idx[1] = INT_MIN; // Maximum value

                for (int j = 0; j < ptr_node->ptr_subzone_size[subzone_idx]; j++)
                {
                    cell_idx = ptr_node->pptr_subzones[subzone_idx][j]; // Cell index in the parent zone
                    if (ptr_node->ptr_aux_idx[0] > ptr_node->ptr_cell_idx_x[cell_idx])
                    {
                        ptr_node->ptr_aux_idx[0] = ptr_node->ptr_cell_idx_x[cell_idx];
                    }
                    if (ptr_node->ptr_aux_idx[1] < ptr_node->ptr_cell_idx_x[cell_idx])
                    {
                        ptr_node->ptr_aux_idx[1] = ptr_node->ptr_cell_idx_x[cell_idx];
                    }
                }

                if (ptr_node->ptr_aux_idx[0] != 0)
                {
                    ptr_node->ptr_aux_idx[0] -= (1 << lv);
                    ptr_node->ptr_aux_idx[1] -= (1 << lv);
                }

                if (ptr_ch->box_min_x > ptr_node->ptr_aux_idx[0])
                {
                    ptr_ch->box_min_x = ptr_node->ptr_aux_idx[0];
                }
                if (ptr_ch->box_max_x < ptr_node->ptr_aux_idx[1])
                {
                    ptr_ch->box_max_x = ptr_node->ptr_aux_idx[1];
                }
            }
            if (ptr_ch->box_max_x - ptr_ch->box_min_x + 1 >= (1 << lv))
            {
                ptr_ch->pbc_crosses_the_whole_simulation_box_x = true;
                //ptr_ch->pbc_crosses_the_boundary_simulation_box_x = false; // Removing the state of crossing the boundary simulation box
                ptr_ch->box_min_x = 0;
                ptr_ch->box_max_x = (1 << lv) - 1;
            }
        }

        if (ptr_node->pbc_crosses_the_boundary_simulation_box_y == false)
        {
            for (int j = 0; j < ptr_node->ptr_zone_size[zone_idx]; j++)
            {
                cell_idx = ptr_node->pptr_zones[zone_idx][j]; // Cell index in the parent zone
                if (ptr_ch->box_min_y > ptr_node->ptr_cell_idx_y[cell_idx])
                {
                    ptr_ch->box_min_y = ptr_node->ptr_cell_idx_y[cell_idx];
                }
                if (ptr_ch->box_max_y < ptr_node->ptr_cell_idx_y[cell_idx])
                {
                    ptr_ch->box_max_y = ptr_node->ptr_cell_idx_y[cell_idx];
                }
            }
        }
        else
        {
            for (int subzone_idx = 0; subzone_idx < ptr_node->subzones_size; subzone_idx++)
            {
                ptr_node->ptr_aux_idx[0] = INT_MAX; // Minimum value
                ptr_node->ptr_aux_idx[1] = INT_MIN; // Maximum value

                for (int j = 0; j < ptr_node->ptr_subzone_size[subzone_idx]; j++)
                {
                    cell_idx = ptr_node->pptr_subzones[subzone_idx][j]; // Cell index in the parent zone
                    if (ptr_node->ptr_aux_idx[0] > ptr_node->ptr_cell_idx_y[cell_idx])
                    {
                        ptr_node->ptr_aux_idx[0] = ptr_node->ptr_cell_idx_y[cell_idx];
                    }
                    if (ptr_node->ptr_aux_idx[1] < ptr_node->ptr_cell_idx_y[cell_idx])
                    {
                        ptr_node->ptr_aux_idx[1] = ptr_node->ptr_cell_idx_y[cell_idx];
                    }
                }

                if (ptr_node->ptr_aux_idx[0] != 0)
                {
                    ptr_node->ptr_aux_idx[0] -= (1 << lv);
                    ptr_node->ptr_aux_idx[1] -= (1 << lv);
                }

                if (ptr_ch->box_min_y > ptr_node->ptr_aux_idx[0])
                {
                    ptr_ch->box_min_y = ptr_node->ptr_aux_idx[0];
                }
                if (ptr_ch->box_max_y < ptr_node->ptr_aux_idx[1])
                {
                    ptr_ch->box_max_y = ptr_node->ptr_aux_idx[1];
                }
            }
            if (ptr_ch->box_max_y - ptr_ch->box_min_y + 1 >= (1 << lv))
            {
                ptr_ch->pbc_crosses_the_whole_simulation_box_y = true;
                //ptr_ch->pbc_crosses_the_boundary_simulation_box_y = false; // Removing the state of crossing the boundary simulation box
                ptr_ch->box_min_y = 0;
                ptr_ch->box_max_y = (1 << lv) - 1;
            }
        }

        if (ptr_node->pbc_crosses_the_boundary_simulation_box_z == false)
        {
            for (int j = 0; j < ptr_node->ptr_zone_size[zone_idx]; j++)
            {
                cell_idx = ptr_node->pptr_zones[zone_idx][j]; // Cell index in the parent zone
                if (ptr_ch->box_min_z > ptr_node->ptr_cell_idx_z[cell_idx])
                {
                    ptr_ch->box_min_z = ptr_node->ptr_cell_idx_z[cell_idx];
                }
                if (ptr_ch->box_max_z < ptr_node->ptr_cell_idx_z[cell_idx])
                {
                    ptr_ch->box_max_z = ptr_node->ptr_cell_idx_z[cell_idx];
                }
            }
        }
        else
        {
            for (int subzone_idx = 0; subzone_idx < ptr_node->subzones_size; subzone_idx++)
            {
                ptr_node->ptr_aux_idx[0] = INT_MAX; // Minimum value
                ptr_node->ptr_aux_idx[1] = INT_MIN; // Maximum value

                for (int j = 0; j < ptr_node->ptr_subzone_size[subzone_idx]; j++)
                {
                    cell_idx = ptr_node->pptr_subzones[subzone_idx][j]; // Cell index in the parent zone
                    if (ptr_node->ptr_aux_idx[0] > ptr_node->ptr_cell_idx_z[cell_idx])
                    {
                        ptr_node->ptr_aux_idx[0] = ptr_node->ptr_cell_idx_z[cell_idx];
                    }
                    if (ptr_node->ptr_aux_idx[1] < ptr_node->ptr_cell_idx_z[cell_idx])
                    {
                        ptr_node->ptr_aux_idx[1] = ptr_node->ptr_cell_idx_z[cell_idx];
                    }
                }

                if (ptr_node->ptr_aux_idx[0] != 0)
                {
                    ptr_node->ptr_aux_idx[0] -= (1 << lv);
                    ptr_node->ptr_aux_idx[1] -= (1 << lv);
                }

                if (ptr_ch->box_min_z > ptr_node->ptr_aux_idx[0])
                {
                    ptr_ch->box_min_z = ptr_node->ptr_aux_idx[0];
                }
                if (ptr_ch->box_max_z < ptr_node->ptr_aux_idx[1])
                {
                    ptr_ch->box_max_z = ptr_node->ptr_aux_idx[1];
                }
            }
            if (ptr_ch->box_max_z - ptr_ch->box_min_z + 1 >= (1 << lv))
            {
                ptr_ch->pbc_crosses_the_whole_simulation_box_z = true;
                //ptr_ch->pbc_crosses_the_boundary_simulation_box_z = false; // Removing the state of crossing the boundary simulation box
                ptr_ch->box_min_z = 0;
                ptr_ch->box_max_z = (1 << lv) - 1;
            }
        }



        // Changing the min and max of the "minimal box" from parent units to child units
        ptr_ch->box_min_x = 2 * ptr_ch->box_min_x;
        ptr_ch->box_min_y = 2 * ptr_ch->box_min_y;
        ptr_ch->box_min_z = 2 * ptr_ch->box_min_z;
        ptr_ch->box_max_x = 2 * ptr_ch->box_max_x + 1;
        ptr_ch->box_max_y = 2 * ptr_ch->box_max_y + 1;
        ptr_ch->box_max_z = 2 * ptr_ch->box_max_z + 1;

        // Size of the "minimal box"
        ptr_ch->box_dim_x = ptr_ch->box_max_x - ptr_ch->box_min_x + 1;
        ptr_ch->box_dim_y = ptr_ch->box_max_y - ptr_ch->box_min_y + 1;
        ptr_ch->box_dim_z = ptr_ch->box_max_z - ptr_ch->box_min_z + 1;

        // Real dimensions of the boxcap = ptr_ch->box_cap;
        ptr_ch->box_real_dim_x = 5 > (n_exp - 1) ? (ptr_ch->box_dim_x + 10) : (ptr_ch->box_dim_x + 2 * n_exp - 2);
        ptr_ch->box_real_dim_y = 5 > (n_exp - 1) ? (ptr_ch->box_dim_y + 10) : (ptr_ch->box_dim_y + 2 * n_exp - 2);
        ptr_ch->box_real_dim_z = 5 > (n_exp - 1) ? (ptr_ch->box_dim_z + 10) : (ptr_ch->box_dim_z + 2 * n_exp - 2);

        ptr_ch->box_real_dim_x_old = ptr_ch->box_real_dim_x;
        ptr_ch->box_real_dim_y_old = ptr_ch->box_real_dim_y;
        // ptr_ch->box_real_dim_z_old = ptr_ch->box_real_dim_z;

        // Translations between cell array and box
        pos_x = (ptr_ch->box_real_dim_x - ptr_ch->box_dim_x) / 2; // Half of the distance of the box side less the "minimal box" side
        ptr_ch->box_ts_x = ptr_ch->box_min_x - pos_x;             // Every cell in the level l in the box must be subtracted this value to obtain the box index
        pos_y = (ptr_ch->box_real_dim_y - ptr_ch->box_dim_y) / 2;
        ptr_ch->box_ts_y = ptr_ch->box_min_y - pos_y;
        pos_z = (ptr_ch->box_real_dim_z - ptr_ch->box_dim_z) / 2;
        ptr_ch->box_ts_z = ptr_ch->box_min_z - pos_z;

        ptr_ch->box_ts_x_old = ptr_ch->box_ts_x;
        ptr_ch->box_ts_y_old = ptr_ch->box_ts_y;
        ptr_ch->box_ts_z_old = ptr_ch->box_ts_z;

        //** >> Cells in the node **/
        cap = 2 * 8 * ptr_node->ptr_zone_size[zone_idx]; // Capacity is the double of total child cells
        ptr_ch->cell_cap = cap;
        size = 8 * ptr_node->ptr_zone_size[zone_idx];
        ptr_ch->cell_size = size;
        ptr_ch->ptr_cell_idx_x = (int *)malloc(ptr_ch->cell_cap * sizeof(int));
        ptr_ch->ptr_cell_idx_y = (int *)malloc(ptr_ch->cell_cap * sizeof(int));
        ptr_ch->ptr_cell_idx_z = (int *)malloc(ptr_ch->cell_cap * sizeof(int));
        ptr_ch->ptr_box_idx = (int *)malloc(ptr_ch->cell_cap * sizeof(int));

        box_real_dim_X_ch = ptr_ch->box_real_dim_x;
        box_real_dim_X_times_Y_ch = ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;
        if (ptr_ch->pbc_anomalies_due_to_the_boundary == false)
        {
            for (int j = 0; j < ptr_node->ptr_zone_size[zone_idx]; j++)
            {
                aux_idx_x = ptr_node->ptr_cell_idx_x[ptr_node->pptr_zones[zone_idx][j]] * 2;
                aux_idx_y = ptr_node->ptr_cell_idx_y[ptr_node->pptr_zones[zone_idx][j]] * 2;
                aux_idx_z = ptr_node->ptr_cell_idx_z[ptr_node->pptr_zones[zone_idx][j]] * 2;

                box_idx_ch = (aux_idx_x - ptr_ch->box_ts_x) + (aux_idx_y - ptr_ch->box_ts_y) * box_real_dim_X_ch + (aux_idx_z - ptr_ch->box_ts_z) * box_real_dim_X_times_Y_ch;

                for (int kk = 0; kk < 2; kk++)
                {
                    for (int jj = 0; jj < 2; jj++)
                    {
                        for (int ii = 0; ii < 2; ii++)
                        {
                            cell_idx = ii + jj * 2 + kk * 4 + j * 8; // Child cell index
                            ptr_ch->ptr_cell_idx_x[cell_idx] = aux_idx_x + ii;
                            ptr_ch->ptr_cell_idx_y[cell_idx] = aux_idx_y + jj;
                            ptr_ch->ptr_cell_idx_z[cell_idx] = aux_idx_z + kk;
                            ptr_ch->ptr_box_idx[cell_idx] = box_idx_ch + ii + jj * box_real_dim_X_ch + kk * box_real_dim_X_times_Y_ch;
                        }
                    }
                }
            }
        }
        else
        {
            for (int j = 0; j < ptr_node->ptr_zone_size[zone_idx]; j++)
            {
                aux_idx_x = ptr_node->ptr_cell_idx_x[ptr_node->pptr_zones[zone_idx][j]] * 2;
                aux_idx_y = ptr_node->ptr_cell_idx_y[ptr_node->pptr_zones[zone_idx][j]] * 2;
                aux_idx_z = ptr_node->ptr_cell_idx_z[ptr_node->pptr_zones[zone_idx][j]] * 2;

                box_idx_ch = (aux_idx_x - ptr_ch->box_ts_x) + (aux_idx_y - ptr_ch->box_ts_y) * box_real_dim_X_ch + (aux_idx_z - ptr_ch->box_ts_z) * box_real_dim_X_times_Y_ch;

                if (ptr_ch->pbc_crosses_the_boundary_simulation_box_x == true && ptr_ch->pbc_crosses_the_whole_simulation_box_x == false && aux_idx_x > ptr_ch->box_max_x)
                {
                    box_idx_ch -= (1 << (lv + 1));
                }
                if (ptr_ch->pbc_crosses_the_boundary_simulation_box_y == true && ptr_ch->pbc_crosses_the_whole_simulation_box_y == false && aux_idx_y > ptr_ch->box_max_y)
                {
                    box_idx_ch -= (1 << (lv + 1)) * box_real_dim_X_ch;
                }
                if (ptr_ch->pbc_crosses_the_boundary_simulation_box_z == true && ptr_ch->pbc_crosses_the_whole_simulation_box_z == false && aux_idx_z > ptr_ch->box_max_z)
                {
                    box_idx_ch -= (1 << (lv + 1)) * box_real_dim_X_times_Y_ch;
                }

                for (int kk = 0; kk < 2; kk++)
                {
                    for (int jj = 0; jj < 2; jj++)
                    {
                        for (int ii = 0; ii < 2; ii++)
                        {
                            cell_idx = ii + jj * 2 + kk * 4 + j * 8; // Child cell index
                            ptr_ch->ptr_cell_idx_x[cell_idx] = aux_idx_x + ii;
                            ptr_ch->ptr_cell_idx_y[cell_idx] = aux_idx_y + jj;
                            ptr_ch->ptr_cell_idx_z[cell_idx] = aux_idx_z + kk;
                            ptr_ch->ptr_box_idx[cell_idx] = box_idx_ch + ii + jj * box_real_dim_X_ch + kk * box_real_dim_X_times_Y_ch;
                        }
                    }
                }
            }
        }


        // Filling the box status
        cap = ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y * ptr_ch->box_real_dim_z; // In general, the size of each side must be 3 times bigger than the same side of the "minimal box"
        ptr_ch->box_cap = cap;
        ptr_ch->ptr_box = (int *)malloc(cap * sizeof(int));
        ptr_ch->ptr_box_old = (int *)malloc(cap * sizeof(int));
        // Putting the value of NO-EXIST (-4) in every box index
        for (int j = 0; j < cap; j++)
        {
            ptr_ch->ptr_box[j] = -4;
        }
        // Changing from NO-EXIST (-4) to EXIST (-3) to all cells in the block of the child node
        for (int j = 0; j < ptr_ch->cell_size; j++)
        {
            box_idx_ch = ptr_ch->ptr_box_idx[j];
            ptr_ch->ptr_box[box_idx_ch] = -3;
        }
        //Adding -5 to the boundary when it corresponds
        if (ptr_ch->boundary_simulation_contact == true)
        {
            if (ptr_ch->pbc_crosses_the_whole_simulation_box_x == true || boundary_type != 0)
            {
                for (int k = 0; k < ptr_ch->box_real_dim_z;k++)
                {
                    for (int j = 0; j < ptr_ch->box_real_dim_y; j++)
                    {
                        for (int i = 0; i < (ptr_ch->box_real_dim_x - ptr_ch->box_dim_x)/2; i++)
                        {
                            box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
                            ptr_ch->ptr_box[box_idx_ch] = -5;
                        }

                        for (int i = (ptr_ch->box_real_dim_x + ptr_ch->box_dim_x) / 2; i < ptr_ch->box_real_dim_x; i++)
                        {
                            box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
                            ptr_ch->ptr_box[box_idx_ch] = -5;
                        }
                    }
                }
            }

            if (ptr_ch->pbc_crosses_the_whole_simulation_box_y == true || boundary_type != 0)
            {
                for (int k = 0; k < ptr_ch->box_real_dim_z; k++)
                {
                    for (int i = 0; i < ptr_ch->box_real_dim_x; i++)
                    {
                        for (int j = 0; j < (ptr_ch->box_real_dim_y - ptr_ch->box_dim_y) / 2; j++)
                        {
                            box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
                            ptr_ch->ptr_box[box_idx_ch] = -5;
                        }

                        for (int j = (ptr_ch->box_real_dim_y + ptr_ch->box_dim_y) / 2; j < ptr_ch->box_real_dim_y; j++)
                        {
                            box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
                            ptr_ch->ptr_box[box_idx_ch] = -5;
                        }
                    }
                }
            }

            if (ptr_ch->pbc_crosses_the_whole_simulation_box_z == true)
            {
                for (int j = 0; j < ptr_ch->box_real_dim_y; j++)
                {
                    for (int i = 0; i < ptr_ch->box_real_dim_x; i++)
                    {
                        for (int k = 0; k < (ptr_ch->box_real_dim_z - ptr_ch->box_dim_z) / 2; k++)
                        {
                            box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
                            ptr_ch->ptr_box[box_idx_ch] = -5;
                        }

                        for (int k = (ptr_ch->box_real_dim_z + ptr_ch->box_dim_z) / 2; k < ptr_ch->box_real_dim_z; k++)
                        {
                            box_idx_ch = i + j * box_real_dim_X_ch + k * box_real_dim_X_times_Y_ch;
                            ptr_ch->ptr_box[box_idx_ch] = -5;
                        }
                    }
                }
            }
        }

        //** >> Struct of cells (Particles and cell mass)
        //** >> Total mass in the node **/

        ptr_ch->ptr_cell_struct = (struct cell_struct *)malloc(ptr_ch->box_cap * sizeof(struct cell_struct));
        // ptr_ch->ptr_cell_struct_old = (struct cell_struct *)malloc(ptr_ch->box_cap * sizeof(struct cell_struct));
        for (int j = 0; j < ptr_ch->box_cap; j++)
        {
            initialize_cell_struct(&(ptr_ch->ptr_cell_struct[j]));
            // initialize_cell_struct(&(ptr_ch->ptr_cell_struct_old[j]));
        }

        if (ptr_node->pbc_anomalies_due_to_the_boundary == false)
        {
            for (int j = 0; j < ptr_node->ptr_zone_size[zone_idx]; j++)
            {
                box_idx_node = ptr_node->ptr_box_idx[ptr_node->pptr_zones[zone_idx][j]];
                // box_idx_x_node = ptr_node->ptr_cell_idx_x[ptr_node->pptr_zones[i][j]] - ptr_node->box_ts_x;
                // box_idx_y_node = ptr_node->ptr_cell_idx_y[ptr_node->pptr_zones[i][j]] - ptr_node->box_ts_y;
                // box_idx_z_node = ptr_node->ptr_cell_idx_z[ptr_node->pptr_zones[i][j]] - ptr_node->box_ts_z;
                // box_idx_node = box_idx_x_node + box_idx_y_node * ptr_node->box_real_dim_x + box_idx_z_node * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;
                if (ptr_node->ptr_cell_struct[box_idx_node].ptcl_size > 0)
                {

                    // Allocating space fot the child particles in the 8  cells
                    box_idx_x_ch = ptr_node->ptr_cell_idx_x[ptr_node->pptr_zones[zone_idx][j]] * 2 - ptr_ch->box_ts_x;
                    box_idx_y_ch = ptr_node->ptr_cell_idx_y[ptr_node->pptr_zones[zone_idx][j]] * 2 - ptr_ch->box_ts_y;
                    box_idx_z_ch = ptr_node->ptr_cell_idx_z[ptr_node->pptr_zones[zone_idx][j]] * 2 - ptr_ch->box_ts_z;
                    box_idx_ch = box_idx_x_ch + box_idx_y_ch * box_real_dim_X_ch + box_idx_z_ch * box_real_dim_X_times_Y_ch;

                    for (int kk = 0; kk < 2; kk++)
                    {
                        for (int jj = 0; jj < 2; jj++)
                        {
                            for (int ii = 0; ii < 2; ii++)
                            {
                                box_idxNbr_ch = box_idx_ch + ii + jj * box_real_dim_X_ch + kk * box_real_dim_X_times_Y_ch;
                                cap = ptr_node->ptr_cell_struct[box_idx_node].ptcl_size;
                                ptr_ch->ptr_cell_struct[box_idxNbr_ch].ptcl_cap = cap;
                                ptr_ch->ptr_cell_struct[box_idxNbr_ch].ptr_ptcl = (int *)malloc(cap * sizeof(int));
                            }
                        }
                    }

                    // Putting parent cell particles into the child cell
                    for (int k = 0; k < ptr_node->ptr_cell_struct[box_idx_node].ptcl_size; k++)
                    {
                        ptcl_idx = ptr_node->ptr_cell_struct[box_idx_node].ptr_ptcl[k];

                        box_idx_ptcl_ch = ptcl_idx_to_box_idx(ptr_ch, ptcl_idx);

                        ptr_ch->ptr_cell_struct[box_idx_ptcl_ch].ptr_ptcl[ptr_ch->ptr_cell_struct[box_idx_ptcl_ch].ptcl_size] = ptcl_idx;
                        ptr_ch->ptr_cell_struct[box_idx_ptcl_ch].ptcl_size += 1;
                        ptr_ch->ptr_cell_struct[box_idx_ptcl_ch].cell_mass += GL_ptcl_mass[ptcl_idx];
                    }

                    ptr_ch->local_mass += ptr_node->ptr_cell_struct[box_idx_node].cell_mass; // Local mass
                }
            }
        }
        else
        {
            for (int j = 0; j < ptr_node->ptr_zone_size[zone_idx]; j++)
            {
                box_idx_node = ptr_node->ptr_box_idx[ptr_node->pptr_zones[zone_idx][j]];
                // box_idx_x_node = ptr_node->ptr_cell_idx_x[ptr_node->pptr_zones[i][j]] - ptr_node->box_ts_x;
                // box_idx_y_node = ptr_node->ptr_cell_idx_y[ptr_node->pptr_zones[i][j]] - ptr_node->box_ts_y;
                // box_idx_z_node = ptr_node->ptr_cell_idx_z[ptr_node->pptr_zones[i][j]] - ptr_node->box_ts_z;
                // box_idx_node = box_idx_x_node + box_idx_y_node * ptr_node->box_real_dim_x + box_idx_z_node * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;
                if (ptr_node->ptr_cell_struct[box_idx_node].ptcl_size > 0)
                {
                    // Allocating space fot the child particles in the 8  cells
                    aux_idx_x = ptr_node->ptr_cell_idx_x[ptr_node->pptr_zones[zone_idx][j]] * 2;
                    aux_idx_y = ptr_node->ptr_cell_idx_y[ptr_node->pptr_zones[zone_idx][j]] * 2;
                    aux_idx_z = ptr_node->ptr_cell_idx_z[ptr_node->pptr_zones[zone_idx][j]] * 2;

                    box_idx_ch = (aux_idx_x - ptr_ch->box_ts_x) + (aux_idx_y - ptr_ch->box_ts_y) * box_real_dim_X_ch + (aux_idx_z - ptr_ch->box_ts_z) * box_real_dim_X_times_Y_ch;

                    if (ptr_ch->pbc_crosses_the_boundary_simulation_box_x == true && ptr_ch->pbc_crosses_the_whole_simulation_box_x == false && aux_idx_x > ptr_ch->box_max_x)
                    {
                        box_idx_ch -= (1 << (lv + 1));
                    }
                    if (ptr_ch->pbc_crosses_the_boundary_simulation_box_y == true && ptr_ch->pbc_crosses_the_whole_simulation_box_y == false && aux_idx_y > ptr_ch->box_max_y)
                    {
                        box_idx_ch -= (1 << (lv + 1)) * box_real_dim_X_ch;
                    }
                    if (ptr_ch->pbc_crosses_the_boundary_simulation_box_z == true && ptr_ch->pbc_crosses_the_whole_simulation_box_z == false && aux_idx_z > ptr_ch->box_max_z)
                    {
                        box_idx_ch -= (1 << (lv + 1)) * box_real_dim_X_times_Y_ch;
                    }

                    for (int kk = 0; kk < 2; kk++)
                    {
                        for (int jj = 0; jj < 2; jj++)
                        {
                            for (int ii = 0; ii < 2; ii++)
                            {
                                box_idxNbr_ch = box_idx_ch + ii + jj * box_real_dim_X_ch + kk * box_real_dim_X_times_Y_ch;
                                cap = ptr_node->ptr_cell_struct[box_idx_node].ptcl_size;
                                ptr_ch->ptr_cell_struct[box_idxNbr_ch].ptcl_cap = cap;
                                ptr_ch->ptr_cell_struct[box_idxNbr_ch].ptr_ptcl = (int *)malloc(cap * sizeof(int));
                            }
                        }
                    }

                    // Putting parent cell particles into the child cell
                    for (int k = 0; k < ptr_node->ptr_cell_struct[box_idx_node].ptcl_size; k++)
                    {
                        ptcl_idx = ptr_node->ptr_cell_struct[box_idx_node].ptr_ptcl[k];

                        box_idx_ptcl_ch = ptcl_idx_to_box_idx_PERIODIC_BOUNDARY(ptr_ch, ptcl_idx);

                        ptr_ch->ptr_cell_struct[box_idx_ptcl_ch].ptr_ptcl[ptr_ch->ptr_cell_struct[box_idx_ptcl_ch].ptcl_size] = ptcl_idx;
                        ptr_ch->ptr_cell_struct[box_idx_ptcl_ch].ptcl_size += 1;
                        ptr_ch->ptr_cell_struct[box_idx_ptcl_ch].cell_mass += GL_ptcl_mass[ptcl_idx];
                    }

                    ptr_ch->local_mass += ptr_node->ptr_cell_struct[box_idx_node].cell_mass; // Local mass
                }
            }
        }


        //size = ptr_ch->cell_size / 8 * 18 + 9; // 27 * N - 9 * (N-1)

        // 27 * N - 9 * (N-1) = 18N+9 >= Total of grid points; N = number of cubics with 8 cells and 27 grid points
        // Maximum of border grid points = 9*(N*2+1) - (N*2-1) = 16N+10 (Imagine a row of cubics)
        // Maximum of interior grid points = (2k-1)^3 < 8k^3, where k^3 = N
        size_bder_grid_points_ch = ptr_ch->cell_size / 8 * 16 + 10;
        size_intr_grid_points_ch = ptr_ch->cell_size;

        //** >> Space checking of border grid points array**/
        ptr_ch->grid_intr_cap = size_intr_grid_points_ch;
        ptr_ch->grid_bder_cap = size_bder_grid_points_ch;
        ptr_ch->ptr_intr_grid_cell_idx_x = (int *)malloc(ptr_ch->grid_intr_cap * sizeof(int));
        ptr_ch->ptr_intr_grid_cell_idx_y = (int *)malloc(ptr_ch->grid_intr_cap * sizeof(int));
        ptr_ch->ptr_intr_grid_cell_idx_z = (int *)malloc(ptr_ch->grid_intr_cap * sizeof(int));
        ptr_ch->ptr_intr_grid_idx = (int *)malloc(ptr_ch->grid_intr_cap * sizeof(int));
        ptr_ch->ptr_bder_grid_cell_idx_x = (int *)malloc(ptr_ch->grid_bder_cap * sizeof(int));
        ptr_ch->ptr_bder_grid_cell_idx_y = (int *)malloc(ptr_ch->grid_bder_cap * sizeof(int));
        ptr_ch->ptr_bder_grid_cell_idx_z = (int *)malloc(ptr_ch->grid_bder_cap * sizeof(int));
        ptr_ch->ptr_bder_grid_idx = (int *)malloc(ptr_ch->grid_bder_cap * sizeof(int));

        if(ptr_ch->pbc_anomalies_due_to_the_boundary == true)
        {
            ptr_ch->grid_SIMULATION_BOUNDARY_cap = size_bder_grid_points_ch;

            ptr_ch->ptr_SIMULATION_BOUNDARY_grid_cell_idx_x = (int *)malloc(ptr_ch->grid_SIMULATION_BOUNDARY_cap * sizeof(int));
            ptr_ch->ptr_SIMULATION_BOUNDARY_grid_cell_idx_y = (int *)malloc(ptr_ch->grid_SIMULATION_BOUNDARY_cap * sizeof(int));
            ptr_ch->ptr_SIMULATION_BOUNDARY_grid_cell_idx_z = (int *)malloc(ptr_ch->grid_SIMULATION_BOUNDARY_cap * sizeof(int));
            ptr_ch->ptr_SIMULATION_BOUNDARY_grid_idx = (int *)malloc(ptr_ch->grid_SIMULATION_BOUNDARY_cap * sizeof(int));
        }

        //** >> Grid points **/
        grid_box_real_dim_X_ch = ptr_ch->box_real_dim_x + 1;
        grid_box_real_dim_X_times_Y_ch = (ptr_ch->box_real_dim_x + 1) * (ptr_ch->box_real_dim_y + 1);

        for (int kk = ptr_ch->box_min_z - ptr_ch->box_ts_z; kk < ptr_ch->box_max_z - ptr_ch->box_ts_z + 2; kk++)
        {
            for (int jj = ptr_ch->box_min_y - ptr_ch->box_ts_y; jj < ptr_ch->box_max_y - ptr_ch->box_ts_y + 2; jj++)
            {
                for (int ii = ptr_ch->box_min_x - ptr_ch->box_ts_x; ii < ptr_ch->box_max_x - ptr_ch->box_ts_x + 2; ii++)
                {
                    box_idx_ch = ii + jj * ptr_ch->box_real_dim_x + kk * ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;

                    grid_point_exist = false;
                    grid_point_simulation_boundary = false;

                    box_idxNbr_i0_j0_k0_ch = box_idx_ch;
                    box_idxNbr_im1_j0_k0_ch = box_idx_ch - 1;
                    box_idxNbr_i0_jm1_k0_ch = box_idx_ch - box_real_dim_X_ch;
                    box_idxNbr_im1_jm1_k0_ch = box_idx_ch - 1 - box_real_dim_X_ch;
                    box_idxNbr_i0_j0_km1_ch = box_idx_ch - box_real_dim_X_times_Y_ch;
                    box_idxNbr_im1_j0_km1_ch = box_idx_ch - 1 - box_real_dim_X_times_Y_ch;
                    box_idxNbr_i0_jm1_km1_ch = box_idx_ch - box_real_dim_X_ch - box_real_dim_X_times_Y_ch;
                    box_idxNbr_im1_jm1_km1_ch = box_idx_ch - 1 - box_real_dim_X_ch - box_real_dim_X_times_Y_ch;


                    // First
                    if (ptr_ch->ptr_box[box_idxNbr_i0_j0_k0_ch] > -4)
                    {
                        grid_point_exist = true;
                        is_bder_grid_point = false;
                        //** Connection to the right  **/

                        //** Connection to the left  **/
                        if (ptr_ch->ptr_box[box_idxNbr_im1_j0_k0_ch] < -3 &&
                            ptr_ch->ptr_box[box_idxNbr_im1_jm1_k0_ch] < -3 &&
                            ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] < -3 &&
                            ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] < -3)
                        {
                            is_bder_grid_point = true;
                        }
                        //** Backward connection   **/

                        //** Forward connection   **/
                        else if (ptr_ch->ptr_box[box_idxNbr_i0_jm1_k0_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_jm1_k0_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] < -3)
                        {
                            is_bder_grid_point = true;
                        }
                        //** Upward connection **/

                        //** Down connection **/
                        else if (ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] < -3)
                        {
                            is_bder_grid_point = true;
                        }
                    }
                    // Second
                    else if (ptr_ch->ptr_box[box_idxNbr_im1_j0_k0_ch] > -4)
                    {
                        grid_point_exist = true;
                        is_bder_grid_point = false;
                        //** Connection to the right  **/
                        if (ptr_ch->ptr_box[box_idxNbr_i0_jm1_k0_ch] < -3 &&
                            ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] < -3 &&
                            ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] < -3)
                        {
                            is_bder_grid_point = true;
                        }
                        //** Connection to the left  **/

                        //** Backward connection   **/

                        //** Forward connection   **/
                        else if (ptr_ch->ptr_box[box_idxNbr_i0_jm1_k0_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_jm1_k0_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] < -3)
                        {
                            is_bder_grid_point = true;
                        }
                        //** Upward connection **/

                        //** Down connection **/
                        else if (ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] < -3)
                        {
                            is_bder_grid_point = true;
                        }
                    }
                    // Third
                    else if (ptr_ch->ptr_box[box_idxNbr_i0_jm1_k0_ch] > -4)
                    {
                        grid_point_exist = true;
                        is_bder_grid_point = false;
                        //** Connection to the right  **/

                        //** Connection to the left  **/
                        if (ptr_ch->ptr_box[box_idxNbr_im1_jm1_k0_ch] < -3 &&
                            ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] < -3 &&
                            ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] < -3)
                        {
                            is_bder_grid_point = true;
                        }
                        //** Backward connection   **/
                        else if (ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] < -3)
                        {
                            is_bder_grid_point = true;
                        }
                        //** Forward connection   **/

                        //** Upward connection **/

                        //** Down connection **/
                        else if (ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] < -3)
                        {
                            is_bder_grid_point = true;
                        }
                    }
                    // Fourth
                    else if (ptr_ch->ptr_box[box_idxNbr_im1_jm1_k0_ch] > -4)
                    {
                        grid_point_exist = true;
                        is_bder_grid_point = false;
                        //** Connection to the right  **/
                        if (ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] < -3 &&
                            ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] < -3)
                        {
                            is_bder_grid_point = true;
                        }
                        //** Connection to the left  **/

                        //** Backward connection   **/
                        else if (ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] < -3)
                        {
                            is_bder_grid_point = true;
                        }
                        //** Forward connection   **/

                        //** Upward connection **/

                        //** Down connection **/
                        else if (ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] < -3)
                        {
                            is_bder_grid_point = true;
                        }
                    }
                    // Others
                    else if (ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] > -4 ||
                             ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] > -4 ||
                             ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] > -4 ||
                             ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] > -4)
                    {
                        grid_point_exist = true;
                        is_bder_grid_point = true;
                    }

                    if (ptr_ch->boundary_simulation_contact == true && grid_point_exist == true)
                    {
                        if (boundary_type == 0 || (boundary_type != 0 && is_bder_grid_point == true))
                        {
                            aux_idx_x = ii + ptr_ch->box_ts_x;
                            aux_idx_y = jj + ptr_ch->box_ts_y;
                            aux_idx_z = kk + ptr_ch->box_ts_z;

                            if (aux_idx_x == 0 || aux_idx_x == (1 << (lv + 1)) ||
                                aux_idx_y == 0 || aux_idx_y == (1 << (lv + 1)) ||
                                aux_idx_z == 0 || aux_idx_z == (1 << (lv + 1)))
                            {
                                grid_point_simulation_boundary = true;
                            }
                        }
                    }


                    if (grid_point_exist == true)
                    {
                        box_grid_idx_ch = ii + jj * grid_box_real_dim_X_ch + kk * grid_box_real_dim_X_times_Y_ch;

                        //** >> Adding the grid point **/

                        //** >> Simulation Boundary grid point**/
                        if (grid_point_simulation_boundary == true)
                        {
                            aux_idx_x = ii + ptr_ch->box_ts_x < 0 ? ii + ptr_ch->box_ts_x + (1 << (lv + 1)) : ii + ptr_ch->box_ts_x;
                            aux_idx_y = jj + ptr_ch->box_ts_y < 0 ? jj + ptr_ch->box_ts_y + (1 << (lv + 1)) : jj + ptr_ch->box_ts_y;
                            aux_idx_z = kk + ptr_ch->box_ts_z < 0 ? kk + ptr_ch->box_ts_z + (1 << (lv + 1)) : kk + ptr_ch->box_ts_z;
                            ptr_ch->ptr_SIMULATION_BOUNDARY_grid_cell_idx_x[ptr_ch->grid_SIMULATION_BOUNDARY_size] = aux_idx_x;
                            ptr_ch->ptr_SIMULATION_BOUNDARY_grid_cell_idx_y[ptr_ch->grid_SIMULATION_BOUNDARY_size] = aux_idx_y;
                            ptr_ch->ptr_SIMULATION_BOUNDARY_grid_cell_idx_z[ptr_ch->grid_SIMULATION_BOUNDARY_size] = aux_idx_z;
                            ptr_ch->ptr_SIMULATION_BOUNDARY_grid_idx[ptr_ch->grid_SIMULATION_BOUNDARY_size] = box_grid_idx_ch;
                            ptr_ch->grid_SIMULATION_BOUNDARY_size += 1; // Increasing the number of border grid points in the array
                        }
                        //** >> Border grid point**/
                        else if (is_bder_grid_point == true)
                        {
                            //** >> Adding the grid point to the border array **/
                            ptr_ch->ptr_bder_grid_cell_idx_x[ptr_ch->grid_bder_size] = ii + ptr_ch->box_ts_x;
                            ptr_ch->ptr_bder_grid_cell_idx_y[ptr_ch->grid_bder_size] = jj + ptr_ch->box_ts_y;
                            ptr_ch->ptr_bder_grid_cell_idx_z[ptr_ch->grid_bder_size] = kk + ptr_ch->box_ts_z;
                            ptr_ch->ptr_bder_grid_idx[ptr_ch->grid_bder_size] = box_grid_idx_ch;
                            ptr_ch->grid_bder_size += 1; // Increasing the number of border grid points in the array
                        }
                        //** Interior grid point **/
                        else
                        {

                            //** >> Adding the grid point to the interior array **/
                            ptr_ch->ptr_intr_grid_cell_idx_x[ptr_ch->grid_intr_size] = ii + ptr_ch->box_ts_x;
                            ptr_ch->ptr_intr_grid_cell_idx_y[ptr_ch->grid_intr_size] = jj + ptr_ch->box_ts_y;
                            ptr_ch->ptr_intr_grid_cell_idx_z[ptr_ch->grid_intr_size] = kk + ptr_ch->box_ts_z;
                            ptr_ch->ptr_intr_grid_idx[ptr_ch->grid_intr_size] = box_grid_idx_ch;
                            ptr_ch->grid_intr_size += 1; // Increasing the number of interior grid points in the array
                        }
                    }
                }
            }
        }

        //* Potential, Acceleration and density of the grid **/
        cap = (ptr_ch->box_real_dim_x + 1) * (ptr_ch->box_real_dim_y + 1) * (ptr_ch->box_real_dim_z + 1);
        ptr_ch->grid_properties_cap = cap;
        ptr_ch->ptr_pot = (vtype *)calloc(cap, sizeof(vtype));
        ptr_ch->ptr_pot_old = (vtype *)calloc(cap, sizeof(vtype));
        ptr_ch->ptr_ax = (vtype *)calloc(cap, sizeof(vtype));
        ptr_ch->ptr_ay = (vtype *)calloc(cap, sizeof(vtype));
        ptr_ch->ptr_az = (vtype *)calloc(cap, sizeof(vtype));
        ptr_ch->ptr_d = (vtype *)calloc(cap, sizeof(vtype));

        //** >> Tree structure **/
        ptr_node->pptr_chn[zone_idx] = ptr_ch; // Parent node pointing to child node i
        ptr_ch->ptr_pt = ptr_node;             // Child node i pointig to parent node

        // //** >> Boundary of the simulation box **/
        // if (ptr_ch->box_min_x == 0 || ptr_ch->box_max_x == (1 << ptr_ch->lv) - 1 ||
        //     ptr_ch->box_min_y == 0 || ptr_ch->box_max_y == (1 << ptr_ch->lv) - 1 ||
        //     ptr_ch->box_min_z == 0 || ptr_ch->box_max_z == (1 << ptr_ch->lv) - 1)
        // {
        //     ptr_ch->pbc_anomalies_due_to_the_boundary = true;
        // }
        // else
        // {
        //     ptr_ch->pbc_anomalies_due_to_the_boundary = false;
        // }

    } // End filling child nodes

    return _SUCCESS_;
} // end function fill_child_nodes

static int fill_child_nodes(struct node *ptr_node)
{

    struct node *ptr_ch; // Pointer to the child node

    int cell_idx;  // Cell index
    int aux_idx_x; // Auxiliary indexes
    int aux_idx_y;
    int aux_idx_z;

    int box_idx_x_ch;  // Box index in X direcction of child node
    int box_idx_y_ch;  // Box index in Y direcction of child node
    int box_idx_z_ch;  // Box index in Z direcction of child node
    int box_idx_ch;    // Box index
    int box_idxNbr_ch; // Box index of the neighboring child cell

    int box_idx_node; // Box index

    int box_idxNbr_i0_j0_k0_ch;    // Box index at x=0,y=0,z=0
    int box_idxNbr_im1_j0_k0_ch;   // Box index of the neighbor at x=-1,y=0,z=0
    int box_idxNbr_i0_jm1_k0_ch;   // Box index of the neighbor at x=0,y=-1,z=0
    int box_idxNbr_im1_jm1_k0_ch;  // Box index of the neighbor at x=-1,y=-1,z=0
    int box_idxNbr_i0_j0_km1_ch;   // Box index of the neighbor at x=0,y=0,z=-1
    int box_idxNbr_im1_j0_km1_ch;  // Box index of the neighbor at x=-1,y=0,z=-1
    int box_idxNbr_i0_jm1_km1_ch;  // Box index of the neighbor at x=0,y=-1,z=-1
    int box_idxNbr_im1_jm1_km1_ch; // Box index of the neighbor at x=-1,y=-1,z=-1

    int box_idx_ptcl_ch;

    int box_grid_idx_ch; // Box grid index

    bool is_bder_grid_point; // Ask if the grid point is interior
    bool grid_point_exist;

    int ptcl_idx;

    int pos_x; // Distance between the real box and the minimal box when the last is localized in the middle of the real one
    int pos_y;
    int pos_z;

    int cap;  // Capacity of arrays in child nodes
    int size; // Size or number of elements in some array in child nodes

    int lv = ptr_node->lv;

    int box_real_dim_X_ch;
    int box_real_dim_X_times_Y_ch;

    int grid_box_real_dim_X_ch;
    int grid_box_real_dim_X_times_Y_ch;

    int size_bder_grid_points_ch;
    int size_intr_grid_points_ch;

    //** >> Creating child nodes in parent node **/
    ptr_node->chn_cap = ptr_node->zones_cap; // Same amount than refinement zones
    ptr_node->chn_size = ptr_node->zones_size;

    if (ptr_node->chn_cap > 0)
    {
        ptr_node->pptr_chn = (struct node **)malloc(ptr_node->chn_cap * sizeof(struct node *)); // Allocating children pointers
    }

    for (int i = 0; i < ptr_node->chn_cap; i++)
    {
        ptr_node->pptr_chn[i] = NULL;
    }

    //** >> Filling child nodes **/
    for (int zone_idx = 0; zone_idx < ptr_node->zones_size; zone_idx++)
    {
        //** >> Create and Initialize child node **/
        ptr_ch = (struct node *)malloc(sizeof(struct node));
        initialize_node(ptr_ch);

        //** >> Global properties **/
        ptr_ch->ID = zone_idx;
        ptr_ch->lv = lv + 1;

        //** >> Boxes **/
        // MIN and MAX cell indexes values of the node.
        for (int j = 0; j < ptr_node->ptr_zone_size[zone_idx]; j++)
        {
            cell_idx = ptr_node->pptr_zones[zone_idx][j]; // Cell index in the parent zone
            if (ptr_ch->box_min_x > ptr_node->ptr_cell_idx_x[cell_idx])
            {
                ptr_ch->box_min_x = ptr_node->ptr_cell_idx_x[cell_idx];
            }
            if (ptr_ch->box_min_y > ptr_node->ptr_cell_idx_y[cell_idx])
            {
                ptr_ch->box_min_y = ptr_node->ptr_cell_idx_y[cell_idx];
            }
            if (ptr_ch->box_min_z > ptr_node->ptr_cell_idx_z[cell_idx])
            {
                ptr_ch->box_min_z = ptr_node->ptr_cell_idx_z[cell_idx];
            }
            if (ptr_ch->box_max_x < ptr_node->ptr_cell_idx_x[cell_idx])
            {
                ptr_ch->box_max_x = ptr_node->ptr_cell_idx_x[cell_idx];
            }
            if (ptr_ch->box_max_y < ptr_node->ptr_cell_idx_y[cell_idx])
            {
                ptr_ch->box_max_y = ptr_node->ptr_cell_idx_y[cell_idx];
            }
            if (ptr_ch->box_max_z < ptr_node->ptr_cell_idx_z[cell_idx])
            {
                ptr_ch->box_max_z = ptr_node->ptr_cell_idx_z[cell_idx];
            }
        }

        // Changing the min and max of the "minimal box" from parent units to child units
        ptr_ch->box_min_x = 2 * ptr_ch->box_min_x;
        ptr_ch->box_min_y = 2 * ptr_ch->box_min_y;
        ptr_ch->box_min_z = 2 * ptr_ch->box_min_z;
        ptr_ch->box_max_x = 2 * ptr_ch->box_max_x + 1;
        ptr_ch->box_max_y = 2 * ptr_ch->box_max_y + 1;
        ptr_ch->box_max_z = 2 * ptr_ch->box_max_z + 1;

        // Size of the "minimal box"
        ptr_ch->box_dim_x = ptr_ch->box_max_x - ptr_ch->box_min_x + 1;
        ptr_ch->box_dim_y = ptr_ch->box_max_y - ptr_ch->box_min_y + 1;
        ptr_ch->box_dim_z = ptr_ch->box_max_z - ptr_ch->box_min_z + 1;

        // Real dimensions of the boxcap = ptr_ch->box_cap;
        ptr_ch->box_real_dim_x = 5 > (n_exp - 1) ? (ptr_ch->box_dim_x + 10) : (ptr_ch->box_dim_x + 2 * n_exp - 2);
        ptr_ch->box_real_dim_y = 5 > (n_exp - 1) ? (ptr_ch->box_dim_y + 10) : (ptr_ch->box_dim_y + 2 * n_exp - 2);
        ptr_ch->box_real_dim_z = 5 > (n_exp - 1) ? (ptr_ch->box_dim_z + 10) : (ptr_ch->box_dim_z + 2 * n_exp - 2);

        ptr_ch->box_real_dim_x_old = ptr_ch->box_real_dim_x;
        ptr_ch->box_real_dim_y_old = ptr_ch->box_real_dim_y;
        // ptr_ch->box_real_dim_z_old = ptr_ch->box_real_dim_z;

        // Translations between cell array and box
        pos_x = (ptr_ch->box_real_dim_x - ptr_ch->box_dim_x) / 2; // Half of the distance of the box side less the "minimal box" side
        ptr_ch->box_ts_x = ptr_ch->box_min_x - pos_x;             // Every cell in the level l in the box must be subtracted this value to obtain the box index
        pos_y = (ptr_ch->box_real_dim_y - ptr_ch->box_dim_y) / 2;
        ptr_ch->box_ts_y = ptr_ch->box_min_y - pos_y;
        pos_z = (ptr_ch->box_real_dim_z - ptr_ch->box_dim_z) / 2;
        ptr_ch->box_ts_z = ptr_ch->box_min_z - pos_z;

        ptr_ch->box_ts_x_old = ptr_ch->box_ts_x;
        ptr_ch->box_ts_y_old = ptr_ch->box_ts_y;
        ptr_ch->box_ts_z_old = ptr_ch->box_ts_z;

        //** >> Cells in the node **/
        cap = 2 * 8 * ptr_node->ptr_zone_size[zone_idx]; // Capacity is the double of total child cells
        ptr_ch->cell_cap = cap;
        size = 8 * ptr_node->ptr_zone_size[zone_idx];
        ptr_ch->cell_size = size;
        ptr_ch->ptr_cell_idx_x = (int *)malloc(ptr_ch->cell_cap * sizeof(int));
        ptr_ch->ptr_cell_idx_y = (int *)malloc(ptr_ch->cell_cap * sizeof(int));
        ptr_ch->ptr_cell_idx_z = (int *)malloc(ptr_ch->cell_cap * sizeof(int));
        ptr_ch->ptr_box_idx = (int *)malloc(ptr_ch->cell_cap * sizeof(int));

        box_real_dim_X_ch = ptr_ch->box_real_dim_x;
        box_real_dim_X_times_Y_ch = ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;
        for (int j = 0; j < ptr_node->ptr_zone_size[zone_idx]; j++)
        {
            aux_idx_x = ptr_node->ptr_cell_idx_x[ptr_node->pptr_zones[zone_idx][j]] * 2;
            aux_idx_y = ptr_node->ptr_cell_idx_y[ptr_node->pptr_zones[zone_idx][j]] * 2;
            aux_idx_z = ptr_node->ptr_cell_idx_z[ptr_node->pptr_zones[zone_idx][j]] * 2;

            box_idx_ch = (aux_idx_x - ptr_ch->box_ts_x) + (aux_idx_y - ptr_ch->box_ts_y) * box_real_dim_X_ch + (aux_idx_z - ptr_ch->box_ts_z) * box_real_dim_X_times_Y_ch;
            for (int kk = 0; kk < 2; kk++)
            {
                for (int jj = 0; jj < 2; jj++)
                {
                    for (int ii = 0; ii < 2; ii++)
                    {
                        cell_idx = ii + jj * 2 + kk * 4 + j * 8; // Child cell index
                        ptr_ch->ptr_cell_idx_x[cell_idx] = aux_idx_x + ii;
                        ptr_ch->ptr_cell_idx_y[cell_idx] = aux_idx_y + jj;
                        ptr_ch->ptr_cell_idx_z[cell_idx] = aux_idx_z + kk;
                        ptr_ch->ptr_box_idx[cell_idx] = box_idx_ch + ii + jj * box_real_dim_X_ch + kk * box_real_dim_X_times_Y_ch;
                    }
                }
            }
        }

        // Filling the box status
        cap = ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y * ptr_ch->box_real_dim_z; // In general, the size of each side must be 3 times bigger than the same side of the "minimal box"
        ptr_ch->box_cap = cap;
        ptr_ch->ptr_box = (int *)malloc(cap * sizeof(int));
        ptr_ch->ptr_box_old = (int *)malloc(cap * sizeof(int));
        // Putting the value of NO-EXIST (-4) in every box index
        for (int j = 0; j < cap; j++)
        {
            ptr_ch->ptr_box[j] = -4;
        }
        // Changing from NO-EXIST (-4) to EXIST (-3) to all cells in the block of the child node
        for (int j = 0; j < ptr_ch->cell_size; j++)
        {
            box_idx_ch = ptr_ch->ptr_box_idx[j];
            ptr_ch->ptr_box[box_idx_ch] = -3;
        }

        //** >> Struct of cells (Particles and cell mass)
        //** >> Total mass in the node **/

        ptr_ch->ptr_cell_struct = (struct cell_struct *)malloc(ptr_ch->box_cap * sizeof(struct cell_struct));
        //ptr_ch->ptr_cell_struct_old = (struct cell_struct *)malloc(ptr_ch->box_cap * sizeof(struct cell_struct));
        for (int j = 0; j < ptr_ch->box_cap; j++)
        {
            initialize_cell_struct(&(ptr_ch->ptr_cell_struct[j]));
            //initialize_cell_struct(&(ptr_ch->ptr_cell_struct_old[j]));
        }

        for (int j = 0; j < ptr_node->ptr_zone_size[zone_idx]; j++)
        {
            box_idx_node = ptr_node->ptr_box_idx[ptr_node->pptr_zones[zone_idx][j]];
            // box_idx_x_node = ptr_node->ptr_cell_idx_x[ptr_node->pptr_zones[i][j]] - ptr_node->box_ts_x;
            // box_idx_y_node = ptr_node->ptr_cell_idx_y[ptr_node->pptr_zones[i][j]] - ptr_node->box_ts_y;
            // box_idx_z_node = ptr_node->ptr_cell_idx_z[ptr_node->pptr_zones[i][j]] - ptr_node->box_ts_z;
            // box_idx_node = box_idx_x_node + box_idx_y_node * ptr_node->box_real_dim_x + box_idx_z_node * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;
            if (ptr_node->ptr_cell_struct[box_idx_node].ptcl_size > 0)
            {

                // Allocating space fot the child particles in the 8  cells
                box_idx_x_ch = ptr_node->ptr_cell_idx_x[ptr_node->pptr_zones[zone_idx][j]] * 2 - ptr_ch->box_ts_x;
                box_idx_y_ch = ptr_node->ptr_cell_idx_y[ptr_node->pptr_zones[zone_idx][j]] * 2 - ptr_ch->box_ts_y;
                box_idx_z_ch = ptr_node->ptr_cell_idx_z[ptr_node->pptr_zones[zone_idx][j]] * 2 - ptr_ch->box_ts_z;
                box_idx_ch = box_idx_x_ch + box_idx_y_ch * box_real_dim_X_ch + box_idx_z_ch * box_real_dim_X_times_Y_ch;

                for (int kk = 0; kk < 2; kk++)
                {
                    for (int jj = 0; jj < 2; jj++)
                    {
                        for (int ii = 0; ii < 2; ii++)
                        {
                            box_idxNbr_ch = box_idx_ch + ii + jj * box_real_dim_X_ch + kk * box_real_dim_X_times_Y_ch;
                            cap = ptr_node->ptr_cell_struct[box_idx_node].ptcl_size;
                            ptr_ch->ptr_cell_struct[box_idxNbr_ch].ptcl_cap = cap;
                            ptr_ch->ptr_cell_struct[box_idxNbr_ch].ptr_ptcl = (int *)malloc(cap * sizeof(int));
                        }
                    }
                }

                // Putting parent cell particles into the child cell
                for (int k = 0; k < ptr_node->ptr_cell_struct[box_idx_node].ptcl_size; k++)
                {
                    ptcl_idx = ptr_node->ptr_cell_struct[box_idx_node].ptr_ptcl[k];

                    box_idx_ptcl_ch = ptcl_idx_to_box_idx(ptr_ch, ptcl_idx);

                    ptr_ch->ptr_cell_struct[box_idx_ptcl_ch].ptr_ptcl[ptr_ch->ptr_cell_struct[box_idx_ptcl_ch].ptcl_size] = ptcl_idx;
                    ptr_ch->ptr_cell_struct[box_idx_ptcl_ch].ptcl_size += 1;
                    ptr_ch->ptr_cell_struct[box_idx_ptcl_ch].cell_mass += GL_ptcl_mass[ptcl_idx];
                }

                ptr_ch->local_mass += ptr_node->ptr_cell_struct[box_idx_node].cell_mass; // Local mass
            }
        }

        // size = ptr_ch->cell_size / 8 * 18 + 9; // 27 * N - 9 * (N-1)

        // 27 * N - 9 * (N-1) = 18N+9 >= Total of grid points; N = number of cubics with 8 cells and 27 grid points
        // Maximum of border grid points = 9*(N*2+1) - (N*2-1) = 16N+10 (Imagine a row of cubics)
        // Maximum of interior grid points = (2k-1)^3 < 8k^3, where k^3 = N
        size_bder_grid_points_ch = ptr_ch->cell_size / 8 * 16 + 10;
        size_intr_grid_points_ch = ptr_ch->cell_size;

        //** >> Space checking of border grid points array**/
        ptr_ch->grid_intr_cap = size_intr_grid_points_ch;
        ptr_ch->grid_bder_cap = size_bder_grid_points_ch;
        ptr_ch->ptr_intr_grid_cell_idx_x = (int *)malloc(ptr_ch->grid_intr_cap * sizeof(int));
        ptr_ch->ptr_intr_grid_cell_idx_y = (int *)malloc(ptr_ch->grid_intr_cap * sizeof(int));
        ptr_ch->ptr_intr_grid_cell_idx_z = (int *)malloc(ptr_ch->grid_intr_cap * sizeof(int));
        ptr_ch->ptr_intr_grid_idx = (int *)malloc(ptr_ch->grid_intr_cap * sizeof(int));
        ptr_ch->ptr_bder_grid_cell_idx_x = (int *)malloc(ptr_ch->grid_bder_cap * sizeof(int));
        ptr_ch->ptr_bder_grid_cell_idx_y = (int *)malloc(ptr_ch->grid_bder_cap * sizeof(int));
        ptr_ch->ptr_bder_grid_cell_idx_z = (int *)malloc(ptr_ch->grid_bder_cap * sizeof(int));
        ptr_ch->ptr_bder_grid_idx = (int *)malloc(ptr_ch->grid_bder_cap * sizeof(int));

        //** >> Grid points **/
        grid_box_real_dim_X_ch = ptr_ch->box_real_dim_x + 1;
        grid_box_real_dim_X_times_Y_ch = (ptr_ch->box_real_dim_x + 1) * (ptr_ch->box_real_dim_y + 1);
        for (int kk = ptr_ch->box_min_z - ptr_ch->box_ts_z; kk < ptr_ch->box_max_z - ptr_ch->box_ts_z + 2; kk++)
        {
            for (int jj = ptr_ch->box_min_y - ptr_ch->box_ts_y; jj < ptr_ch->box_max_y - ptr_ch->box_ts_y + 2; jj++)
            {
                for (int ii = ptr_ch->box_min_x - ptr_ch->box_ts_x; ii < ptr_ch->box_max_x - ptr_ch->box_ts_x + 2; ii++)
                {
                    box_idx_ch = ii + jj * ptr_ch->box_real_dim_x + kk * ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;

                    grid_point_exist = false;

                    box_idxNbr_i0_j0_k0_ch = box_idx_ch;
                    box_idxNbr_im1_j0_k0_ch = box_idx_ch - 1;
                    box_idxNbr_i0_jm1_k0_ch = box_idx_ch - box_real_dim_X_ch;
                    box_idxNbr_im1_jm1_k0_ch = box_idx_ch - 1 - box_real_dim_X_ch;
                    box_idxNbr_i0_j0_km1_ch = box_idx_ch - box_real_dim_X_times_Y_ch;
                    box_idxNbr_im1_j0_km1_ch = box_idx_ch - 1 - box_real_dim_X_times_Y_ch;
                    box_idxNbr_i0_jm1_km1_ch = box_idx_ch - box_real_dim_X_ch - box_real_dim_X_times_Y_ch;
                    box_idxNbr_im1_jm1_km1_ch = box_idx_ch - 1 - box_real_dim_X_ch - box_real_dim_X_times_Y_ch;

                    // First
                    if (ptr_ch->ptr_box[box_idxNbr_i0_j0_k0_ch] > -4)
                    {
                        grid_point_exist = true;
                        is_bder_grid_point = false;
                        //** Connection to the right  **/

                        //** Connection to the left  **/
                        if (ptr_ch->ptr_box[box_idxNbr_im1_j0_k0_ch] < -3 &&
                            ptr_ch->ptr_box[box_idxNbr_im1_jm1_k0_ch] < -3 &&
                            ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] < -3 &&
                            ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] < -3)
                        {
                            is_bder_grid_point = true;
                        }
                        //** Backward connection   **/

                        //** Forward connection   **/
                        else if (ptr_ch->ptr_box[box_idxNbr_i0_jm1_k0_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_jm1_k0_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] < -3)
                        {
                            is_bder_grid_point = true;
                        }
                        //** Upward connection **/

                        //** Down connection **/
                        else if (ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] < -3)
                        {
                            is_bder_grid_point = true;
                        }
                    }
                    // Second
                    else if (ptr_ch->ptr_box[box_idxNbr_im1_j0_k0_ch] > -4)
                    {
                        grid_point_exist = true;
                        is_bder_grid_point = false;
                        //** Connection to the right  **/
                        if (ptr_ch->ptr_box[box_idxNbr_i0_jm1_k0_ch] < -3 &&
                            ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] < -3 &&
                            ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] < -3)
                        {
                            is_bder_grid_point = true;
                        }
                        //** Connection to the left  **/

                        //** Backward connection   **/

                        //** Forward connection   **/
                        else if (ptr_ch->ptr_box[box_idxNbr_i0_jm1_k0_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_jm1_k0_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] < -3)
                        {
                            is_bder_grid_point = true;
                        }
                        //** Upward connection **/

                        //** Down connection **/
                        else if (ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] < -3)
                        {
                            is_bder_grid_point = true;
                        }
                    }
                    // Third
                    else if (ptr_ch->ptr_box[box_idxNbr_i0_jm1_k0_ch] > -4)
                    {
                        grid_point_exist = true;
                        is_bder_grid_point = false;
                        //** Connection to the right  **/

                        //** Connection to the left  **/
                        if (ptr_ch->ptr_box[box_idxNbr_im1_jm1_k0_ch] < -3 &&
                            ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] < -3 &&
                            ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] < -3)
                        {
                            is_bder_grid_point = true;
                        }
                        //** Backward connection   **/
                        else if (ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] < -3)
                        {
                            is_bder_grid_point = true;
                        }
                        //** Forward connection   **/

                        //** Upward connection **/

                        //** Down connection **/
                        else if (ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] < -3)
                        {
                            is_bder_grid_point = true;
                        }
                    }
                    // Fourth
                    else if (ptr_ch->ptr_box[box_idxNbr_im1_jm1_k0_ch] > -4)
                    {
                        grid_point_exist = true;
                        is_bder_grid_point = false;
                        //** Connection to the right  **/
                        if (ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] < -3 &&
                            ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] < -3)
                        {
                            is_bder_grid_point = true;
                        }
                        //** Connection to the left  **/

                        //** Backward connection   **/
                        else if (ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] < -3)
                        {
                            is_bder_grid_point = true;
                        }
                        //** Forward connection   **/

                        //** Upward connection **/

                        //** Down connection **/
                        else if (ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] < -3)
                        {
                            is_bder_grid_point = true;
                        }
                    }
                    // Others
                    else if (ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] > -4 ||
                             ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] > -4 ||
                             ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] > -4 ||
                             ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] > -4)
                    {
                        grid_point_exist = true;
                        is_bder_grid_point = true;
                    }

                    if (grid_point_exist == true)
                    {
                        box_grid_idx_ch = ii + jj * grid_box_real_dim_X_ch + kk * grid_box_real_dim_X_times_Y_ch;
                        //** >> Adding the grid point **/

                        //** >> Border grid point**/
                        if (is_bder_grid_point == true)
                        {
                            //** >> Adding the grid point to the border array **/
                            ptr_ch->ptr_bder_grid_cell_idx_x[ptr_ch->grid_bder_size] = ii + ptr_ch->box_ts_x;
                            ptr_ch->ptr_bder_grid_cell_idx_y[ptr_ch->grid_bder_size] = jj + ptr_ch->box_ts_y;
                            ptr_ch->ptr_bder_grid_cell_idx_z[ptr_ch->grid_bder_size] = kk + ptr_ch->box_ts_z;
                            ptr_ch->ptr_bder_grid_idx[ptr_ch->grid_bder_size] = box_grid_idx_ch;
                            ptr_ch->grid_bder_size += 1; // Increasing the number of border grid points in the array
                        }
                        //** Interior grid point **/
                        else
                        {

                            //** >> Adding the grid point to the interior array **/
                            ptr_ch->ptr_intr_grid_cell_idx_x[ptr_ch->grid_intr_size] = ii + ptr_ch->box_ts_x;
                            ptr_ch->ptr_intr_grid_cell_idx_y[ptr_ch->grid_intr_size] = jj + ptr_ch->box_ts_y;
                            ptr_ch->ptr_intr_grid_cell_idx_z[ptr_ch->grid_intr_size] = kk + ptr_ch->box_ts_z;
                            ptr_ch->ptr_intr_grid_idx[ptr_ch->grid_intr_size] = box_grid_idx_ch;
                            ptr_ch->grid_intr_size += 1; // Increasing the number of interior grid points in the array
                        }
                    }
                }
            }
        }

        //* Potential, Acceleration and density of the grid **/
        cap = (ptr_ch->box_real_dim_x + 1) * (ptr_ch->box_real_dim_y + 1) * (ptr_ch->box_real_dim_z + 1);
        ptr_ch->grid_properties_cap = cap;
        ptr_ch->ptr_pot = (vtype *)calloc(cap, sizeof(vtype));
        ptr_ch->ptr_pot_old = (vtype *)calloc(cap, sizeof(vtype));
        ptr_ch->ptr_ax = (vtype *)calloc(cap, sizeof(vtype));
        ptr_ch->ptr_ay = (vtype *)calloc(cap, sizeof(vtype));
        ptr_ch->ptr_az = (vtype *)calloc(cap, sizeof(vtype));
        ptr_ch->ptr_d = (vtype *)calloc(cap, sizeof(vtype));

        //** >> Tree structure **/
        ptr_node->pptr_chn[zone_idx] = ptr_ch; // Parent node pointing to child node i
        ptr_ch->ptr_pt = ptr_node;             // Child node i pointig to parent node

        //** >> Boundary of the simulation box **/
        if (ptr_ch->box_min_x == 0 || ptr_ch->box_max_x == (1 << ptr_ch->lv) - 1 ||
            ptr_ch->box_min_y == 0 || ptr_ch->box_max_y == (1 << ptr_ch->lv) - 1 ||
            ptr_ch->box_min_z == 0 || ptr_ch->box_max_z == (1 << ptr_ch->lv) - 1)
        {
            ptr_ch->boundary_simulation_contact = true;
        }
        else
        {
            ptr_ch->boundary_simulation_contact = false;
        }

    } // End filling child nodes

    return _SUCCESS_;
} // end function fill_child_nodes

static int fill_tentacles(const struct node *ptr_node)
{
    int lv = ptr_node->lv - lmin + 1; // Children level

    int size = GL_tentacles_size[lv] + ptr_node->zones_size;

    if (space_check(&(GL_tentacles_cap[lv]), size, 4.0f, "p1n2", &(GL_tentacles[lv])) == _FAILURE_)
    {
        printf("Error, in space_check function\n");
        return _FAILURE_;
    }

    for (int i = 0; i < ptr_node->zones_size; i++)
    {
        //** >> Putting elements in the new tentacles **/
        GL_tentacles[lv][GL_tentacles_size[lv] + i] = ptr_node->pptr_chn[i];
    }

    //** Increasing the number of structs in the level lv **/
    GL_tentacles_size[lv] = size;

    //** >> Increasing the deepth of levels of the tentacles **/
    if (GL_tentacles_level_max < lv)
    {
        GL_tentacles_level_max++;
    }

    return _SUCCESS_;
} // end function fill_tentacles

int tree_construction(void)
{
    struct node *ptr_node; // Parent node

    int no_pts; // Number of parents in the cycle

    int lv = 0;

    while (lv < lmax - lmin && lv <= GL_tentacles_level_max)
    {

        //periodic boundary conditions
        if (boundary_type == 0)
        {
            while (GL_max_lv_ref_crosses_the_entire_simulation <= lmin + lv)
            {
                no_pts = GL_tentacles_size[lv];
                for (int i = 0; i < no_pts; i++)
                {
                    ptr_node = GL_tentacles[lv][i];
                    //** >> Filling the refinement cells array **/
                    if (fill_cell_ref_PERIODIC_BOUNDARY(ptr_node) == _FAILURE_)
                    {
                        printf("Error at function fill_cell_ref()\n");
                        return _FAILURE_;
                    }

                    //** >> Filling the different zones of refinement **/
                    if (fill_zones_ref_PERIODIC_BOUNDARY(ptr_node) == _FAILURE_)
                    {
                        printf("Error at function fill_zones_ref()\n");
                        return _FAILURE_;
                    }

                    //** >> Creating child nodes **/
                    if (ptr_node->zones_size > 0)
                    {
                        if (fill_child_nodes_PERIODIC_BOUNDARY(ptr_node) == _FAILURE_)
                        {
                            printf("Error at function create_child_nodes()\n");
                            return _FAILURE_;
                        }

                        //** >> Filling Tentacles for the next cycle at level of refinement
                        if (fill_tentacles(ptr_node) == _FAILURE_)
                        {
                            printf("Error at function fill_tentacles()\n");
                            return _FAILURE_;
                        }
                    }
                }

                lv++;
            }
        }
        else
        {
            no_pts = GL_tentacles_size[lv];
            //** >> For cycle over parent nodes **/
            for (int i = 0; i < no_pts; i++)
            {
                ptr_node = GL_tentacles[lv][i];

                //** >> Filling the refinement cells array **/
                if (fill_cell_ref(ptr_node) == _FAILURE_)
                {
                    printf("Error at function fill_cell_ref()\n");
                    return _FAILURE_;
                }

                //** >> Filling the different zones of refinement **/
                if (fill_zones_ref(ptr_node) == _FAILURE_)
                {
                    printf("Error at function fill_zones_ref()\n");
                    return _FAILURE_;
                }

                //** >> Creating child nodes **/
                if (ptr_node->zones_size > 0)
                {
                    if (fill_child_nodes(ptr_node) == _FAILURE_)
                    {
                        printf("Error at function create_child_nodes()\n");
                        return _FAILURE_;
                    }

                    //** >> Filling Tentacles for the next cycle at level of refinement
                    if (fill_tentacles(ptr_node) == _FAILURE_)
                    {
                        printf("Error at function fill_tentacles()\n");
                        return _FAILURE_;
                    }
                }

            } // End of cycle over parent nodes
            lv++;
        }
        

    } // End of cycle over refinement levels

    return _SUCCESS_;
}