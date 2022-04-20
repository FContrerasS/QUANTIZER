/*
 * tree_adaptation.c
 *
 * Adapting the tree to the next time-step iteration
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

#include "tree_adaptation.h"

static void updating_box_mass(struct node *ptr_node)
{

    int no_chn; // Number of child nodes
    struct node *ptr_ch = NULL; //child node
    vtype aux_box_mass; //Total mass in the 8 child cells

    int box_idx_x_ch; // Box index in X direcction of the child cell 
    int box_idx_y_ch; // Box index in Y direcction of the child cell
    int box_idx_z_ch; // Box index in Z direcction of the child cell
    int box_idx_ch;  // Box index of the child cell

    int box_idx_x_node;   // Box index in X direcction of the node cell
    int box_idx_y_node;   // Box index in Y direcction of the node cell
    int box_idx_z_node;   // Box index in Z direcction of the node cell
    int box_idx_node;  // Box index of the node cell

    no_chn = ptr_node->chn_size;

    for (int i = 0; i < no_chn; i++)    //Cycle over children
    {
        ptr_ch = ptr_node->pptr_chn[i];
        for (int j = 0; j < ptr_ch->cell_size; j += 8) // Cycle over packeges of 8 cells
        {
            //** >> Computing the mass of the 8 child cells **/
            box_idx_x_ch = ptr_ch->ptr_cell_idx_x[j] - ptr_ch->box_ts_x;
            box_idx_y_ch = ptr_ch->ptr_cell_idx_y[j] - ptr_ch->box_ts_y;
            box_idx_z_ch = ptr_ch->ptr_cell_idx_z[j] - ptr_ch->box_ts_z;

            aux_box_mass = 0;

            for (int kk = 0; kk < 2; kk++)
            {
                for (int jj = 0; jj < 2; jj++)
                {
                    for (int ii = 0; ii < 2; ii++)
                    {
                        box_idx_ch = (box_idx_x_ch + ii) + (box_idx_y_ch + jj) * ptr_ch->box_real_dim_x + (box_idx_z_ch + kk) * ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;
                        aux_box_mass += ptr_ch->ptr_box_mass[box_idx_ch];
                    }
                }
            }
            //** Re-defining the mass of the node cell corresponding to those 8 child cells
            box_idx_x_node = (ptr_ch->ptr_cell_idx_x[j] >> 1) - ptr_node->box_ts_x;
            box_idx_y_node = (ptr_ch->ptr_cell_idx_y[j] >> 1) - ptr_node->box_ts_y;
            box_idx_z_node = (ptr_ch->ptr_cell_idx_z[j] >> 1) - ptr_node->box_ts_z;
            box_idx_node = box_idx_x_node + box_idx_y_node * ptr_node->box_real_dim_x + box_idx_z_node * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;
            ptr_node->ptr_box_mass[box_idx_node] = aux_box_mass;
            total_mass += aux_box_mass;
        }
    }
    ptr_ch = NULL;
}

static void initialization_box_aux(struct node *ptr_node)
{

    int cap; // Maximum amount of cells in the box
    int size; // Number of cells in the node

    int box_idx_x; // Box index in X direcction
    int box_idx_y; // Box index in Y direcction
    int box_idx_z; // Box index in Z direcction
    int box_idx;   // Box index
    int box_idxNbr; // Box index in the neigborhood

    // Filling the box status
    cap = ptr_node->box_real_dim_x * ptr_node->box_real_dim_y * ptr_node->box_real_dim_z; // In general, the size of each side must be 3 times bigger than the same side of the "minimal box"
    size = ptr_node->cell_size;                                                                 // Number of cells in the block
    // Putting the value of NO-EXIST (-4) in every box index
    for (int i = 0; i < cap; i++)
    {
        ptr_node->ptr_box_aux[i] = -4;
    }
    // Changing from NO-EXIST (-4) to EXIST (-3) to all cells in the block of the child node
    for (int i = 0; i < size; i++)
    {
        box_idx_x = ptr_node->ptr_cell_idx_x[i] - ptr_node->box_ts_x;
        box_idx_y = ptr_node->ptr_cell_idx_y[i] - ptr_node->box_ts_y;
        box_idx_z = ptr_node->ptr_cell_idx_z[i] - ptr_node->box_ts_z;
        box_idx = box_idx_x + box_idx_y * ptr_node->box_real_dim_x + box_idx_z * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;
        ptr_node->ptr_box_aux[box_idx] = -3;
    }
    // Changing the EXIST (-3) TO BORDER (-2) to all border cells in the block
    for (int i = 0; i < size; i++)
    {
        box_idx_x = ptr_node->ptr_cell_idx_x[i] - ptr_node->box_ts_x;
        box_idx_y = ptr_node->ptr_cell_idx_y[i] - ptr_node->box_ts_y;
        box_idx_z = ptr_node->ptr_cell_idx_z[i] - ptr_node->box_ts_z;
        box_idx = box_idx_x + box_idx_y * ptr_node->box_real_dim_x + box_idx_z * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;

        for (int kk = -1; kk < 2; kk++)
        {
            for (int jj = -1; jj < 2; jj++)
            {
                for (int ii = -1; ii < 2; ii++)
                {
                    box_idxNbr = box_idx + ii + jj * ptr_node->box_real_dim_x + kk * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;
                    if (ptr_node->ptr_box_aux[box_idxNbr] == -4) // Border cells are those such that at least on of their first neighbors are NO-EXIST cells.
                    {
                        ptr_node->ptr_box_aux[box_idx] = -2;
                        ii = 2;
                        jj = 2;
                        kk = 2;
                    }
                }
            }
        }
    }
}

static void initialization_ref_aux(struct node *ptr_node)
{
    //ptr_node->cell_ref_size = 0;    //Number of cell to refine is equal to 0
    for (int i = 0; i < ptr_node->zones_size; i++)
    {
        ptr_node->ptr_zone_size[i] = 0;    //The zone i contains 0 elements
    }

    //ptr_node->zones_size = 0;   //The total number of zones of refinement is 0
}

static int fill_cell_ref(struct node *ptr_node)
{
    //** >> Adding cells which satisfy the refinement criteria or contains grandchild cells to the array ptr_cell_ref and chaning the box_aux the status of refinement -1 **/

    struct node *ptr_ch = NULL; // Child node

    int size; // Size of the refinement cells array
    
    int box_idx_x; //Box index in X direcction of the node
    int box_idx_y; // Box index in X direcction of the node
    int box_idx_z; // Box index in X direcction of the node
    int box_idx; // Box index of the node
    int box_idxNbr;   // Box index in the neigborhood

    int cell_idx_ch; //Cell index of the child node

    int cntr; // Counter used to add cell index to the ptr_cell_ref


    size = 0; // Initial size of the cell refined array must be 0
    //** >> Adding to refine (-1) to cells in the node wich contains grandchildren refinement cells
    for (int ch = 0; ch < ptr_node->chn_size; ch++) // Cycle over children
    {
        ptr_ch = ptr_node->pptr_chn[ch];
        for (int i = 0; i < ptr_ch->cell_ref_size; i++) // Cycle over refined cells of the child ch
        {
            cell_idx_ch = ptr_ch->ptr_cell_ref[i];
            box_idx_x = (ptr_ch->ptr_cell_idx_x[cell_idx_ch] >> 1) - ptr_node->box_ts_x;
            box_idx_y = (ptr_ch->ptr_cell_idx_y[cell_idx_ch] >> 1) - ptr_node->box_ts_y;
            box_idx_z = (ptr_ch->ptr_cell_idx_z[cell_idx_ch] >> 1) - ptr_node->box_ts_z;
            box_idx = box_idx_x + box_idx_y * ptr_node->box_real_dim_x + box_idx_z * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;

            if (ptr_node->ptr_box_aux[box_idx] == -3 || ptr_node->ptr_box_aux[box_idx] == -2) // Cell has not been added yet
            {
                size++; // +1 to the total number of cells to be refined
                //** >> Chaning the cell box status from EXIST (-3) to REFINEMENT REQUIRED (-1) **/
                ptr_node->ptr_box_aux[box_idx] = -1;
            }

            //** >> Changing the neighboring cell status **/
            for (int kk = -n_exp; kk < n_exp + 1; kk++)
            {
                for (int jj = -n_exp; jj < n_exp + 1; jj++)
                {
                    for (int ii = -n_exp; ii < n_exp + 1; ii++)
                    {
                        box_idxNbr = box_idx + ii + jj * ptr_node->box_real_dim_x + kk * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;
                        //** >> Asking if the neighboring cell has not been changed yet **/
                        // The border (-2) of the simulation can be added
                        if (ptr_node->ptr_box_aux[box_idxNbr] == -3 || ptr_node->ptr_box_aux[box_idxNbr] == -2) // Cell has not been added yet
                        {
                            size++; // +1 to the total number of cells to be refined
                            //** >> Chaning the cell box status from EXIST (-3) to REFINEMENT REQUIRED (-1) **/
                            ptr_node->ptr_box_aux[box_idxNbr] = -1;
                        }
                    }
                }
            }
        }
    }

    //** >> Changing the box status from EXIST (-3) to REFINEMENT REQUIRED (-1) of cells satisfying the refinement criterion **/
    for (int i = 0; i < ptr_node->cell_size; i++)
    {
        box_idx_x = ptr_node->ptr_cell_idx_x[i] - ptr_node->box_ts_x;
        box_idx_y = ptr_node->ptr_cell_idx_y[i] - ptr_node->box_ts_y;
        box_idx_z = ptr_node->ptr_cell_idx_z[i] - ptr_node->box_ts_z;
        box_idx = box_idx_x + box_idx_y * ptr_node->box_real_dim_x + box_idx_z * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;

        // Refinement criterion in the box_mass in no border box points
        if (ptr_node->ptr_box_aux[box_idx] != -2 && ptr_node->ptr_box_mass[box_idx] >= ref_criterion_mass) // No border (-2)
        {
            if (ptr_node->ptr_box_aux[box_idx] == -3) // Cell has not been added yet
            {
                size++;
                //** >> Chaning the cell box status from EXIST (-3) to REFINEMENT REQUIRED (-1) **/
                ptr_node->ptr_box_aux[box_idx] = -1;
            }

            //** >> Changing the neighboring cell status **/
            for (int kk = -n_exp; kk < n_exp + 1; kk++)
            {
                for (int jj = -n_exp; jj < n_exp + 1; jj++)
                {
                    for (int ii = -n_exp; ii < n_exp + 1; ii++)
                    {
                        box_idxNbr = box_idx + ii + jj * ptr_node->box_real_dim_x + kk * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;
                        //** >> Asking if the neighboring cell has not been changed yet **/
                        // The border (-2) of the simulation can be added
                        if (ptr_node->ptr_box_aux[box_idxNbr] == -3 || ptr_node->ptr_box_aux[box_idxNbr] == -2) // Cell has not been added yet
                        {
                            size++;
                            //** >> Chaning the cell box status from EXIST (-3) to REFINEMENT REQUIRED (-1) **/
                            ptr_node->ptr_box_aux[box_idxNbr] = -1;
                        }
                    }
                }
            }
        }
    }

    //** >> Space checking of the capacity of the refined cells **/
    if (space_check(&(ptr_node->cell_ref_cap), size, "p1i1", &(ptr_node->ptr_cell_ref)) == _FAILURE_)
    {
        printf("Error, in space_check function\n");
        return _FAILURE_;
    }

    //** >> Adding the infomation about size of the ptr_cell_ref array **/
    ptr_node->cell_ref_size = size;

    //** >> Adding cells refined to the array ptr_cell_ref **/
    cntr = 0; // Counter for the position in the cell refined array
    for (int i = 0; i < ptr_node->cell_size; i++)
    {
        box_idx_x = ptr_node->ptr_cell_idx_x[i] - ptr_node->box_ts_x;
        box_idx_y = ptr_node->ptr_cell_idx_y[i] - ptr_node->box_ts_y;
        box_idx_z = ptr_node->ptr_cell_idx_z[i] - ptr_node->box_ts_z;
        box_idx = box_idx_x + box_idx_y * ptr_node->box_real_dim_x + box_idx_z * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;

        if (ptr_node->ptr_box_aux[box_idx] == -1) // Cell require refinement
        {
            ptr_node->ptr_cell_ref[cntr] = i;
            cntr++;
        }
    }

    ptr_ch = NULL;

    return _SUCCESS_;
} // end function fill_cell_ref

static int fill_zones_ref(struct node *ptr_node)
{

    //** >> Filling the auxiliary diferent refinement zones **/
    int cntr;  // Counter Number of cells added to any refinement zone
    int cntr2; // Counter number of cells added per inspection
    int insp;  // Numer of inspected cells in the zone

    int zone_size;    // Number of cells in the zone
    int zone_idx_max; // Maximum id of the zone. It is equal to to the capacity in the number of zones
    int zone_idx;     // Index of the zone

    int box_idx_x; // Box index in X direcction
    int box_idx_y; // Box index in Y direcction
    int box_idx_z; // Box index in Z direcction
    int box_idx;   // Box index

    int box_idxNbr_x_plus;  // Box index in the neigborhood on the right
    int box_idxNbr_x_minus; // Box index in the neigborhood on the left
    int box_idxNbr_y_plus;  // Box index in the neigborhood behind
    int box_idxNbr_y_minus; // Box index in the neigborhood in front
    int box_idxNbr_z_plus;  // Box index in the neigborhood up
    int box_idxNbr_z_minus; // Box index in the neigborhood down

    int cell_idx;     // The cell index is simply i of the for loop
    int cell_ref_idx; // The index in the cell refined array ptr_cell_ref

    // The auxiliary array ptr_aux_idx will be used to store the box indexes of the zone

    //** >> Space checking of auxiliary array ptr_aux_idx **/
    // The size will always be bigger or equal than the number of refined cells

    if (space_check(&(ptr_node->aux_idx_cap), ptr_node->cell_ref_cap, "p1i1", &(ptr_node->ptr_aux_idx)) == _FAILURE_)
    {
        printf("Error, in space_check function\n");
        return _FAILURE_;
    }

    //** >>  Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0) **/
    cntr = 0;                              // Counter Number of cells added to any refinement zone
    cell_ref_idx = 0;                      // The index in the cell refined array ptr_cell_ref
    zone_idx = 0;                          // Initial zone corresponds to the cero zone
    zone_idx_max = ptr_node->zones_cap;    // Maximum id of the zone. It is equal to to the capacity in the number of zones
    while (ptr_node->cell_ref_size > cntr) // The loop while works as long as the number of cell addeed is less than the total refined cells
    {
        // Notes that the initiality we inspect the elements in the refined cell array until an element has been found that is not found in any of the current refinement zones
        cell_idx = ptr_node->ptr_cell_ref[cell_ref_idx]; // Index of the cells array in the node
        box_idx_x = ptr_node->ptr_cell_idx_x[cell_idx] - ptr_node->box_ts_x;
        box_idx_y = ptr_node->ptr_cell_idx_y[cell_idx] - ptr_node->box_ts_y;
        box_idx_z = ptr_node->ptr_cell_idx_z[cell_idx] - ptr_node->box_ts_z;
        box_idx = box_idx_x + box_idx_y * ptr_node->box_real_dim_x + box_idx_z * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;

        if (ptr_node->ptr_box_aux[box_idx] == -1) // A cell without zone has been founded
        {
            zone_size = 0; // Initial number of element in the zone

            //** >> Including the first element of the box to the auxiliary array ptr_aux_idx **/
            ptr_node->ptr_aux_idx[0] = box_idx;

            //** >>  Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0) **/
            ptr_node->ptr_box_aux[box_idx] = zone_idx;

            zone_size++; // +1 to the number of cells in the zone
            cntr++;      // +1 to the number of cells added in total

            //** >> Building of the zone of ID = zone_idx **/
            insp = 0;                                                  // Counter for the number of elements inspectioned in the current zone array
            while (zone_size > insp && ptr_node->cell_ref_size > cntr) // The cycle ends when all the elements of the zone have been inspected or when the number of cells refined is equal to the total number of cells added.
            {
                // Note that the number of elements in the zone increases if the neighbors of the inspected cell must be added to the zone

                cntr2 = 0; // Counter number of cells added per cell inspected
                box_idx = ptr_node->ptr_aux_idx[insp];

                box_idxNbr_x_plus = box_idx + 1;
                box_idxNbr_x_minus = box_idx - 1;
                box_idxNbr_y_plus = box_idx + ptr_node->box_real_dim_x;
                box_idxNbr_y_minus = box_idx - ptr_node->box_real_dim_x;
                box_idxNbr_z_plus = box_idx + ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;
                box_idxNbr_z_minus = box_idx - ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;

                //** Checking the nearest 6 neighbors of face
                // First neighbor
                if (ptr_node->ptr_box_aux[box_idxNbr_x_plus] == -1)
                {
                    ptr_node->ptr_aux_idx[zone_size + cntr2] = box_idxNbr_x_plus; // Including the neighboring element of the box to the auxiliary array
                    ptr_node->ptr_box_aux[box_idxNbr_x_plus] = zone_idx;              // Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0)
                    cntr2++;                                                      // +1 to the number of cells added in the current inspection
                }
                // Second neighbor
                if (ptr_node->ptr_box_aux[box_idxNbr_x_minus] == -1)
                {
                    ptr_node->ptr_aux_idx[zone_size + cntr2] = box_idxNbr_x_minus; // Including the neighboring element of the box to the auxiliary array
                    ptr_node->ptr_box_aux[box_idxNbr_x_minus] = zone_idx;              // Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0)
                    cntr2++;                                                       // +1 to the number of cells added in the current inspection
                }
                // Third neighbor
                if (ptr_node->ptr_box_aux[box_idxNbr_y_plus] == -1)
                {
                    ptr_node->ptr_aux_idx[zone_size + cntr2] = box_idxNbr_y_plus; // Including the neighboring element of the box to the auxiliary array
                    ptr_node->ptr_box_aux[box_idxNbr_y_plus] = zone_idx;              // Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0)
                    cntr2++;                                                      // +1 to the number of cells added in the current inspection
                }
                // Fourth neighbor
                if (ptr_node->ptr_box_aux[box_idxNbr_y_minus] == -1)
                {
                    ptr_node->ptr_aux_idx[zone_size + cntr2] = box_idxNbr_y_minus; // Including the neighboring element of the box to the auxiliary array
                    ptr_node->ptr_box_aux[box_idxNbr_y_minus] = zone_idx;              // Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0)
                    cntr2++;                                                       // +1 to the number of cells added in the current inspection
                }
                // Fifth neighbor
                if (ptr_node->ptr_box_aux[box_idxNbr_z_plus] == -1)
                {
                    ptr_node->ptr_aux_idx[zone_size + cntr2] = box_idxNbr_z_plus; // Including the neighboring element of the box to the auxiliary array
                    ptr_node->ptr_box_aux[box_idxNbr_z_plus] = zone_idx;              // Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0)
                    cntr2++;                                                      // +1 to the number of cells added in the current inspection
                }
                // Sixth neighbor
                if (ptr_node->ptr_box_aux[box_idxNbr_z_minus] == -1)
                {
                    ptr_node->ptr_aux_idx[zone_size + cntr2] = box_idxNbr_z_minus; // Including the neighboring element of the box to the auxiliary array
                    ptr_node->ptr_box_aux[box_idxNbr_z_minus] = zone_idx;              // Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0)
                    cntr2++;                                                       // +1 to the number of cells added in the current inspection
                }
                zone_size += cntr2; // Increasing the number of cells in the zone
                cntr += cntr2;      // Increasing the number of cells added to the zones
                insp++;             // Increasing the number of inspections
            }                       // End while cycle, now the box contains the the information about all cells of the zone "zone_idx"

            //** >> Space checking of refinement zones arrays, refinement capacity array and refinement size array **/
            if (zone_idx_max < zone_idx + 1)
            {
                if (space_check(&(zone_idx_max), zone_idx + 1, "p3i2i1i1", &(ptr_node->pptr_zones), &(ptr_node->ptr_zone_cap), &(ptr_node->ptr_zone_size)) == _FAILURE_)
                {
                    printf("Error, in space_check function\n");
                    return _FAILURE_;
                }

                // //** >> Initiazling the new allocated values **/
                // for (int i = 0; i < zone_idx_max - zone_idx; i++)
                // {
                //     ptr_node->pptr_zones[zone_idx_max - i - 1] = NULL;
                //     ptr_node->ptr_zone_cap[zone_idx_max - i - 1] = 0;
                // }
            }
            ptr_node->ptr_zone_size[zone_idx] = zone_size;
            zone_idx++; // Increasing the zone number
        }               // Zone defined in the box
        cell_ref_idx++;
    } // At this point the box_aux contains the information of all new refinement zones

    ptr_node->zones_cap = zone_idx_max; // Maximum amount of zones
    ptr_node->zones_size = zone_idx;    // Total amount of zones

    //** >> Space checking in each zone **/
    for (int i = 0; i < ptr_node->zones_size; i++)
    {
        if (space_check(&(ptr_node->ptr_zone_cap[i]), ptr_node->ptr_zone_size[i], "p1i1", &(ptr_node->pptr_zones[i])) == _FAILURE_)
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
        box_idx_x = ptr_node->ptr_cell_idx_x[cell_idx] - ptr_node->box_ts_x;
        box_idx_y = ptr_node->ptr_cell_idx_y[cell_idx] - ptr_node->box_ts_y;
        box_idx_z = ptr_node->ptr_cell_idx_z[cell_idx] - ptr_node->box_ts_z;
        box_idx = box_idx_x + box_idx_y * ptr_node->box_real_dim_x + box_idx_z * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;
        zone_idx = ptr_node->ptr_box_aux[box_idx];
        cntr = ptr_node->ptr_aux_idx[zone_idx];          // Counter the element in the zone "zone_idx"
        ptr_node->pptr_zones[zone_idx][cntr] = cell_idx; // Adding the index of the cell array in the block to the zone
        ptr_node->ptr_aux_idx[zone_idx] += 1;            // Counter the number of elements added in the zone "zone_idx"
    }
    return _SUCCESS_;
} // End function fill_zones_ref

static int create_links(struct node *ptr_node, int *links_old_ord_old, int *links_new_ord_old, int *links_old_ord_new, int *links_new_ord_new, int *ptr_links_cap)
{
    struct node *ptr_ch = NULL; // child node

    int cntr_links; // Counter the number of links between the new and old zones of refinement
    // int cntr_old_ref_cells; // Counter the total number of old cell refined in the node
    int cntr_cell_ch;       // Counter the cells of the child node
    int cntr_old_ref_zones; // total number of zones refined in the old box
    int cntr_links_plus;    // Counte links
    int cntr_link_elements; // Counter linked elements
    bool check_link;        // Check if the element is linked
    // int links_old_ord_old[ptr_node->zones_size];    // Storing the old zone of refinement id using old order
    // int links_new_ord_old[ptr_node->zones_size];    // Storing the new zone of refinement id using old order
    // int links_old_ord_new[ptr_node->zones_size];    // Storing the old zone of refinement id using new order
    // int links_new_ord_new[ptr_node->zones_size];    // Storing the new zone of refinement id using new order
    // int cell_idx_ch; // Index of the child cell
    int box_value_new; // Value in the new box

    bool create_link; // Decide if the link need to be created
    // int no_ptcl;      // Total number of particles in the node
    // int lv;

    int box_idx_x; // Box index at X direction
    int box_idx_y; // Box index at Y direction
    int box_idx_z; // Box index at Z direction
    int box_idx;   // Box index

    // int ptcl_idx;

    // int box_idx_ptcl;   // Box index

    int aux_min;
    int element_idx;
    int aux_int;

    // lv = ptr_node->lv;

    //** >> Creating the links between old and new refinement zones IDs **/

    cntr_links = 0;
    cntr_old_ref_zones = 0;

    // //** >> Space checking of the capacity of links order arrays **/
    if (space_check(ptr_links_cap, ptr_node->zones_size, "p4i1i1i1i1", &(links_old_ord_old), &(links_new_ord_old), &(links_old_ord_new), &(links_new_ord_new)) == _FAILURE_)
    {
        printf("Error, in space_check function\n");
        return _FAILURE_;
    }

    for (int i = 0; i < ptr_node->zones_size; i++)
    {
        links_new_ord_old[i] = -1;
        links_old_ord_old[i] = -1;
    }
    while (cntr_links < ptr_node->zones_size && cntr_old_ref_zones < ptr_node->chn_size)
    {
        ptr_ch = ptr_node->pptr_chn[cntr_old_ref_zones];
        cntr_cell_ch = 0;
        while (cntr_cell_ch < ptr_ch->cell_size)
        {
            box_idx_x = (ptr_ch->ptr_cell_idx_x[cntr_cell_ch] >> 1) - ptr_node->box_ts_x;
            box_idx_y = (ptr_ch->ptr_cell_idx_y[cntr_cell_ch] >> 1) - ptr_node->box_ts_y;
            box_idx_z = (ptr_ch->ptr_cell_idx_z[cntr_cell_ch] >> 1) - ptr_node->box_ts_z;
            box_idx = box_idx_x + box_idx_y * ptr_node->box_real_dim_x + box_idx_z * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;
            box_value_new = ptr_node->ptr_box_aux[box_idx];
            if (box_value_new < 0) // The old cell no longer requires refinement
            {
                cntr_cell_ch += 8;
            }
            else
            {
                create_link = true;
                for (int i = 0; i < cntr_links; i++)
                {
                    if (box_value_new == links_new_ord_old[i]) // The link already exist
                    {
                        create_link = false;
                        i = cntr_links;
                        cntr_cell_ch += 8;
                    }
                }
                if (create_link == true)
                {
                    links_new_ord_old[cntr_links] = box_value_new;
                    links_old_ord_old[cntr_links] = cntr_old_ref_zones;
                    cntr_links++;
                    cntr_cell_ch = ptr_ch->cell_size; // Breiking the while cycle
                }
            }
        } // End while cycle over the child node cells
        cntr_old_ref_zones++;
    } // End while cycle over all old refinement zones

    // printf("\nlv = %d, parent = %d\n", ptr_node->lv,ptr_node->ID);
    // printf("zones new = %d, zones old = %d\n", ptr_node->zones_size, ptr_node->chn_size);
    printf("cntr_links = %d\n", cntr_links);

    if (cntr_links < ptr_node->zones_size)
    {

        //** >> Old links **/
        if (cntr_links == ptr_node->chn_size)
        {
            for (int i = cntr_links; i < ptr_node->zones_size; i++)
            {
                // ptr_node->ptr_aux_idx[i - cntr_links] = i;
                links_old_ord_old[i] = i;
            }
            // cntr_links = ptr_node->zones_size;
        }
        else
        {
            cntr_link_elements = 0; // Counter linked elements
            cntr_links_plus = 0;    // Counte links old

            while (cntr_link_elements < cntr_links && cntr_links + cntr_links_plus < ptr_node->zones_size && cntr_links + cntr_links_plus < ptr_node->chn_size)
            {
                if (cntr_link_elements + cntr_links_plus == links_old_ord_old[cntr_link_elements])
                {
                    cntr_link_elements++;
                }
                else
                {
                    links_old_ord_old[cntr_links + cntr_links_plus] = cntr_link_elements + cntr_links_plus;
                    cntr_links_plus++;
                }
            }

            if (cntr_links + cntr_links_plus < ptr_node->zones_size)
            {
                for (int i = cntr_links + cntr_links_plus; i < ptr_node->zones_size; i++)
                {
                    // ptr_node->ptr_aux_idx[i - cntr_links] = i;
                    links_old_ord_old[i] = i;
                }
            }
        }

        //** >> New links **/
        cntr_link_elements = -1;
        for (int i = cntr_links; i < ptr_node->zones_size; i++)
        {
            check_link = true;
            while (check_link == true)
            {
                check_link = false;
                cntr_link_elements++;
                for (int j = 0; j < cntr_links; j++)
                {
                    if (links_new_ord_old[j] == cntr_link_elements)
                    {
                        check_link = true;
                        j = cntr_links;
                    }
                }
            }
            links_new_ord_old[i] = cntr_link_elements;
        }
    }

    printf("\n\n");
    printf("Links_old_ord_old = [ ");
    for (int i = 0; i < ptr_node->zones_size; i++)
    {
        printf("%d ", links_old_ord_old[i]);
    }
    printf("]\n");
    printf("Links_new_ord_old = [ ");
    for (int i = 0; i < ptr_node->zones_size; i++)
    {
        printf("%d ", links_new_ord_old[i]);
    }
    printf("]\n");

    //** >> Ordering links arrays **/
    //** >> Old order **/
    cntr_links_plus = 0;
    for (int i = 0; i < ptr_node->zones_size - 1; i++)
    {
        aux_min = ptr_node->chn_size + ptr_node->zones_size;
        element_idx = -1;
        cntr_link_elements = i;
        while (aux_min > cntr_links_plus && cntr_link_elements < ptr_node->zones_size)
        {
            if (aux_min > links_old_ord_old[cntr_link_elements])
            {
                aux_min = links_old_ord_old[cntr_link_elements];
                element_idx = cntr_link_elements;
            }
            cntr_link_elements++;
        }

        if (links_old_ord_old[i] != aux_min)
        {
            //** >> Old **/
            aux_int = links_old_ord_old[i];
            links_old_ord_old[i] = aux_min;
            links_old_ord_old[element_idx] = aux_int;
            //** >> New **/
            aux_int = links_new_ord_old[i];
            links_new_ord_old[i] = links_new_ord_old[element_idx];
            links_new_ord_old[element_idx] = aux_int;
        }
        cntr_links_plus = aux_min + 1;
    }

    //** >> New order **/
    for (int i = 0; i < ptr_node->zones_size; i++)
    {
        links_new_ord_new[i] = i;
    }
    for (int i = 0; i < ptr_node->zones_size; i++)
    {
        for (int j = 0; j < ptr_node->zones_size; j++)
        {
            if (links_new_ord_old[j] == i)
            {
                links_old_ord_new[i] = links_old_ord_old[j];
                j = ptr_node->zones_size;
            }
        }
    }
    // printf("\n\n");
    // printf("Links_old_ord_old = [   ");
    // for (int i = 0; i< ptr_node->zones_size; i++)
    // {
    //     printf("%d  ", links_old_ord_old[i]);
    // }
    // printf("]\n");
    // printf("Links_new_ord_old = [   ");
    // for (int i = 0; i< ptr_node->zones_size; i++)
    // {
    //     printf("%d  ", links_new_ord_old[i]);
    // }
    // printf("]\n");
    // printf("Links_old_ord_new = [   ");
    // for (int i = 0; i < ptr_node->zones_size; i++)
    // {
    //     printf("%d  ", links_old_ord_new[i]);
    // }
    // printf("]\n");
    // printf("Links_new_ord_new = [   ");
    // for (int i = 0; i < ptr_node->zones_size; i++)
    // {
    //     printf("%d  ", links_new_ord_new[i]);
    // }
    // printf("]\n");

    // //** >> Initializing the ptr_aux_idx array to be used as a counter of cells in each current child node**/
    // // Notes that the number of refined cells is always bigger or equal than the number of refined zones, so that the capacity of ptr_aux_idx is always enough to counter the number of refined zones
    // for (int i = 0; i < ptr_node->chn_size; i++)
    // {
    //     ptr_node->ptr_aux_idx[i] = 0;
    // }

    ptr_ch = NULL;

    return _SUCCESS_;
}

static int adapt_child_box_and_cells(struct node *ptr_node, const int *links_old_ord_old, const int *links_new_ord_old)
{
    struct node *ptr_ch = NULL;

    vtype *ptr_aux;

    int new_box_min_x;
    int new_box_min_y;
    int new_box_min_z;
    int new_box_max_x;
    int new_box_max_y;
    int new_box_max_z;

    int new_box_dim_x;
    int new_box_dim_y;
    int new_box_dim_z;

    int new_box_real_dim_x;
    int new_box_real_dim_y;
    int new_box_real_dim_z;

    int new_box_ts_x;
    int new_box_ts_y;
    int new_box_ts_z;

    int box_old_idx_x;
    int box_old_idx_y;
    int box_old_idx_z;
    int box_old_idx;
    int box_new_idx_x;
    int box_new_idx_y;
    int box_new_idx_z;
    int box_new_idx;

    int cell_idx; // The cell index is simply i of the for loop

    int new_zone_idx; // ID of the new zone of refinement

    int cap;
    int size;

    bool check_fit;

    int aux_int;

    int pos_x; // Distance between the real box and the minimal box when the last is localized in the middle of the real one
    int pos_y;
    int pos_z;

    int bder_box = 1 > n_exp ? 1 : n_exp;

    // Cycle over new refinement zones
    for (int zone_idx = 0; zone_idx < ptr_node->zones_size; zone_idx++)
    {
        if (links_old_ord_old[zone_idx] < ptr_node->chn_size)
        {
            new_zone_idx = links_new_ord_old[zone_idx];
            ptr_ch = ptr_node->pptr_chn[links_old_ord_old[zone_idx]];
            
            //** CELLS **/
            size = 8 * ptr_node->ptr_zone_size[new_zone_idx];
            //** >> Space checking of cells indexes of the child node
            if (space_check(&(ptr_ch->cell_cap), size, "p3i1i1i1", &(ptr_ch->ptr_cell_idx_x), &(ptr_ch->ptr_cell_idx_y), &(ptr_ch->ptr_cell_idx_z)) == _FAILURE_)
            {
                printf("Error, in space_check function\n");
                return _FAILURE_;
            }

            //** BOXES **/
            check_fit = true;

            new_box_min_x = INT_MAX;
            new_box_min_y = INT_MAX;
            new_box_min_z = INT_MAX;
            new_box_max_x = 0;
            new_box_max_y = 0;
            new_box_max_z = 0;
            //** >> The MIN and MAX of the set containig only the new cells to be refined
            for (int i = 0; i < ptr_node->ptr_zone_size[new_zone_idx]; i++)
            {
                cell_idx = ptr_node->pptr_zones[new_zone_idx][i];

                if (new_box_min_x > ptr_node->ptr_cell_idx_x[cell_idx])
                {
                    new_box_min_x = ptr_node->ptr_cell_idx_x[cell_idx];
                }
                if (new_box_min_y > ptr_node->ptr_cell_idx_y[cell_idx])
                {
                    new_box_min_y = ptr_node->ptr_cell_idx_y[cell_idx];
                }
                if (new_box_min_z > ptr_node->ptr_cell_idx_z[cell_idx])
                {
                    new_box_min_z = ptr_node->ptr_cell_idx_z[cell_idx];
                }
                if (new_box_max_x < ptr_node->ptr_cell_idx_x[cell_idx])
                {
                    new_box_max_x = ptr_node->ptr_cell_idx_x[cell_idx];
                }
                if (new_box_max_y < ptr_node->ptr_cell_idx_y[cell_idx])
                {
                    new_box_max_y = ptr_node->ptr_cell_idx_y[cell_idx];
                }
                if (new_box_max_z < ptr_node->ptr_cell_idx_z[cell_idx])
                {
                    new_box_max_z = ptr_node->ptr_cell_idx_z[cell_idx];
                }
            }

            // Changing the min and max of the "minimal box" from parent units to child units
            new_box_min_x = 2 * new_box_min_x;
            new_box_min_y = 2 * new_box_min_y;
            new_box_min_z = 2 * new_box_min_z;
            new_box_max_x = 2 * new_box_max_x + 1;
            new_box_max_y = 2 * new_box_max_y + 1;
            new_box_max_z = 2 * new_box_max_z + 1;

            //** >> Updating the The MIN and MAX of the set containig only the new cells to be refined
            ptr_ch->box_min_x = new_box_min_x;
            ptr_ch->box_min_y = new_box_min_y;
            ptr_ch->box_min_z = new_box_min_z;
            ptr_ch->box_max_x = new_box_max_x;
            ptr_ch->box_max_y = new_box_max_y;
            ptr_ch->box_max_z = new_box_max_z;

            //** >> Checking if the new dimensions fit in the old box **/
            aux_int = new_box_min_x - ptr_ch->box_ts_x;
            aux_int = aux_int < new_box_min_y - ptr_ch->box_ts_y ? aux_int : new_box_min_y - ptr_ch->box_ts_y;
            aux_int = aux_int < new_box_min_z - ptr_ch->box_ts_z ? aux_int : new_box_min_z - ptr_ch->box_ts_z;
            if (aux_int < bder_box)
            {
                check_fit = false;
            }
            aux_int = new_box_max_x - ptr_ch->box_ts_x;
            if (aux_int > ptr_ch->box_real_dim_x - bder_box)
            {
                check_fit = false;
            }
            aux_int = new_box_max_y - ptr_ch->box_ts_y;
            if (aux_int > ptr_ch->box_real_dim_y - bder_box)
            {
                check_fit = false;
            }
            aux_int = new_box_max_z - ptr_ch->box_ts_z;
            if (aux_int > ptr_ch->box_real_dim_z - bder_box)
            {
                check_fit = false;
            }
            //** >> The new box does not fit in the old box **/
            if (check_fit == false)
            {

                // Size of the "minimal new box"
                new_box_dim_x = new_box_max_x - new_box_min_x + 1;;
                new_box_dim_y = new_box_max_y - new_box_min_y + 1;;
                new_box_dim_z = new_box_max_z - new_box_min_z + 1;;

                // Real dimensions of the new box
                new_box_real_dim_x = 10 > 2 * n_exp - 2 ? ptr_ch->box_dim_x + 10 : ptr_ch->box_dim_x + 2 * n_exp - 2;
                new_box_real_dim_y = 10 > 2 * n_exp - 2 ? ptr_ch->box_dim_y + 10 : ptr_ch->box_dim_y + 2 * n_exp - 2;
                new_box_real_dim_z = 10 > 2 * n_exp - 2 ? ptr_ch->box_dim_z + 10 : ptr_ch->box_dim_z + 2 * n_exp - 2;

                // Translations between cell array and new box
                pos_x = (new_box_real_dim_x - new_box_dim_x) / 2; // Half of the distance of the box side less the "minimal box" side
                new_box_ts_x = new_box_min_x - pos_x;             // Every cell in the level l in the box must be subtracted this value to obtain the box index
                pos_y = (new_box_real_dim_y - new_box_dim_y) / 2;
                new_box_ts_y = new_box_min_y - pos_y;
                pos_z = (new_box_real_dim_z - new_box_dim_z) / 2;
                new_box_ts_z = new_box_min_z - pos_z;

                cap = ptr_ch->box_cap;
                size = new_box_real_dim_x * new_box_real_dim_y * new_box_real_dim_z;

                //** >> Space checking of boxes: box, box_aux and box_mass_aux
                if (space_check(&(cap), size, "p3i1i1v1", &(ptr_ch->ptr_box), &(ptr_ch->ptr_box_aux), &(ptr_ch->ptr_box_mass_aux)) == _FAILURE_)
                {
                    printf("Error, in space_check function\n");
                    return _FAILURE_;
                }

                // Initializing box at NO-EXIST (-4), and mass box in 0
                for (int j = 0; j < size; j++)
                {
                    ptr_ch->ptr_box[j] = -4;
                    ptr_ch->ptr_box_mass_aux[j] = 0;
                }

                //** Adapting box mass **/
                for (int j = 0; j < ptr_ch->cell_size; j++)
                {
                    //Old
                    box_old_idx_x = ptr_ch->ptr_cell_idx_x[j] - ptr_ch->box_ts_x;
                    box_old_idx_y = ptr_ch->ptr_cell_idx_y[j] - ptr_ch->box_ts_y;
                    box_old_idx_z = ptr_ch->ptr_cell_idx_z[j] - ptr_ch->box_ts_z;
                    box_old_idx = box_old_idx_x + box_old_idx_y * ptr_ch->box_real_dim_x + box_old_idx_z * ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;

                    //New
                    box_new_idx_x = ptr_ch->ptr_cell_idx_x[j] - new_box_ts_x;
                    box_new_idx_y = ptr_ch->ptr_cell_idx_y[j] - new_box_ts_y;
                    box_new_idx_z = ptr_ch->ptr_cell_idx_z[j] - new_box_ts_z;
                    box_new_idx = box_new_idx_x + box_new_idx_y * new_box_real_dim_x + box_new_idx_z * new_box_real_dim_x * new_box_real_dim_y;
                    
                    //Transfer the mass cell from old to new auxiliary mass box cell
                    ptr_ch->ptr_box_mass_aux[box_new_idx] = ptr_ch->ptr_box_mass[box_old_idx];
                }

                cap = ptr_ch->box_cap;
                //** >> Reallocating the new box if it is necessary
                if (space_check(&(cap), size, "p1v1", &(ptr_ch->ptr_box_mass)) == _FAILURE_)
                {
                    printf("Error, in space_check function\n");
                    return _FAILURE_;
                }

                // Transfer auxiliary mass box to box mass
                ptr_aux = ptr_ch->ptr_box_mass_aux;
                ptr_ch->ptr_box_mass_aux = ptr_ch->ptr_box_mass;
                ptr_ch->ptr_box_mass = ptr_aux;

                //** Updating new box values

                // Size of the "minimal box"
                ptr_ch->box_dim_x = new_box_dim_x;
                ptr_ch->box_dim_y = new_box_dim_y;
                ptr_ch->box_dim_z = new_box_dim_z;

                // Real dimensions of the new box
                ptr_ch->box_real_dim_x = new_box_real_dim_x;
                ptr_ch->box_real_dim_y = new_box_real_dim_y;
                ptr_ch->box_real_dim_z = new_box_real_dim_z;

                // Translations between cell array and box
                ptr_ch->box_ts_x = new_box_ts_x; // Every cell in the level l in the box must be subtracted this value to obtain the box index
                ptr_ch->box_ts_y = new_box_ts_y;
                ptr_ch->box_ts_z = new_box_ts_z;

                //New box capacity
                ptr_ch->box_cap = cap;
            }

            // //** >> Density initialization **/
            // size = (ptr_ch->box_real_dim_x + 1) * (ptr_ch->box_real_dim_y + 1) * (ptr_ch->box_real_dim_z + 1);
            // for (int j = 0; j < size; j++)
            // {
            //     ptr_ch->ptr_d[j] = 0;
            // }
        }
    }

    ptr_ch = NULL;

    return _SUCCESS_;
}

static int create_new_child_nodes(struct node *ptr_node, const int *links_old_ord_old, const int *links_new_ord_old)
{

    struct node *ptr_ch = NULL; // child node
    int cell_idx;        // The cell index is simply i of the for loop

    int pos_x;  // Distance between the real box and the minimal box when the last is localized in the middle of the real one
    int pos_y;
    int pos_z;

    int cap;  // Capacity of arrays in child nodes
    //int size; // Size or number of elements in some array in child nodes

    //int box_idx_x; // Box index at X direction
    //int box_idx_y; // Box index at Y direction
    //int box_idx_z; // Box index at Z direction
    //int box_idx;   // Box index
    //int box_idxNbr; // Box index in the neigborhood

    //** >> Space checking in the number of child nodes of ptr_node
    if (space_check(&(ptr_node->chn_cap), ptr_node->zones_size, "p1n2", &(ptr_node->pptr_chn)) == _FAILURE_)
    {
        printf("Error, in space_check function\n");
        return _FAILURE_;
    }

    for (int i = ptr_node->chn_size; i < ptr_node->zones_size; i++)
    {
        ptr_ch = (struct node *)malloc(sizeof(struct node));
        initialize_node(ptr_ch);

        //** >> Global node properties **/
        ptr_ch->ID = i;
        ptr_ch->lv = ptr_node->lv + 1;

        //** >> Cells in the node **/
        ptr_ch->cell_cap = 8 * ptr_node->ptr_zone_cap[links_new_ord_old[i]];
        ptr_ch->ptr_cell_idx_x = (int *)malloc(ptr_ch->cell_cap * sizeof(int));
        ptr_ch->ptr_cell_idx_y = (int *)malloc(ptr_ch->cell_cap * sizeof(int));
        ptr_ch->ptr_cell_idx_z = (int *)malloc(ptr_ch->cell_cap * sizeof(int));

        //** >> Particles in the node **/
        ptr_ch->ptcl_cap = 8;
        ptr_ch->ptr_ptcl = (int *)malloc(ptr_node->ptcl_cap * sizeof(int));

        //** >> Boxes **/
        // MIN and MAX cell indexes values of the node.
        for (int j = 0; j < ptr_node->ptr_zone_size[links_new_ord_old[i]]; j++)
        {
            cell_idx = ptr_node->pptr_zones[links_new_ord_old[i]][j]; // Cell index in the zone
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

        // Real dimensions of the box, or real capacity
        // ptr_ch->box_real_dim_x = fmax(ptr_ch->box_dim_x * 3, ptr_ch->box_dim_x + 2 * n_exp - 2);
        // ptr_ch->box_real_dim_y = fmax(ptr_ch->box_dim_y * 3, ptr_ch->box_dim_y + 2 * n_exp - 2);
        // ptr_ch->box_real_dim_z = fmax(ptr_ch->box_dim_z * 3, ptr_ch->box_dim_z + 2 * n_exp - 2);
        // ptr_ch->box_real_dim_x = fmax(ptr_ch->box_dim_x + 10, ptr_ch->box_dim_x + 2 * n_exp - 2);
        // ptr_ch->box_real_dim_y = fmax(ptr_ch->box_dim_y + 10, ptr_ch->box_dim_y + 2 * n_exp - 2);
        // ptr_ch->box_real_dim_z = fmax(ptr_ch->box_dim_z + 10, ptr_ch->box_dim_z + 2 * n_exp - 2);
        ptr_ch->box_real_dim_x = ptr_ch->box_dim_x + 10 > ptr_ch->box_dim_x + 2 * n_exp - 2 ? ptr_ch->box_dim_x + 10 : ptr_ch->box_dim_x + 2 * n_exp - 2;
        ptr_ch->box_real_dim_y = ptr_ch->box_dim_y + 10 > ptr_ch->box_dim_y + 2 * n_exp - 2 ? ptr_ch->box_dim_y + 10 : ptr_ch->box_dim_y + 2 * n_exp - 2;
        ptr_ch->box_real_dim_z = ptr_ch->box_dim_z + 10 > ptr_ch->box_dim_z + 2 * n_exp - 2 ? ptr_ch->box_dim_z + 10 : ptr_ch->box_dim_z + 2 * n_exp - 2;

        // Translations between cell array and box
        pos_x = (ptr_ch->box_real_dim_x - ptr_ch->box_dim_x) / 2; // Half of the distance of the box side less the "minimal box" side
        ptr_ch->box_ts_x = ptr_ch->box_min_x - pos_x;             // Every cell in the level l in the box must be subtracted this value to obtain the box index
        pos_y = (ptr_ch->box_real_dim_y - ptr_ch->box_dim_y) / 2;
        ptr_ch->box_ts_y = ptr_ch->box_min_y - pos_y;
        pos_z = (ptr_ch->box_real_dim_z - ptr_ch->box_dim_z) / 2;
        ptr_ch->box_ts_z = ptr_ch->box_min_z - pos_z;

        // Filling the box status
        cap = ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y * ptr_ch->box_real_dim_z; // In general, the size of each side must be 3 times bigger than the same side of the "minimal box"
        ptr_ch->box_cap = cap;
        //size = 8 * ptr_node->ptr_zone_size[links_new_ord_old[i]]; // Number of cells in the block
        ptr_ch->ptr_box = (int *)malloc(cap * sizeof(int));
        ptr_ch->ptr_box_aux = (int *)malloc(cap * sizeof(int));

        // Putting the value of NO-EXIST (-4) in every box index
        for (int j = 0; j < cap; j++)
        {
            ptr_ch->ptr_box[j] = -4;
        }
        // // Changing from NO-EXIST (-4) to EXIST (-3) to all cells in the block of the child node
        // for (int j = 0; j < size; j++)
        // {
        //     box_idx_x = ptr_ch->ptr_cell_idx_x[j] - ptr_ch->box_ts_x;
        //     box_idx_y = ptr_ch->ptr_cell_idx_y[j] - ptr_ch->box_ts_y;
        //     box_idx_z = ptr_ch->ptr_cell_idx_z[j] - ptr_ch->box_ts_z;
        //     box_idx = box_idx_x + box_idx_y * ptr_ch->box_real_dim_x + box_idx_z * ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;
        //     ptr_ch->ptr_box[box_idx] = -3;
        // }
        // // Changing the EXIST (-3) TO BORDER (-2) to all border cells in the block
        // for (int j = 0; j < size; j++)
        // {
        //     box_idx_x = ptr_ch->ptr_cell_idx_x[j] - ptr_ch->box_ts_x;
        //     box_idx_y = ptr_ch->ptr_cell_idx_y[j] - ptr_ch->box_ts_y;
        //     box_idx_z = ptr_ch->ptr_cell_idx_z[j] - ptr_ch->box_ts_z;
        //     box_idx = box_idx_x + box_idx_y * ptr_ch->box_real_dim_x + box_idx_z * ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;

        //     for (int kk = -1; kk < 2; kk++)
        //     {
        //         for (int jj = -1; jj < 2; jj++)
        //         {
        //             for (int ii = -1; ii < 2; ii++)
        //             {
        //                 box_idxNbr = box_idx + ii + jj * ptr_ch->box_real_dim_x + kk * ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;
        //                 if (ptr_ch->ptr_box[box_idxNbr] == -4) // Border cells are those such that at least on of their first neighbors are NO-EXIST cells.
        //                 {
        //                     ptr_ch->ptr_box[box_idx] = -2;
        //                     ii = 2;
        //                     jj = 2;
        //                     kk = 2;
        //                 }
        //             }
        //         }
        //     }
        // }

        //** >> Grid points **/
        //They are defined after the box definition

        //** >> Refinement Criterion **/
        cap = ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y * ptr_ch->box_real_dim_z;
        ptr_ch->ptr_box_mass = (vtype *)calloc(cap, sizeof(vtype));
        ptr_ch->ptr_box_mass_aux = (vtype *)calloc(cap, sizeof(vtype));

        //** >> Potential, acceleration and density of the grid **/
        cap = (ptr_ch->box_real_dim_x + 1) * (ptr_ch->box_real_dim_y + 1) * (ptr_ch->box_real_dim_z + 1);
        ptr_ch->grid_properties_cap = cap;
        ptr_ch->ptr_pot = (vtype *)calloc(cap, sizeof(vtype));
        ptr_ch->ptr_ax = (vtype *)calloc(cap, sizeof(vtype));
        ptr_ch->ptr_ay = (vtype *)calloc(cap, sizeof(vtype));
        ptr_ch->ptr_az = (vtype *)calloc(cap, sizeof(vtype));
        ptr_ch->ptr_d = (vtype *)calloc(cap, sizeof(vtype));

        //** >> Tree structure **/
        ptr_ch->ptr_pt = ptr_node;
        ptr_node->pptr_chn[i] = ptr_ch;
    }

    ptr_ch = NULL;

    return _SUCCESS_;
}

static void remov_cells_nolonger_require_refinement(struct node *ptr_node, const int *links_old_ord_old, const int *links_new_ord_old)
{

    struct node *ptr_ch = NULL;

    int box_idx_x_ch; // Box index in X direcction of the child cell
    int box_idx_y_ch; // Box index in Y direcction of the child cell
    int box_idx_z_ch; // Box index in Z direcction of the child cell
    int box_idx_ch;   // Box index of the child cell

    int box_idx_x_node; // Box index in X direcction of the node cell
    int box_idx_y_node; // Box index in Y direcction of the node cell
    int box_idx_z_node; // Box index in Z direcction of the node cell
    int box_idx_node;   // Box index of the node cell

    int box_idx_x_ptcl; // Box index in X direcction of the particles
    int box_idx_y_ptcl; // Box index in Y direcction of the particles
    int box_idx_z_ptcl; // Box index in Z direcction of the particles
    int box_idx_ptcl;   // Box index

    int ptcl_idx;

    int no_ptcl;      // Total number of particles in the node
    int no_cells;      // Total number of cells in the node

    

    

    int lv;

    lv = ptr_node->lv;

    // Cycle over new refinement zones
    for (int zone_idx = 0; zone_idx < ptr_node->zones_size; zone_idx++)
    {
        if (links_old_ord_old[zone_idx] < ptr_node->chn_size)
        {
            ptr_ch = ptr_node->pptr_chn[links_old_ord_old[zone_idx]];

            no_ptcl = ptr_ch->ptcl_size;
            no_cells = ptr_ch->cell_size;
            for (int cell_idx = 0; cell_idx < no_cells; cell_idx += 8)
            {
                box_idx_x_node = (ptr_ch->ptr_cell_idx_x[cell_idx] >> 1) - ptr_node->box_ts_x;
                box_idx_y_node = (ptr_ch->ptr_cell_idx_y[cell_idx] >> 1) - ptr_node->box_ts_y;
                box_idx_z_node = (ptr_ch->ptr_cell_idx_z[cell_idx] >> 1) - ptr_node->box_ts_z;
                box_idx_node = box_idx_x_node + box_idx_y_node * ptr_node->box_real_dim_x + box_idx_z_node * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;
                //** >> The child cell no longer requires refinement **/
                if (ptr_node->ptr_box_aux[box_idx_node] < 0)
                {
                    
                    //** >> Removing particles in the cell from the child node **/
                    //Updating local mass and particles array
                    for (int i = 0; i < no_ptcl; i++)
                    {
                        ptcl_idx = ptr_ch->ptr_ptcl[i];
                        box_idx_x_ptcl = GL_ptcl_x[ptcl_idx] * (1 << lv) - ptr_node->box_ts_x; // Particle indexes in the level
                        box_idx_y_ptcl = GL_ptcl_y[ptcl_idx] * (1 << lv) - ptr_node->box_ts_y;
                        box_idx_z_ptcl = GL_ptcl_z[ptcl_idx] * (1 << lv) - ptr_node->box_ts_z;
                        box_idx_ptcl = box_idx_x_ptcl + box_idx_y_ptcl * ptr_node->box_real_dim_x + box_idx_z_ptcl * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;
                        if (box_idx_ptcl == box_idx_node)
                        {
                            ptr_ch->ptr_ptcl[i] = ptr_ch->ptr_ptcl[no_ptcl - 1];
                            ptr_ch->local_mass -= GL_ptcl_mass[ptcl_idx];
                            no_ptcl--; // The total number of particle decrease
                            i--;       // The last element that was just moved to the current position should also must be analized
                        }
                    }

                    box_idx_x_ch = ptr_ch->ptr_cell_idx_x[cell_idx] - ptr_ch->box_ts_x;
                    box_idx_y_ch = ptr_ch->ptr_cell_idx_y[cell_idx] - ptr_ch->box_ts_y;
                    box_idx_z_ch = ptr_ch->ptr_cell_idx_z[cell_idx] - ptr_ch->box_ts_z;
                    
                    //** >> Updating the mass box **/
                    for (int kk = 0; kk < 2; kk++)
                    {
                        for (int jj = 0; jj < 2; jj++)
                        {
                            for(int ii = 0; ii < 2; ii++)
                            {
                                box_idx_ch = (box_idx_x_ch + ii) + (box_idx_y_ch + jj) * ptr_ch->box_real_dim_x + (box_idx_z_ch + kk) * ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;
                                ptr_ch->ptr_box_mass[box_idx_ch] = 0;
                            }
                        }
                    }

                    
                    //** >> Removing the cells **/
                    for (int j = 0; j < 8; j++)
                    {
                        ptr_ch->ptr_cell_idx_x[cell_idx + j] = ptr_ch->ptr_cell_idx_x[no_cells - 8 + j];
                        ptr_ch->ptr_cell_idx_y[cell_idx + j] = ptr_ch->ptr_cell_idx_y[no_cells - 8 + j];
                        ptr_ch->ptr_cell_idx_z[cell_idx + j] = ptr_ch->ptr_cell_idx_z[no_cells - 8 + j];
                    }
                    no_cells -= 8;
                    cell_idx -= 8;
                }
            }
            ptr_ch->ptcl_size = no_ptcl;
            ptr_ch->cell_size = no_cells;
        }
    }
    ptr_ch = NULL;
}

static int moving_old_child_to_new_child(struct node *ptr_node, const int *links_old_ord_old, const int *links_new_ord_old, const int *links_old_ord_new, const int *links_new_ord_new)
{
    struct node *ptr_ch_A = NULL;
    struct node *ptr_ch_B = NULL;
    int no_cells_ch_A; // Total number of cells in the child node A
    int no_cells_ch_B; // Total number of cells in the child node A

    int no_ptcl_ch_A; // Total number of particles in the child node A
    int no_ptcl_ch_B; // Total number of particles in the child node B

    int cntr_ch; //Counter the number of children

    int new_zone_idx; // ID of the new zone of refinement

    int box_idx_x_node; // Box index in X direcction of the node cell
    int box_idx_y_node; // Box index in Y direcction of the node cell
    int box_idx_z_node; // Box index in Z direcction of the node cell
    int box_idx_node;   // Box index of the node cell

    int box_idx_x_ch_A; // Box index in X direcction of the child cell
    int box_idx_y_ch_A; // Box index in Y direcction of the child cell
    int box_idx_z_ch_A; // Box index in Z direcction of the child cell
    int box_idx_x_ch_B; // Box index in X direcction of the child cell
    int box_idx_y_ch_B; // Box index in Y direcction of the child cell
    int box_idx_z_ch_B; // Box index in Z direcction of the child cell
    int box_idx_ch_A;   // Box index of the child cell
    int box_idx_ch_B;   // Box index of the child cell

    int box_idxNbr_ch; // Box index in the neigborhood

    int ptcl_idx;

    int box_idx_x_ptcl_ch_A; // Box index in X direcction of the particles
    int box_idx_y_ptcl_ch_A; // Box index in Y direcction of the particles
    int box_idx_z_ptcl_ch_A; // Box index in Z direcction of the particles
    int box_idx_ptcl_ch_A;   // Box index

    int box_idx_x_ptcl_ch_B; // Box index in X direcction of the particles
    int box_idx_y_ptcl_ch_B; // Box index in Y direcction of the particles
    int box_idx_z_ptcl_ch_B; // Box index in Z direcction of the particles
    int box_idx_ptcl_ch_B;   // Box index

    int box_idx_x_ptcl_node; // Box index in X direcction of the particles
    int box_idx_y_ptcl_node; // Box index in Y direcction of the particles
    int box_idx_z_ptcl_node; // Box index in Z direcction of the particles
    int box_idx_ptcl_node;   // Box index

    int lv;

    lv = ptr_node->lv;

    printf("A\n");


for(int i = 0; i< ptr_node->zones_size; i++)
{
    printf("%d, zone cell-size = %d\n",i, ptr_node->ptr_zone_size[i]);
}

printf("\n\n");

for (int i = 0; i < ptr_node->chn_size; i++)
{
    printf("%d, child cell-size = %d, cap = %d\n", i, ptr_node->pptr_chn[i]->cell_size, ptr_node->pptr_chn[i]->cell_cap);
}

    // Cycle over new refinement zones,
    // Case of old child nodes to be reused
    for (int zone_idx = 0; zone_idx < ptr_node->zones_size; zone_idx++)
    {
        if (links_old_ord_old[zone_idx] < ptr_node->chn_size)
        {
            new_zone_idx = links_new_ord_old[zone_idx];
            //printf("%snew zone = %d, old zone = %d, zone_idx = %d%s\n", KBLU, new_zone_idx, links_old_ord_old[zone_idx], zone_idx, KNRM);
            ptr_ch_A = ptr_node->pptr_chn[links_old_ord_old[zone_idx]];
            // printf("ID ch A = %d\n", ptr_ch_A->ID);
            no_ptcl_ch_A = ptr_ch_A->ptcl_size;
            no_cells_ch_A = ptr_ch_A->cell_size;
            //printf("no_cell_ch_A = %d\n", no_cells_ch_A);

            // printf("no_cells_ch_A = %d, no_ptcl_ch_A = %d\n", no_cells_ch_A, no_ptcl_ch_A);

            for (int cell_idx = 0; cell_idx < no_cells_ch_A; cell_idx += 8)
            {
                // if (ptr_ch_A->ID == 8)
                // {
                //     printf("ptr_ch_A->ptr_cell_idx[cell_idx], x = %d, y = %d, z = %d\n", ptr_ch_A->ptr_cell_idx_x[cell_idx], ptr_ch_A->ptr_cell_idx_y[cell_idx], ptr_ch_A->ptr_cell_idx_z[cell_idx]);
                // }

                // printf("\nno_cells ch A = %d, cell_idx ch A = %d\n", no_cells_ch_A, cell_idx);
                box_idx_x_node = (ptr_ch_A->ptr_cell_idx_x[cell_idx] >> 1) - ptr_node->box_ts_x;
                box_idx_y_node = (ptr_ch_A->ptr_cell_idx_y[cell_idx] >> 1) - ptr_node->box_ts_y;
                box_idx_z_node = (ptr_ch_A->ptr_cell_idx_z[cell_idx] >> 1) - ptr_node->box_ts_z;
                box_idx_node = box_idx_x_node + box_idx_y_node * ptr_node->box_real_dim_x + box_idx_z_node * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;
                // printf("before box_idx_node = %d, ptr_node->ptr_box_aux[box_idx_node] = %d\n", box_idx_node, ptr_node->ptr_box_aux[box_idx_node]);

                if (ptr_node->ptr_box_aux[box_idx_node] == new_zone_idx)
                {
                    box_idx_x_ch_A = ptr_ch_A->ptr_cell_idx_x[cell_idx] - ptr_ch_A->box_ts_x;
                    box_idx_y_ch_A = ptr_ch_A->ptr_cell_idx_y[cell_idx] - ptr_ch_A->box_ts_y;
                    box_idx_z_ch_A = ptr_ch_A->ptr_cell_idx_z[cell_idx] - ptr_ch_A->box_ts_z;
                    box_idx_ch_A = box_idx_x_ch_A + box_idx_y_ch_A * ptr_ch_A->box_real_dim_x + box_idx_z_ch_A * ptr_ch_A->box_real_dim_x * ptr_ch_A->box_real_dim_y;

                    if (ptr_ch_A->ptr_box[box_idx_ch_A] == -4)
                    {
                        for (int kk = 0; kk < 2; kk++)
                        {
                            for (int jj = 0; jj < 2; jj++)
                            {
                                for (int ii = 0; ii < 2; ii++)
                                {
                                    box_idxNbr_ch = box_idx_ch_A + ii + jj * ptr_ch_A->box_real_dim_x + kk * ptr_ch_A->box_real_dim_x * ptr_ch_A->box_real_dim_y;
                                    ptr_ch_A->ptr_box[box_idxNbr_ch] = -3; // Putting the status of EXIST (-3) in the child node cell
                                }
                            }
                        }
                    }

                }
                else 
                {
                    // printf("A2\n");
                    // printf("ptr_node->ptr_box_aux[box_idx_node] = %d\n", ptr_node->ptr_box_aux[box_idx_node]);
                    //  ptr_ch_B: Is the child node where the information from the child A need to be sent
                    // printf("id old child b = %d\n", links_old_ord_new[ptr_node->ptr_box_aux[box_idx_node]]);
                    ptr_ch_B = ptr_node->pptr_chn[links_old_ord_new[ptr_node->ptr_box_aux[box_idx_node]]];
                    // printf("ID child b = %d\n", ptr_ch_B->ID);
                    //  printf("A21\n");
                    //** >> CELLS **/
                    no_cells_ch_B = ptr_ch_B->cell_size;
                    // printf("no_cells_ch B = %d\n", no_cells_ch_B);
                    //  printf("A22\n");
                    for (int j = 0; j < 8; j++)
                    {
                        //Moving from child A to child B
                        ptr_ch_B->ptr_cell_idx_x[no_cells_ch_B + j] = ptr_ch_A->ptr_cell_idx_x[cell_idx + j];
                        ptr_ch_B->ptr_cell_idx_y[no_cells_ch_B + j] = ptr_ch_A->ptr_cell_idx_y[cell_idx + j];
                        ptr_ch_B->ptr_cell_idx_z[no_cells_ch_B + j] = ptr_ch_A->ptr_cell_idx_z[cell_idx + j];
                        //Moving from end of child A to removed cells of child A
                        ptr_ch_A->ptr_cell_idx_x[cell_idx + j] = ptr_ch_A->ptr_cell_idx_x[no_cells_ch_A - 8 + j];
                        ptr_ch_A->ptr_cell_idx_y[cell_idx + j] = ptr_ch_A->ptr_cell_idx_y[no_cells_ch_A - 8 + j];
                        ptr_ch_A->ptr_cell_idx_z[cell_idx + j] = ptr_ch_A->ptr_cell_idx_z[no_cells_ch_A - 8 + j];

                        //Updating the box status of child node B
                        box_idx_x_ch_B = ptr_ch_B->ptr_cell_idx_x[no_cells_ch_B + j] - ptr_ch_B->box_ts_x;
                        box_idx_y_ch_B = ptr_ch_B->ptr_cell_idx_y[no_cells_ch_B + j] - ptr_ch_B->box_ts_y;
                        box_idx_z_ch_B = ptr_ch_B->ptr_cell_idx_z[no_cells_ch_B + j] - ptr_ch_B->box_ts_z;
                        box_idx_ch_B = box_idx_x_ch_B + box_idx_y_ch_B * ptr_ch_B->box_real_dim_x + box_idx_z_ch_B * ptr_ch_B->box_real_dim_x * ptr_ch_B->box_real_dim_y;
                        ptr_ch_B->ptr_box[box_idx_ch_B] = -3; // Putting the status of EXIST (-3) in the child node cell
                    }

                    // if (ptr_ch_B->ID == 8)
                    // {
                    //     printf("ptr_ch_B->ptr_cell_idx[no_cells_ch_B]: x = %d, y = %d, z = %d\n", ptr_ch_B->ptr_cell_idx_x[no_cells_ch_B], ptr_ch_B->ptr_cell_idx_y[no_cells_ch_B], ptr_ch_B->ptr_cell_idx_z[no_cells_ch_B]);
                    // }
                    ptr_ch_B->cell_size += 8;
                    no_cells_ch_A -= 8;
                    cell_idx -= 8;

                    no_ptcl_ch_B = ptr_ch_B->ptcl_size;
                    //** >> MASS and PARTICLES **/
                    //printf("A3\n");
                    for (int j = 0; j < no_ptcl_ch_A; j++)
                    {
                        ptcl_idx = ptr_ch_A->ptr_ptcl[j];
                        box_idx_x_ptcl_node = GL_ptcl_x[ptcl_idx] * (1 << lv) - ptr_node->box_ts_x; // Particle indexes in the level
                        box_idx_y_ptcl_node = GL_ptcl_y[ptcl_idx] * (1 << lv) - ptr_node->box_ts_y;
                        box_idx_z_ptcl_node = GL_ptcl_z[ptcl_idx] * (1 << lv) - ptr_node->box_ts_z;
                        box_idx_ptcl_node = box_idx_x_ptcl_node + box_idx_y_ptcl_node * ptr_node->box_real_dim_x + box_idx_z_ptcl_node * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;
                        if (box_idx_ptcl_node == box_idx_node)
                        {

                            //** >> MASS **//

                            // Updating the mass status of child node B
                            box_idx_x_ptcl_ch_B = GL_ptcl_x[ptcl_idx] * (1 << (lv+1)) - ptr_ch_B->box_ts_x; // Particle indexes in the level
                            box_idx_y_ptcl_ch_B = GL_ptcl_y[ptcl_idx] * (1 << (lv+1)) - ptr_ch_B->box_ts_y;
                            box_idx_z_ptcl_ch_B = GL_ptcl_z[ptcl_idx] * (1 << (lv+1)) - ptr_ch_B->box_ts_z;
                            box_idx_ptcl_ch_B = box_idx_x_ptcl_ch_B + box_idx_y_ptcl_ch_B * ptr_ch_B->box_real_dim_x + box_idx_z_ptcl_ch_B * ptr_ch_B->box_real_dim_x * ptr_ch_B->box_real_dim_y;

                            //Box mass and local mass
                            ptr_ch_B->ptr_box_mass[box_idx_ptcl_ch_B] += GL_ptcl_mass[ptcl_idx];
                            ptr_ch_B->local_mass += GL_ptcl_mass[ptcl_idx];

                            // Updating the mass status of child node A
                            box_idx_x_ptcl_ch_A = GL_ptcl_x[ptcl_idx] * (1 << (lv + 1)) - ptr_ch_A->box_ts_x; // Particle indexes in the level
                            box_idx_y_ptcl_ch_A = GL_ptcl_y[ptcl_idx] * (1 << (lv + 1)) - ptr_ch_A->box_ts_y;
                            box_idx_z_ptcl_ch_A = GL_ptcl_z[ptcl_idx] * (1 << (lv + 1)) - ptr_ch_A->box_ts_z;
                            box_idx_ptcl_ch_A = box_idx_x_ptcl_ch_A + box_idx_y_ptcl_ch_A * ptr_ch_A->box_real_dim_x + box_idx_z_ptcl_ch_A * ptr_ch_A->box_real_dim_x * ptr_ch_A->box_real_dim_y;

                            // Box mass and local mass
                            ptr_ch_A->ptr_box_mass[box_idx_ptcl_ch_A] -= GL_ptcl_mass[ptcl_idx];
                            ptr_ch_A->local_mass -= GL_ptcl_mass[ptcl_idx];


                            //** >> PARTICLES **//

                            //** >> Space checking of the capacity of the particles in the child node B **/
                            if (space_check(&(ptr_ch_B->ptcl_cap), no_ptcl_ch_B + 1, "p1i1", &(ptr_ch_B->ptr_ptcl)) == _FAILURE_)
                            {
                                printf("Error, in space_check function\n");
                                return _FAILURE_;
                            }
                            ptr_ch_B->ptr_ptcl[no_ptcl_ch_B] = ptr_ch_A->ptr_ptcl[j];  // Adding the particle j of child node A to child node B
                            ptr_ch_A->ptr_ptcl[j] = ptr_ch_A->ptr_ptcl[no_ptcl_ch_A - 1]; // Replacing the particle j of the child node A to the last particle of the child node A
                            no_ptcl_ch_B++;                                               // The total number of particle increse
                            no_ptcl_ch_A--; // The total number of particle decrease
                            j--;       // The last element that was just moved to the current position should also must be analized
                        }
                    }
                    ptr_ch_B->ptcl_size = no_ptcl_ch_B;
                }
            }
            ptr_ch_A->cell_size = no_cells_ch_A;
            ptr_ch_A->ptcl_size = no_ptcl_ch_A;
        }
    }

    //printf("B\n");

    //Case of old child nodes that will not be used
    cntr_ch = 0;
    for (int ch = 0; ch < ptr_node->chn_size - ptr_node->zones_size; ch++)
    {
        while (links_old_ord_old[cntr_ch] == cntr_ch)
        {
            cntr_ch++;
        }

        ptr_ch_A = ptr_node->pptr_chn[cntr_ch];

        no_ptcl_ch_A = ptr_ch_A->ptcl_size;
        no_cells_ch_A = ptr_ch_A->cell_size;
        for (int cell_idx = 0; cell_idx < no_cells_ch_A; cell_idx += 8)
        {
            box_idx_x_node = (ptr_ch_A->ptr_cell_idx_x[cell_idx] >> 1) - ptr_node->box_ts_x;
            box_idx_y_node = (ptr_ch_A->ptr_cell_idx_y[cell_idx] >> 1) - ptr_node->box_ts_y;
            box_idx_z_node = (ptr_ch_A->ptr_cell_idx_z[cell_idx] >> 1) - ptr_node->box_ts_z;
            box_idx_node = box_idx_x_node + box_idx_y_node * ptr_node->box_real_dim_x + box_idx_z_node * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;
            if (ptr_node->ptr_box_aux[box_idx_node] >= 0)
            {
                // ptr_ch_B: Is the child node where the information from the child A need to be sent
                ptr_ch_B = ptr_node->pptr_chn[links_old_ord_new[ptr_node->ptr_box_aux[box_idx_node]]];

                //** >> CELLS **/
                no_cells_ch_B = ptr_ch_B->cell_size;
                for (int j = 0; j < 8; j++)
                {
                    // Moving from child A to child B
                    ptr_ch_B->ptr_cell_idx_x[no_cells_ch_B + j] = ptr_ch_A->ptr_cell_idx_x[cell_idx + j];
                    ptr_ch_B->ptr_cell_idx_y[no_cells_ch_B + j] = ptr_ch_A->ptr_cell_idx_y[cell_idx + j];
                    ptr_ch_B->ptr_cell_idx_z[no_cells_ch_B + j] = ptr_ch_A->ptr_cell_idx_z[cell_idx + j];

                    // Updating the box status of child node B
                    box_idx_x_ch_B = ptr_ch_B->ptr_cell_idx_x[no_cells_ch_B + j] - ptr_ch_B->box_ts_x;
                    box_idx_y_ch_B = ptr_ch_B->ptr_cell_idx_y[no_cells_ch_B + j] - ptr_ch_B->box_ts_y;
                    box_idx_z_ch_B = ptr_ch_B->ptr_cell_idx_z[no_cells_ch_B + j] - ptr_ch_B->box_ts_z;
                    box_idx_ch_B = box_idx_x_ch_B + box_idx_y_ch_B * ptr_ch_B->box_real_dim_x + box_idx_z_ch_B * ptr_ch_B->box_real_dim_x * ptr_ch_B->box_real_dim_y;
                    ptr_ch_B->ptr_box[box_idx_ch_B] = -3; // Putting the status of EXIST (-3) in the child node cell
                }
                ptr_ch_B->cell_size += 8;

                no_ptcl_ch_B = ptr_ch_B->ptcl_size;
                //** >> MASS and PARTICLES **/
                for (int j = 0; j < no_ptcl_ch_A; j++)
                {
                    ptcl_idx = ptr_ch_A->ptr_ptcl[j];
                    box_idx_x_ptcl_node = GL_ptcl_x[ptcl_idx] * (1 << lv) - ptr_node->box_ts_x; // Particle indexes in the level
                    box_idx_y_ptcl_node = GL_ptcl_y[ptcl_idx] * (1 << lv) - ptr_node->box_ts_y;
                    box_idx_z_ptcl_node = GL_ptcl_z[ptcl_idx] * (1 << lv) - ptr_node->box_ts_z;
                    box_idx_ptcl_node = box_idx_x_ptcl_node + box_idx_y_ptcl_node * ptr_node->box_real_dim_x + box_idx_z_ptcl_node * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;
                    if (box_idx_ptcl_node == box_idx_node)
                    {

                        //** >> MASS **//

                        // Updating the mass status of child node B
                        box_idx_x_ptcl_ch_B = GL_ptcl_x[ptcl_idx] * (1 << (lv + 1)) - ptr_ch_B->box_ts_x; // Particle indexes in the level
                        box_idx_y_ptcl_ch_B = GL_ptcl_y[ptcl_idx] * (1 << (lv + 1)) - ptr_ch_B->box_ts_y;
                        box_idx_z_ptcl_ch_B = GL_ptcl_z[ptcl_idx] * (1 << (lv + 1)) - ptr_ch_B->box_ts_z;
                        box_idx_ptcl_ch_B = box_idx_x_ptcl_ch_B + box_idx_y_ptcl_ch_B * ptr_ch_B->box_real_dim_x + box_idx_z_ptcl_ch_B * ptr_ch_B->box_real_dim_x * ptr_ch_B->box_real_dim_y;

                        // Box mass and local mass
                        ptr_ch_B->ptr_box_mass[box_idx_ptcl_ch_B] += GL_ptcl_mass[ptcl_idx];
                        ptr_ch_B->local_mass += GL_ptcl_mass[ptcl_idx];

                        //** >> PARTICLES **//

                        //** >> Space checking of the capacity of the particles in the child node B **/
                        if (space_check(&(ptr_ch_B->ptcl_cap), no_ptcl_ch_B + 1, "p1i1", &(ptr_ch_B->ptr_ptcl)) == _FAILURE_)
                        {
                            printf("Error, in space_check function\n");
                            return _FAILURE_;
                        }
                        ptr_ch_B->ptr_ptcl[no_ptcl_ch_B] = ptcl_idx;                  // Adding the particle j of child node A to child node B
                        ptr_ch_A->ptr_ptcl[j] = ptr_ch_A->ptr_ptcl[no_ptcl_ch_A - 1]; // Replacing the particle j of the child node A to the last particle of the child node A
                        no_ptcl_ch_B++;                                               // The total number of particle increse
                        no_ptcl_ch_A--;                                            // The total number of particle decrease
                        j--;                                                       // The last element that was just moved to the current position should also must be analized
                    }
                }
                ptr_ch_B->ptcl_size = no_ptcl_ch_B;
            }
        }
        cntr_ch++;
    }

    ptr_ch_A = NULL;
    ptr_ch_B = NULL;

    return _SUCCESS_;
}

static int moving_new_zones_to_new_child(struct node *ptr_node, int *links_old_ord_new)
{
    struct node *ptr_ch = NULL; // child node

    int lv;

    int box_idx_x_node; // Box index in X direcction of the node cell
    int box_idx_y_node; // Box index in Y direcction of the node cell
    int box_idx_z_node; // Box index in Z direcction of the node cell
    int box_idx_node;   // Box index of the node cell

    int box_idx_x_ch; // Box index in X direcction of the node cell
    int box_idx_y_ch; // Box index in Y direcction of the node cell
    int box_idx_z_ch; // Box index in Z direcction of the node cell
    int box_idx_ch;   // Box index of the node cell

    int box_idx_x_ptcl_node; // Box index in X direcction of the particles in the node
    int box_idx_y_ptcl_node; // Box index in Y direcction of the particles in the node
    int box_idx_z_ptcl_node; // Box index in Z direcction of the particles in the node
    int box_idx_ptcl_node;   // Box index of particles in the node

    int box_idx_x_ptcl_ch; // Box index in X direcction of the particles in the child node 
    int box_idx_y_ptcl_ch; // Box index in Y direcction of the particles in the child node
    int box_idx_z_ptcl_ch; // Box index in Z direcction of the particles in the child node
    int box_idx_ptcl_ch;   // Box index of particles in the child node

    int cell_idx_x_aux;
    int cell_idx_y_aux;
    int cell_idx_z_aux;

    int cell_idx_node;
    int cell_idx_ch;

    int ptcl_idx;

    int no_ptcl_node; // Number of particles in the node
    int no_cells_ch;    //Number of cell in the child node
    int no_ptcl_ch; //Number of particles in the child node

    lv = ptr_node->lv;

    no_ptcl_node = ptr_node->ptcl_size;
    for (int zone_idx = 0; zone_idx < ptr_node->zones_size; zone_idx++)
    {
        ptr_ch = ptr_node->pptr_chn[links_old_ord_new[zone_idx]];
        no_cells_ch = ptr_ch->cell_size;
        no_ptcl_ch = ptr_ch->ptcl_size;
        for (int zone_element = 0; zone_element < ptr_node->ptr_zone_size[zone_idx]; zone_element++)
        {
            cell_idx_node = ptr_node->pptr_zones[zone_idx][zone_element];
            box_idx_x_node = ptr_node->ptr_cell_idx_x[cell_idx_node] - ptr_node->box_ts_x;
            box_idx_y_node = ptr_node->ptr_cell_idx_y[cell_idx_node] - ptr_node->box_ts_y;
            box_idx_z_node = ptr_node->ptr_cell_idx_z[cell_idx_node] - ptr_node->box_ts_z;
            box_idx_node = box_idx_x_node + box_idx_y_node * ptr_node->box_real_dim_x + box_idx_z_node * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;

            // Case cell is a new refinement cell. If it was not new, the past child nodes would have already added it to the corresponding child node
            if(ptr_node->ptr_box[box_idx_node] < 0)
            {
                //** >> CELLS AND BOX STATUS**/
                cell_idx_x_aux = ptr_node->ptr_cell_idx_x[cell_idx_node] * 2;
                cell_idx_y_aux = ptr_node->ptr_cell_idx_y[cell_idx_node] * 2;
                cell_idx_z_aux = ptr_node->ptr_cell_idx_z[cell_idx_node] * 2;
                for (int kk = 0; kk < 2; kk++)
                {
                    for (int jj = 0; jj < 2; jj++)
                    {
                        for (int ii = 0; ii < 2; ii++)
                        {
                            cell_idx_ch = no_cells_ch + ii + jj * 2 + kk * 4;
                            //Cells
                            ptr_ch->ptr_cell_idx_x[cell_idx_ch] = cell_idx_x_aux + ii;
                            ptr_ch->ptr_cell_idx_y[cell_idx_ch] = cell_idx_y_aux + jj;
                            ptr_ch->ptr_cell_idx_z[cell_idx_ch] = cell_idx_z_aux + kk;

                            //Box status
                            box_idx_x_ch = ptr_ch->ptr_cell_idx_x[cell_idx_ch] - ptr_ch->box_ts_x;
                            box_idx_y_ch = ptr_ch->ptr_cell_idx_y[cell_idx_ch] - ptr_ch->box_ts_y;
                            box_idx_z_ch = ptr_ch->ptr_cell_idx_z[cell_idx_ch] - ptr_ch->box_ts_z;
                            box_idx_ch = box_idx_x_ch + box_idx_y_ch * ptr_ch->box_real_dim_x + box_idx_z_ch * ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;
                            ptr_ch->ptr_box[box_idx_ch] = -3;
                        }
                    }
                }
                no_cells_ch += 8;
                //** >> MASS and PARTICLES **/
                //Cycle over all particles in the node
                for (int j = 0; j < no_ptcl_node; j++)
                {
                    ptcl_idx = ptr_node->ptr_ptcl[j];
                    box_idx_x_ptcl_node = GL_ptcl_x[ptcl_idx] * (1 << lv) - ptr_node->box_ts_x; // Particle indexes in the level
                    box_idx_y_ptcl_node = GL_ptcl_y[ptcl_idx] * (1 << lv) - ptr_node->box_ts_y;
                    box_idx_z_ptcl_node = GL_ptcl_z[ptcl_idx] * (1 << lv) - ptr_node->box_ts_z;
                    box_idx_ptcl_node = box_idx_x_ptcl_node + box_idx_y_ptcl_node * ptr_node->box_real_dim_x + box_idx_z_ptcl_node * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;
                    if (box_idx_ptcl_node == box_idx_node)
                    {

                        //** >> MASS **//
                        // Updating the mass status of child node
                        box_idx_x_ptcl_ch = GL_ptcl_x[ptcl_idx] * (1 << (lv + 1)) - ptr_ch->box_ts_x; // Particle indexes in the level
                        box_idx_y_ptcl_ch = GL_ptcl_y[ptcl_idx] * (1 << (lv + 1)) - ptr_ch->box_ts_y;
                        box_idx_z_ptcl_ch = GL_ptcl_z[ptcl_idx] * (1 << (lv + 1)) - ptr_ch->box_ts_z;
                        box_idx_ptcl_ch = box_idx_x_ptcl_ch + box_idx_y_ptcl_ch * ptr_ch->box_real_dim_x + box_idx_z_ptcl_ch * ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;
                        // Box mass and local mass
                        ptr_ch->ptr_box_mass[box_idx_ptcl_ch] += GL_ptcl_mass[ptcl_idx];
                        ptr_ch->local_mass += GL_ptcl_mass[ptcl_idx];
                        //** >> PARTICLES **//

                        //** >> Space checking of the capacity of the particles in the child node B **/
                        if (space_check(&(ptr_ch->ptcl_cap), no_ptcl_ch + 1, "p1i1", &(ptr_ch->ptr_ptcl)) == _FAILURE_)
                        {
                            printf("Error, in space_check function\n");
                            return _FAILURE_;
                        }
                        ptr_ch->ptr_ptcl[no_ptcl_ch] = ptcl_idx;                   // Adding the particle j of child node to child node
                        no_ptcl_ch++;                                               // The total number of particle increse
                    }
                }
            }
        }
        ptr_ch->cell_size = no_cells_ch;
        ptr_ch->ptcl_size = no_ptcl_ch;
    }

    ptr_ch = NULL;

    return _SUCCESS_;
}

static void update_border_child_boxes(struct node *ptr_node, const int *links_old_ord_old)
{

    struct node *ptr_ch = NULL;

    int box_idx_x_ch; // Box index in X direcction of the child cell
    int box_idx_y_ch; // Box index in Y direcction of the child cell
    int box_idx_z_ch; // Box index in Z direcction of the child cell
    int box_idx_ch;   // Box index of the child cell

    int box_idxNbr_ch; // Box index in the neigborhood

    for (int zone_idx = 0; zone_idx < ptr_node->zones_size; zone_idx++)
    {
        ptr_ch = ptr_node->pptr_chn[links_old_ord_old[zone_idx]];

        // Changing the EXIST (-3) TO BORDER (-2) to all border cells in the box
        for (int j = 0; j < ptr_ch->cell_size; j++)
        {
            box_idx_x_ch = ptr_ch->ptr_cell_idx_x[j] - ptr_ch->box_ts_x;
            box_idx_y_ch = ptr_ch->ptr_cell_idx_y[j] - ptr_ch->box_ts_y;
            box_idx_z_ch = ptr_ch->ptr_cell_idx_z[j] - ptr_ch->box_ts_z;
            box_idx_ch = box_idx_x_ch + box_idx_y_ch * ptr_ch->box_real_dim_x + box_idx_z_ch * ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;

            for (int kk = -1; kk < 2; kk++)
            {
                for (int jj = -1; jj < 2; jj++)
                {
                    for (int ii = -1; ii < 2; ii++)
                    {
                        box_idxNbr_ch = box_idx_ch + ii + jj * ptr_ch->box_real_dim_x + kk * ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;
                        if (ptr_ch->ptr_box[box_idxNbr_ch] == -4) // Border cells are those such that at least on of their first neighbors are NO-EXIST cells.
                        {
                            ptr_ch->ptr_box[box_idx_ch] = -2;
                            ii = 2;
                            jj = 2;
                            kk = 2;
                        }
                    }
                }
            }
        }
    }

    ptr_ch = NULL;
}

static int reset_grid_properties(struct node *ptr_node, const int *links_old_ord_old)
{
    struct node *ptr_ch = NULL; // child node

    int size;

    // Cycle over new refinement zones
    for (int zone_idx = 0; zone_idx < ptr_node->zones_size; zone_idx++)
    {
        if (links_old_ord_old[zone_idx] < ptr_node->chn_size)
        {
            ptr_ch = ptr_node->pptr_chn[links_old_ord_old[zone_idx]];

            //* Potential, Acceleration and density capacity of the grid **/
            size = (ptr_ch->box_real_dim_x + 1) * (ptr_ch->box_real_dim_y + 1) * (ptr_ch->box_real_dim_z + 1);
            //** >> Space checking
            if (space_check(&(ptr_ch->grid_properties_cap), size, "p5v1v1v1v1v1", &(ptr_ch->ptr_pot), &(ptr_ch->ptr_ax), &(ptr_ch->ptr_ay), &(ptr_ch->ptr_az), &(ptr_ch->ptr_d)) == _FAILURE_)
            {
                printf("Error, in space_check function\n");
                return _FAILURE_;
            }

            //** >> Density initialization **/
            for (int j = 0; j < size; j++)
            {
                ptr_ch->ptr_d[j] = 0;
            }
        }
    }

    ptr_ch = NULL;

    return _SUCCESS_;
}

static void reorganization_child_node(struct node *ptr_node, const int *links_old_ord_old, const int *links_new_ord_old,  const int *links_old_ord_new)
{
    printf("A00\n");

    struct node *ptr_ch = NULL;
    struct node **pptr_aux = NULL;
    struct node **pptr_aux_free = NULL;

    int cntr_ch;


    // if(chn_cap > 0)
    // {
    //     pptr_aux = (struct node **)malloc(ptr_node->chn_cap * sizeof(struct node *));
    // }

    pptr_aux = (struct node **) malloc(ptr_node->chn_cap * sizeof(struct node *));

    // Updating the ID of the child nodes
    for (int zone_idx = 0; zone_idx < ptr_node->zones_size; zone_idx++)
    {
        ptr_ch = ptr_node->pptr_chn[links_old_ord_old[zone_idx]];
        ptr_ch->ID = links_new_ord_old[zone_idx];
    }

    //Reorganization of child nodes
    for (int i = 0; i < ptr_node->zones_size; i++)
    {
        pptr_aux[i] = ptr_node->pptr_chn[links_old_ord_new[i]];
    }

    cntr_ch = 0;
    for (int i = 0; i < ptr_node->chn_size - ptr_node->zones_size; i++)
    {
        while (links_old_ord_old[cntr_ch] == cntr_ch)
        {
            cntr_ch++;
        }
        pptr_aux[ptr_node->zones_size + i] = ptr_node->pptr_chn[cntr_ch];
    }

    pptr_aux_free = ptr_node->pptr_chn;
    ptr_node->pptr_chn = pptr_aux;

    pptr_aux = NULL;
    free(pptr_aux_free);
    pptr_aux_free = NULL;
    ptr_ch = NULL;
}

static int reorganization_grandchild_node(struct node *ptr_node)
{
     

    int no_grandchildren = 0;

    //Computing the total number of granchildren and reseting its value to 0
    for (int ch = 0; ch < ptr_node->chn_size; ch++)
    {
        no_grandchildren += ptr_node->pptr_chn[ch]->chn_size;
        ptr_node->pptr_chn[ch]->chn_size = 0;   //Reseting the number of grandchildren
    }

    if(no_grandchildren > 0)
    {
        int box_idx_x_node; // Box index in X direcction of the node cell
        int box_idx_y_node; // Box index in Y direcction of the node cell
        int box_idx_z_node; // Box index in Z direcction of the node cell
        int box_idx_node;   // Box index of the node cell

        int child_ID;

        int size;

        struct node **pptr_aux = NULL;
        int cntr_ch = 0; //Counter the number of grandchildren
        pptr_aux = (struct node **)malloc(no_grandchildren * sizeof(struct node *));

        //Storing the grandchildren in the auxiliary array pptr_aux
        for (int ch = 0; ch < ptr_node->chn_size; ch++)
        {
            for (int grandch = 0; grandch < ptr_node->pptr_chn[ch]->chn_size; grandch++)
            {
                pptr_aux[cntr_ch] = ptr_node->pptr_chn[ch]->pptr_chn[grandch];
                cntr_ch++;
            }
        }

        for (int grandch = 0; grandch < no_grandchildren; grandch++)
        {
            box_idx_x_node = (pptr_aux[grandch]->ptr_cell_idx_x[0] >> 2) - ptr_node->box_ts_x;
            box_idx_y_node = (pptr_aux[grandch]->ptr_cell_idx_y[0] >> 2) - ptr_node->box_ts_y;
            box_idx_z_node = (pptr_aux[grandch]->ptr_cell_idx_z[0] >> 2) - ptr_node->box_ts_z;
            box_idx_node = box_idx_x_node + box_idx_y_node * ptr_node->box_real_dim_x + box_idx_z_node * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;
            child_ID = ptr_node->ptr_box_aux[box_idx_node];

            size = ptr_node->pptr_chn[child_ID]->chn_size;
            //** >> Space checking of the child capacity of the child node "child_ID" **/
            if (space_check(&(ptr_node->pptr_chn[child_ID]->chn_cap), size + 1, "p1n2", &(ptr_node->pptr_chn[child_ID]->pptr_chn)) == _FAILURE_)
            {
                printf("Error, in space_check function\n");
                return _FAILURE_;
            }

            ptr_node->pptr_chn[child_ID]->pptr_chn[size] = pptr_aux[grandch];
            ptr_node->pptr_chn[child_ID]->chn_size = size + 1;
        }

        free(pptr_aux);
        pptr_aux = NULL;
    }


    return _SUCCESS_;
}

static int tentacles_updating(struct node *ptr_node, int tentacle_lv)
{
    int no_tentacles = GL_tentacles_size[tentacle_lv];
    int size = no_tentacles + ptr_node->zones_size;

    //** >> Space checking of the capacity of the refined cells **/
    if (space_check(&(GL_tentacles_cap[tentacle_lv]), size, "p1n2", &(GL_tentacles_new[tentacle_lv])) == _FAILURE_)
    {
        printf("Error, in space_check function\n");
        return _FAILURE_;
    }

    for (int i = 0; i < ptr_node->zones_size; i++)
    {
        GL_tentacles_new[tentacle_lv][no_tentacles + i] = ptr_node->pptr_chn[i];
    }
    GL_tentacles_size[tentacle_lv] = size;

    return _SUCCESS_;
}

static void updating_tentacles_max_lv()
{
    int no_lvs = GL_tentacles_level_max < lmax - lmin ? GL_tentacles_level_max : GL_tentacles_level_max - 1;

    for (int lv = no_lvs + 1; lv > -1; lv--)
    {
        if(GL_tentacles_size[lv] > 0)
        {
            GL_tentacles_level_max = lv;
            lv = -1;
        }
    }
}

int tree_adaptation()
{

    //** >> Working in the refinement zones **/
    if (lmin < lmax)
    {
        struct node *ptr_node = NULL;
        int no_pts; // Number of parents in the cycle
        int no_lvs; // Number of level of refinement to adapt

        int *links_old_ord_old = NULL; // Storing the old zone of refinement id using old order
        int *links_new_ord_old = NULL; // Storing the new zone of refinement id using old order
        int *links_old_ord_new = NULL; // Storing the old zone of refinement id using new order
        int *links_new_ord_new = NULL; // Storing the new zone of refinement id using new order
        int links_cap;

        links_cap = 64;
        links_old_ord_old = (int *) malloc (64 * sizeof(int));
        links_new_ord_old = (int *) malloc (64 * sizeof(int));
        links_old_ord_new = (int *) malloc (64 * sizeof(int));
        links_new_ord_new = (int *) malloc (64 * sizeof(int));

        no_lvs = GL_tentacles_level_max < lmax - lmin ? GL_tentacles_level_max : GL_tentacles_level_max - 1;

        for (int lv = no_lvs; lv > -1; lv--)
        {
            GL_tentacles_size[lv + 1] = 0;
            no_pts = GL_tentacles_size[lv];
            //** >> For cycle over parent nodes **/
            for (int i = 0; i < no_pts; i++)
            {
                ptr_node = GL_tentacles_old[lv][i];
                printf("\nLv = %d, parent = %d\n", lv, i);

                //** Updating the box mass information **/
                printf("Updating box mass\n");
                updating_box_mass(ptr_node);

                //** Initialization of the box_aux **/
                printf("Initialization box aux\n");
                initialization_box_aux(ptr_node);

                //** Initialization of the auiliary refinement arrays**/
                printf("Initialization ref aux\n");
                initialization_ref_aux(ptr_node);

                //** >> Filling the refinement cells array **/
                printf("Fill cell ref\n");
                if (fill_cell_ref(ptr_node) == _FAILURE_)
                {
                    printf("Error at function fill_cell_ref()\n");
                    return _FAILURE_;
                }

                //** >> Filling the different zones of refinement **/
                printf("Fill zones ref\n");
                if (fill_zones_ref(ptr_node) == _FAILURE_)
                {
                    printf("Error at function fill_zones_ref()\n");
                    return _FAILURE_;
                }

                //** >> Adapting to new refinement zones **/
                printf("Create links\n");
                if (create_links(ptr_node, links_old_ord_old, links_new_ord_old, links_old_ord_new, links_new_ord_new, &links_cap) == _FAILURE_)
                {
                    printf("Error at function fill_zones_ref()\n");
                    return _FAILURE_;
                }

                //** >> Adapting child boxes to the new space **/
                printf("Adapt child box\n");
                if (0 < ptr_node->zones_size && 0 < ptr_node->chn_size)
                {
                    if (adapt_child_box_and_cells(ptr_node, links_old_ord_old, links_new_ord_old) == _FAILURE_)
                    {
                        printf("Error at function adapt_child_box()\n");
                        return _FAILURE_;
                    }
                }

                //** >> Creating and defining new children nodes for excess in refinement zones and and linking them to the parent node ptr_node **/
                printf("Create new child nodes\n");
                if (ptr_node->zones_size > ptr_node->chn_size)
                {

                    if (create_new_child_nodes(ptr_node, links_old_ord_old, links_new_ord_old) == _FAILURE_)
                    {
                        printf("Error, in new_child_nodes function\n");
                        return _FAILURE_;
                    }
                }

                //** >> Removing cells that no longer require refinement **/
                printf("Removing cells nolonger requiere refinement\n");
                if (0 < ptr_node->zones_size && 0 < ptr_node->chn_size)
                {
                    remov_cells_nolonger_require_refinement(ptr_node, links_old_ord_old, links_new_ord_old);
                }

                printf("chn size = %d, zones size = %d, chn cap = %d\n", ptr_node->chn_size, ptr_node->zones_size, ptr_node->chn_cap);
                printf("links_old_ord_old = [ ");
                for (int i = 0; i < ptr_node->zones_size; i++)
                {
                    printf("%d ", links_old_ord_old[i]);
                }
                printf("]\n");
                printf("links_new_ord_old = [ ");
                for (int i = 0; i < ptr_node->zones_size; i++)
                {
                    printf("%d ", links_new_ord_old[i]);
                }
                printf("]\n");

                //** >> Adapting the information from old child nodes to new child nodes **/
                printf("Moving old child to new child\n");
                if (0 < ptr_node->zones_size && 0 < ptr_node->chn_size)
                {
                    if (moving_old_child_to_new_child(ptr_node, links_old_ord_old, links_new_ord_old, links_old_ord_new, links_new_ord_new) == _FAILURE_)
                    {
                        printf("Error at function adapt_old_child_cells_to_new_child()\n");
                        return _FAILURE_;
                    }
                }

                //** >> Moving new zones of refienemnt information to all child nodes **/
                printf("Moving new zones to new childd\n");
                if (0 < ptr_node->zones_size)
                {
                    if (moving_new_zones_to_new_child(ptr_node, links_old_ord_new) == _FAILURE_)
                    {
                        printf("Error at function adapt_old_child_cells_to_new_child()\n");
                        return _FAILURE_;
                    }
                }

                

                //** >> Update border of the child boxes **//
                printf("Update border child boxes\n");
                if (ptr_node->zones_size > 0)
                {
                    update_border_child_boxes(ptr_node, links_old_ord_old);
                }
                    

                //** >> Reseting the density of the reused child nodes **/
                printf("Reset grid properties\n");
                if (0 < ptr_node->zones_size)
                {
                    if (reset_grid_properties(ptr_node, links_old_ord_old) == _FAILURE_)
                    {
                        printf("Error at function reset_grid_properties()\n");
                        return _FAILURE_;
                    }
                }

                //** >> Reorganization child nodes **/
                printf("Reorganization child nodes\n");
                if(ptr_node->zones_size > 0)
                {
                    reorganization_child_node(ptr_node, links_old_ord_old, links_new_ord_old, links_old_ord_new);
                }

                //** >> Reorganization grandchild nodes **/
                printf("Reorganization grandchild nodes\n");
                if (ptr_node->zones_size > 0 && ptr_node->chn_size > 0 && lv < no_lvs)
                {
                    if (reorganization_grandchild_node(ptr_node) == _FAILURE_)
                    {
                        printf("Error at function reorganization_grandchild_node()\n");
                        return _FAILURE_;
                    }
                }

                //** Tentacles updating **/
                printf("Tentacles Updating\n");
                if (0 < ptr_node->zones_size)
                {
                    if (tentacles_updating(ptr_node, lv + 1) == _FAILURE_)
                    {
                        printf("Error at function tentacles_updating()\n");
                        return _FAILURE_;
                    }
                }  
            }
        }

        //** >> Tentacles Updating lv max **/
        printf("Tentacles updating lv max\n");
        updating_tentacles_max_lv();

        free(links_old_ord_old);
        free(links_new_ord_old);
        free(links_old_ord_new);
        free(links_new_ord_new);
        links_old_ord_old = NULL;
        links_new_ord_old = NULL;
        links_old_ord_new = NULL;
        links_new_ord_new = NULL;
        ptr_node = NULL;
    }

    return _SUCCESS_;
}