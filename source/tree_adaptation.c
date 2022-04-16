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
    struct node *ptr_node_ch; //child node
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
        ptr_node_ch = ptr_node->pptr_chn[i];
        for (int j = 0; j < ptr_node_ch->cell_size; j += 8) // Cycle over packeges of 8 cells
        {
            //** >> Computing the mass of the 8 child cells **/
            box_idx_x_ch = ptr_node_ch->ptr_cell_idx_x[j] - ptr_node_ch->box_ts_x;
            box_idx_y_ch = ptr_node_ch->ptr_cell_idx_y[j] - ptr_node_ch->box_ts_y;
            box_idx_z_ch = ptr_node_ch->ptr_cell_idx_z[j] - ptr_node_ch->box_ts_z;

            aux_box_mass = 0;

            for (int kk = 0; kk < 2; kk++)
            {
                for (int jj = 0; jj < 2; jj++)
                {
                    for (int ii = 0; ii < 2; ii++)
                    {
                        box_idx_ch = (box_idx_x_ch + ii) + (box_idx_y_ch + jj) * ptr_node_ch->box_real_dim_x + (box_idx_z_ch + kk) * ptr_node_ch->box_real_dim_x * ptr_node_ch->box_real_dim_y;
                        aux_box_mass += ptr_node_ch->ptr_box_mass[box_idx_ch];
                    }
                }
            }
            //** Re-defining the mass of the node cell corresponding to those 8 child cells
            box_idx_x_node = (ptr_node_ch->ptr_cell_idx_x[j] >> 1) - ptr_node->box_ts_x;
            box_idx_y_node = (ptr_node_ch->ptr_cell_idx_y[j] >> 1) - ptr_node->box_ts_y;
            box_idx_z_node = (ptr_node_ch->ptr_cell_idx_z[j] >> 1) - ptr_node->box_ts_z;
            box_idx_node = box_idx_x_node + box_idx_y_node * ptr_node->box_real_dim_x + box_idx_z_node * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;
            ptr_node->ptr_box_mass[box_idx_node] = aux_box_mass;
            total_mass += aux_box_mass;
        }
    }
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

    struct node *ptr_ch; // Child node

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

    // //** >> Space checking of the capacity of the refined cells **/
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

                //** >> Initiazling the new allocated values **/
                for (int i = 0; i < zone_idx_max - zone_idx; i++)
                {
                    ptr_node->pptr_zones[zone_idx_max - i - 1] = NULL;
                    ptr_node->ptr_zone_cap[zone_idx_max - i - 1] = 0;
                }
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

static int adapt_new_ref_zones(struct node *ptr_node)
{
    struct node *ptr_node_ch; // child node

    int cntr_links; // Counter the number of links between the new and old zones of refinement
    //int cntr_old_ref_cells; // Counter the total number of old cell refined in the node
    int cntr_cell_ch; //Counter the cells of the child node
    int cntr_old_ref_zones;             // total number of zones refined in the old box
    int cntr_old_links_plus;      // Counte links old
    int cntr_link_elements;         // Counter linked elements
    bool check_link; //Check if the element is linked
    int links_old[ptr_node->zones_size];    // Storing the old zone of refinement id
    int links_new[ptr_node->zones_size];    // Storing the new zone of refinement id
    //int cell_idx_ch; // Index of the child cell
    int box_value_new; // Value in the new box

    bool create_link; //Decide if the link need to be created

    int box_idx_x; // Box index at X direction
    int box_idx_y; // Box index at Y direction 
    int box_idx_z; // Box index at Z direction 
    int box_idx; // Box index 

    //** >> Creating the links between old and new refinement zones IDs **/

    cntr_links = 0;
    cntr_old_ref_zones = 0;

    for (int i = 0; i < ptr_node->zones_size; i++)
    {
        links_new[i] = -1;
        links_old[i] = -1;
    }

    while (cntr_links < ptr_node->zones_size && cntr_old_ref_zones < ptr_node->chn_size)
    {
        ptr_node_ch = ptr_node->pptr_chn[cntr_old_ref_zones];
        cntr_cell_ch = 0;
        while (cntr_cell_ch < ptr_node_ch->cell_size)
        {
            box_idx_x = (ptr_node_ch->ptr_cell_idx_x[cntr_cell_ch] >> 1) - ptr_node->box_ts_x;
            box_idx_y = (ptr_node_ch->ptr_cell_idx_y[cntr_cell_ch] >> 1) - ptr_node->box_ts_y;
            box_idx_z = (ptr_node_ch->ptr_cell_idx_z[cntr_cell_ch] >> 1) - ptr_node->box_ts_z;
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
                    if (box_value_new == links_new[i]) // The link already exist
                    {
                        create_link = false;
                        i = cntr_links;
                        cntr_cell_ch += 8;
                    }
                }
                if (create_link == true)
                {
                    links_new[cntr_links] = box_value_new;
                    links_old[cntr_links] = cntr_old_ref_zones;
                    cntr_links++;
                    cntr_cell_ch = ptr_node_ch->cell_size; //Breiking the while cycle
                }
            }
        }//End while cycle over the child node cells
        cntr_old_ref_zones++;
    } // End while cycle over all old refinement zones

    // printf("\nlv = %d, parent = %d\n", ptr_node->lv,ptr_node->ID);
    // printf("zones new = %d, zones old = %d\n", ptr_node->zones_size, ptr_node->chn_size);
    // printf("cntr_links = %d\n", cntr_links);

    if (cntr_links < ptr_node->zones_size)
    {
        
        //** >> Old links **/
        if (cntr_links == ptr_node->chn_size)
        {
            for (int i = cntr_links; i < ptr_node->zones_size; i++)
            {
                //ptr_node->ptr_aux_idx[i - cntr_links] = i;
                links_old[i] = i;
            }
            //cntr_links = ptr_node->zones_size;
        }
        else
        {
            cntr_link_elements = 0;  // Counter linked elements
            cntr_old_links_plus = 0; // Counte links old

            while (cntr_link_elements < cntr_links && cntr_links + cntr_old_links_plus < ptr_node->zones_size && cntr_links + cntr_old_links_plus < ptr_node->chn_size)
            {
                if (cntr_link_elements + cntr_old_links_plus == links_old[cntr_link_elements])
                {
                    cntr_link_elements++;
                }
                else
                {
                    links_old[cntr_links + cntr_old_links_plus] = cntr_link_elements + cntr_old_links_plus;
                    cntr_old_links_plus++;
                }

            }

            if (cntr_links + cntr_old_links_plus < ptr_node->zones_size)
            {
                for (int i = cntr_links + cntr_old_links_plus; i < ptr_node->zones_size; i++)
                {
                    // ptr_node->ptr_aux_idx[i - cntr_links] = i;
                    links_old[i] = i;
                }
            }
        }

        //** >> New links **/
        cntr_link_elements = -1;
        for (int i = cntr_links; i < ptr_node->zones_size; i++)
        {
            check_link = true;
            while(check_link == true)
            {
                check_link = false;
                cntr_link_elements++;
                for (int j = 0; j < cntr_links; j++)
                {
                    if (links_new[j] == cntr_link_elements)
                    {
                        check_link = true;
                        j = cntr_links;
                    }
                }
            }
            links_new[i] = cntr_link_elements;
        }
    }

    // printf("Links_old = [   ");
    // for (int i = 0; i< ptr_node->zones_size; i++)
    // {
    //     printf("%d  ", links_old[i]);
    // }
    // printf("]\n");
    // printf("Links_new = [   ");
    // for (int i = 0; i< ptr_node->zones_size; i++)
    // {
    //     printf("%d  ", links_new[i]);
    // }
    // printf("]\n");



    //** >> Creating new children nodes for excess in refinement zones and and linking them to the parent node ptr_node **/
    if (ptr_node->zones_size > ptr_node->chn_size)
    {
        //** >> Space checking in the number of child nodes of ptr_node
        if (space_check(&(ptr_node->chn_cap), ptr_node->zones_size, "p1n2", &(ptr_node->pptr_chn)) == _FAILURE_)
        {
            printf("Error, in space_check function\n");
            return _FAILURE_;
        }

        for (int i = ptr_node->chn_size; i < ptr_node->zones_size; i++)
        {
            ptr_node_ch = (struct node *)malloc(sizeof(struct node));
            initialize_node(ptr_node_ch);
            ptr_node_ch->ID = i;
            ptr_node_ch->lv = ptr_node->lv + 1;
            ptr_node_ch->ptr_pt = ptr_node;
            ptr_node->pptr_chn[i] = ptr_node_ch;
        }
    }

    //** >> Initializing the ptr_aux_idx array to be used as a counter of cells in each current child node**/
    // Notes that the number of refined cells is always bigger or equal than the number of refined zones, so that the capacity of ptr_aux_idx is always enough to counter the number of refined zones
    for (int i = 0; i < ptr_node->chn_size; i++)
    {
        ptr_node->ptr_aux_idx[i] = 0;
    }

    return _SUCCESS_;
}

int tree_adaptation()
{

    //** >> Working in the refinement zones **/
    if (lmin < lmax)
    {
        struct node *ptr_node;
        int no_pts; // Number of parents in the cycle
        int no_lvs; // Number of level of refinement to adapt

        no_lvs = GL_tentacles_level_max < lmax - lmin ? GL_tentacles_level_max : GL_tentacles_level_max - 1;

        for (int lv = no_lvs; lv > -1; lv--)
        {
            no_pts = GL_tentacles_size[lv];
            //** >> For cycle over parent nodes **/
            for (int i = 0; i < no_pts; i++)
            {
                ptr_node = GL_tentacles_old[lv][i];
                //** Updating the box mass information **/
                updating_box_mass(ptr_node);

                //** Initialization of the box_aux **/
                initialization_box_aux(ptr_node);

                //** Initialization of the auiliary refinement arrays**/
                initialization_ref_aux(ptr_node);

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

                //** >> Adapting to new refinement zones **/
                if (adapt_new_ref_zones(ptr_node) == _FAILURE_)
                {
                    printf("Error at function fill_zones_ref()\n");
                    return _FAILURE_;
                }

                // if (node_adaptation(ptr_node) == _FAILURE_)
                // {
                //     printf("Error at function node_adaptation()\n");
                //     return _FAILURE_;
                // }
            }
        }
    }

        return _SUCCESS_;
}

// static int node_adaptation(struct node *ptr_node)
// {

//     return _SUCCESS_;
// }