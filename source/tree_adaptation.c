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

static int updating_cell_struct(struct node *ptr_node)
{

    int no_chn; // Number of child nodes
    struct node *ptr_ch; //child node
    vtype aux_cell_mass; //Total mass in the 8 child cells
    int aux_no_ptcl; // Total mass in the 8 child cells
    int aux_cntr_no_ptcl;

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

            box_idx_x_node = (ptr_ch->ptr_cell_idx_x[j] >> 1) - ptr_node->box_ts_x;
            box_idx_y_node = (ptr_ch->ptr_cell_idx_y[j] >> 1) - ptr_node->box_ts_y;
            box_idx_z_node = (ptr_ch->ptr_cell_idx_z[j] >> 1) - ptr_node->box_ts_z;
            box_idx_node = box_idx_x_node + box_idx_y_node * ptr_node->box_real_dim_x + box_idx_z_node * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;

            aux_cell_mass = 0;
            aux_no_ptcl = 0;

            //** >> Counting mass and number of particles **/
            for (int kk = 0; kk < 2; kk++)
            {
                for (int jj = 0; jj < 2; jj++)
                {
                    for (int ii = 0; ii < 2; ii++)
                    {
                        box_idx_ch = (box_idx_x_ch + ii) + (box_idx_y_ch + jj) * ptr_ch->box_real_dim_x + (box_idx_z_ch + kk) * ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;
                        aux_cell_mass += ptr_ch->ptr_cell_struct[box_idx_ch].cell_mass;
                        aux_no_ptcl += ptr_ch->ptr_cell_struct[box_idx_ch].ptcl_size;
                    }
                }
            }

            ptr_node->ptr_cell_struct[box_idx_node].cell_mass = aux_cell_mass;
            ptr_node->ptr_cell_struct[box_idx_node].ptcl_size = aux_no_ptcl;


            //** >> Space checking of the capacity of the number of particles in the cell **/
            printf("cap = %d, aux_no_ptcl = %d\n",ptr_node->ptr_cell_struct[box_idx_node].ptcl_cap, aux_no_ptcl);
            if (space_check(&(ptr_node->ptr_cell_struct[box_idx_node].ptcl_cap), aux_no_ptcl, 4.0f, "p1i1", &(ptr_node->ptr_cell_struct[box_idx_node].ptr_ptcl)) == _FAILURE_)
            {
                printf("Error, in space_check function\n");
                return _FAILURE_;
            }

            aux_cntr_no_ptcl = 0;

            //** >> Transfering particles **/
            for (int kk = 0; kk < 2; kk++)
            {
                for (int jj = 0; jj < 2; jj++)
                {
                    for (int ii = 0; ii < 2; ii++)
                    {
                        box_idx_ch = (box_idx_x_ch + ii) + (box_idx_y_ch + jj) * ptr_ch->box_real_dim_x + (box_idx_z_ch + kk) * ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;
                        for (int k = 0; k < ptr_ch->ptr_cell_struct[box_idx_ch].ptcl_size; k++)
                        {
                            ptr_node->ptr_cell_struct[box_idx_node].ptr_ptcl[aux_cntr_no_ptcl] = ptr_ch->ptr_cell_struct[box_idx_ch].ptr_ptcl[k];
                            aux_cntr_no_ptcl++;
                        }
                    }
                }
            }

        }
    }

    return _SUCCESS_;
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
    ptr_node->cell_ref_size = 0;    //Number of cell to refine is equal to 0
    for (int i = 0; i < ptr_node->zones_cap; i++)
    {
        ptr_node->ptr_zone_size[i] = 0;    //The zone i contains 0 elements
    }
    ptr_node->zones_size = 0;   //The total number of zones of refinement is 0
}

static int fill_cell_ref(struct node *ptr_node)
{
    //** >> Adding cells which satisfy the refinement criteria or contains grandchild cells to the array ptr_cell_ref and chaning the box_aux the status of refinement -1 **/

    struct node *ptr_ch; // Child node
    struct node *ptr_grandch; // Grandchild node

    int size; // Size of the refinement cells array
    
    int box_idxNbr;   // Box index in the neigborhood

    int cntr; // Counter used to add cell index to the ptr_cell_ref

    int box_idx_x_node;
    int box_idx_y_node;
    int box_idx_z_node;
    int box_idx_node;

    size = 0; // Initial size of the cell refined array must be 0
    //** >> Adding to refine (-1) to cells in the node wich contains grandchildren refinement cells
    for (int ch = 0; ch < ptr_node->chn_size; ch++) // Cycle over children
    {
        ptr_ch = ptr_node->pptr_chn[ch];

        for(int grandch = 0; grandch < ptr_ch->chn_size; grandch++)
        {
            ptr_grandch = ptr_ch->pptr_chn[grandch];

            for (int cell_grandch = 0; cell_grandch < ptr_grandch->cell_size; cell_grandch +=8)
            {
                box_idx_x_node = (ptr_grandch->ptr_cell_idx_x[cell_grandch] >> 2) - ptr_node->box_ts_x;
                box_idx_y_node = (ptr_grandch->ptr_cell_idx_y[cell_grandch] >> 2) - ptr_node->box_ts_y;
                box_idx_z_node = (ptr_grandch->ptr_cell_idx_z[cell_grandch] >> 2) - ptr_node->box_ts_z;
                box_idx_node = box_idx_x_node + box_idx_y_node * ptr_node->box_real_dim_x + box_idx_z_node * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;
                if (ptr_node->ptr_box_aux[box_idx_node] == -3 || ptr_node->ptr_box_aux[box_idx_node] == -2) // Cell has not been added yet
                {
                    size++; // +1 to the total number of cells to be refined
                    //** >> Chaning the cell box status from EXIST (-3) to REFINEMENT REQUIRED (-1) **/
                    ptr_node->ptr_box_aux[box_idx_node] = -1;
                }

                for (int kk = -n_exp; kk < n_exp + 1; kk++)
                {
                    for (int jj = -n_exp; jj < n_exp + 1; jj++)
                    {
                        for (int ii = -n_exp; ii < n_exp + 1; ii++)
                        {
                            box_idxNbr = box_idx_node + ii + jj * ptr_node->box_real_dim_x + kk * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;
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
    }

    //** >> Changing the box status from EXIST (-3) to REFINEMENT REQUIRED (-1) of cells satisfying the refinement criterion **/
    for (int i = 0; i < ptr_node->cell_size; i++)
    {
        box_idx_x_node = ptr_node->ptr_cell_idx_x[i] - ptr_node->box_ts_x;
        box_idx_y_node = ptr_node->ptr_cell_idx_y[i] - ptr_node->box_ts_y;
        box_idx_z_node = ptr_node->ptr_cell_idx_z[i] - ptr_node->box_ts_z;
        box_idx_node = box_idx_x_node + box_idx_y_node * ptr_node->box_real_dim_x + box_idx_z_node * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;

        // Refinement criterion in the box_mass in no border box points
        if ((ptr_node->ptr_box_aux[box_idx_node] == -3 || ptr_node->ptr_box_aux[box_idx_node] == -1) && ptr_node->ptr_cell_struct[box_idx_node].cell_mass >= ref_criterion_mass) // No border (-2)
        {
            if (ptr_node->ptr_box_aux[box_idx_node] == -3) // Cell has not been added yet
            {
                size++;
                //** >> Chaning the cell box status from EXIST (-3) to REFINEMENT REQUIRED (-1) **/
                ptr_node->ptr_box_aux[box_idx_node] = -1;
            }

            //** >> Changing the neighboring cell status **/
            for (int kk = -n_exp; kk < n_exp + 1; kk++)
            {
                for (int jj = -n_exp; jj < n_exp + 1; jj++)
                {
                    for (int ii = -n_exp; ii < n_exp + 1; ii++)
                    {
                        box_idxNbr = box_idx_node + ii + jj * ptr_node->box_real_dim_x + kk * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;
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
    if (space_check(&(ptr_node->cell_ref_cap), size, 2.0f, "p1i1", &(ptr_node->ptr_cell_ref)) == _FAILURE_)
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
        box_idx_x_node = ptr_node->ptr_cell_idx_x[i] - ptr_node->box_ts_x;
        box_idx_y_node = ptr_node->ptr_cell_idx_y[i] - ptr_node->box_ts_y;
        box_idx_z_node = ptr_node->ptr_cell_idx_z[i] - ptr_node->box_ts_z;
        box_idx_node = box_idx_x_node + box_idx_y_node * ptr_node->box_real_dim_x + box_idx_z_node * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;

        if (ptr_node->ptr_box_aux[box_idx_node] == -1) // Cell require refinement
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

    if (space_check(&(ptr_node->aux_idx_cap), ptr_node->cell_ref_cap, 2.0f, "p1i1", &(ptr_node->ptr_aux_idx)) == _FAILURE_)
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
    } // At this point the box_aux contains the information of all new refinement zones

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

int create_links(struct node *ptr_node, int **links_old_ord_old, int **links_new_ord_old, int **links_old_ord_new, int **links_new_ord_new, int *ptr_links_cap)
{
    struct node *ptr_ch; // child node

    int cntr_links; // Counter the number of links between the new and old zones of refinement

    int cntr_links_plus;    // Counte links
    int cntr_link_elements; // Counter linked elements
    bool check_link;        // Check if the element is linked

    int box_value_new; // Value in the new box

    int box_idx_x; // Box index at X direction
    int box_idx_y; // Box index at Y direction
    int box_idx_z; // Box index at Z direction
    int box_idx;   // Box index

    int aux_min;
    int element_idx;
    int aux_int;

    //** >> Creating the links between old and new refinement zones IDs **/

    cntr_links = 0;

    //** >> Space checking of the capacity of links order arrays **/
    if (space_check(ptr_links_cap, ptr_node->zones_size, 4.0f, "p4i1i1i1i1", links_old_ord_old, links_new_ord_old, links_old_ord_new, links_new_ord_new) == _FAILURE_)
    {
        printf("Error, in space_check function\n");
        return _FAILURE_;
    }

    for (int i = 0; i < ptr_node->zones_size; i++)
    {
        (*links_new_ord_old)[i] = -1;
        (*links_old_ord_old)[i] = -1;
    }

    if(ptr_node->chn_size > 0)
    {
        int **pptr_cntr_zones;
        pptr_cntr_zones = (int **)malloc(ptr_node->chn_size * sizeof(int *));
        for (int i = 0; i< ptr_node->chn_size; i++)
        {
            pptr_cntr_zones[i] = (int *)calloc(ptr_node->zones_size, sizeof(int));
        }

        //Filling the repetitions zone number of each child
        for (int ch = 0; ch < ptr_node->chn_size; ch++)
        {
            ptr_ch = ptr_node->pptr_chn[ch];
            for(int cell_idx = 0; cell_idx < ptr_ch->cell_size; cell_idx+=8)
            {
                box_idx_x = (ptr_ch->ptr_cell_idx_x[cell_idx] >> 1) - ptr_node->box_ts_x;
                box_idx_y = (ptr_ch->ptr_cell_idx_y[cell_idx] >> 1) - ptr_node->box_ts_y;
                box_idx_z = (ptr_ch->ptr_cell_idx_z[cell_idx] >> 1) - ptr_node->box_ts_z;
                box_idx = box_idx_x + box_idx_y * ptr_node->box_real_dim_x + box_idx_z * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;
                box_value_new = ptr_node->ptr_box_aux[box_idx];
                if(box_value_new >= 0)
                {
                    pptr_cntr_zones[ch][box_value_new] += 1;
                }
            }
        }

        int *ptr_higher_repetitions, *ptr_zone_higher_repetitions;
        ptr_higher_repetitions = (int *)calloc(ptr_node->chn_size , sizeof(int));
        ptr_zone_higher_repetitions = (int *) calloc(ptr_node->chn_size , sizeof(int));

        //Finding the higher number of zone repetitions of each child
        for(int ch = 0; ch < ptr_node->chn_size; ch++)
        {
            ptr_higher_repetitions[ch] = pptr_cntr_zones[ch][0];

            for (int zone_idx = 1; zone_idx < ptr_node->zones_size; zone_idx++)
            {
                if (ptr_higher_repetitions[ch] < pptr_cntr_zones[ch][zone_idx])
                {
                    ptr_higher_repetitions[ch] = pptr_cntr_zones[ch][zone_idx];
                    ptr_zone_higher_repetitions[ch] = zone_idx;
                }
            }
        }

        //Reorganization of pptr_cntr_zones according to the biggest number of zone repetitions
        int aux_pos;
        int aux_int;
        int *ptr_ch_ID;
        ptr_ch_ID = (int *) malloc(ptr_node->chn_size * sizeof(int));

        for (int ch = 0; ch < ptr_node->chn_size; ch ++)
        {
            ptr_ch_ID[ch] = ch;
        }

        for (int ch = 0; ch < ptr_node->chn_size; ch++)
        {
            aux_pos = ch;
            for (int i = ch + 1; i < ptr_node->chn_size; i++)
            {
                if (ptr_higher_repetitions[ptr_ch_ID[aux_pos]] < ptr_higher_repetitions[ptr_ch_ID[i]])
                {
                    aux_pos = i;
                }
            }
            aux_int = ptr_ch_ID[ch];
            ptr_ch_ID[ch] = ptr_ch_ID[aux_pos];
            ptr_ch_ID[aux_pos] = aux_int;

            for(int i = ch + 1; i < ptr_node->chn_size; i++ )
            {
                pptr_cntr_zones[ptr_ch_ID[i]][ptr_zone_higher_repetitions[ptr_ch_ID[ch]]] = -1;
                if (ptr_zone_higher_repetitions[ptr_ch_ID[i]] == ptr_zone_higher_repetitions[ptr_ch_ID[ch]])
                {
                    ptr_higher_repetitions[ptr_ch_ID[i]] = pptr_cntr_zones[ptr_ch_ID[i]][0];
                    ptr_zone_higher_repetitions[ptr_ch_ID[i]] = 0;

                    for (int zone_idx = 1; zone_idx < ptr_node->zones_size; zone_idx++)
                    {
                        if (ptr_higher_repetitions[ptr_ch_ID[i]] < pptr_cntr_zones[ptr_ch_ID[i]][zone_idx])
                        {
                            ptr_higher_repetitions[ptr_ch_ID[i]] = pptr_cntr_zones[ptr_ch_ID[i]][zone_idx];
                            ptr_zone_higher_repetitions[ptr_ch_ID[i]] = zone_idx;
                        }
                    }
                }
            }


            if (ptr_node->zones_size < ch + 1)
            {
                ch = ptr_node->chn_size;
            }
        }

        while (cntr_links < ptr_node->zones_size && cntr_links < ptr_node->chn_size)
        {
            (*links_old_ord_old)[cntr_links] = ptr_ch_ID[cntr_links];
            (*links_new_ord_old)[cntr_links] = ptr_zone_higher_repetitions[ptr_ch_ID[cntr_links]];
            cntr_links++;
        }

        for (int i = 0; i < ptr_node->chn_size; i++)
        {
            free(pptr_cntr_zones[i]);
        }
        free(pptr_cntr_zones);
        free(ptr_higher_repetitions);
        free(ptr_zone_higher_repetitions);
        free(ptr_ch_ID);
    }

    if (ptr_node->chn_size < ptr_node->zones_size)
    {
        cntr_link_elements = -1;
        for (int ch = ptr_node->chn_size; ch < ptr_node->zones_size; ch++)
        {
            (*links_old_ord_old)[ch] = ch;

            check_link = true;
            while (check_link == true)
            {
                check_link = false;
                cntr_link_elements++;
                for (int j = 0; j < cntr_links; j++)
                {
                    if ((*links_new_ord_old)[j] == cntr_link_elements)
                    {
                        check_link = true;
                        j = cntr_links;
                    }
                }
            }
            (*links_new_ord_old)[ch] = cntr_link_elements;
        }
    }

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
            if (aux_min > (*links_old_ord_old)[cntr_link_elements])
            {
                aux_min = (*links_old_ord_old)[cntr_link_elements];
                element_idx = cntr_link_elements;
            }
            cntr_link_elements++;
        }

        if ((*links_old_ord_old)[i] != aux_min)
        {
            //** >> Old **/
            aux_int = (*links_old_ord_old)[i];
            (*links_old_ord_old)[i] = aux_min;
            (*links_old_ord_old)[element_idx] = aux_int;
            //** >> New **/
            aux_int = (*links_new_ord_old)[i];
            (*links_new_ord_old)[i] = (*links_new_ord_old)[element_idx];
            (*links_new_ord_old)[element_idx] = aux_int;
        }
        cntr_links_plus = aux_min + 1;
    }

    //** >> New order **/
    for (int i = 0; i < ptr_node->zones_size; i++)
    {
        (*links_new_ord_new)[i] = i;
    }
    for (int i = 0; i < ptr_node->zones_size; i++)
    {
        for (int j = 0; j < ptr_node->zones_size; j++)
        {
            if ((*links_new_ord_old)[j] == i)
            {
                (*links_old_ord_new)[i] = (*links_old_ord_old)[j];
                j = ptr_node->zones_size;
            }
        }
    }

        return _SUCCESS_;
}

static void remov_cells_nolonger_require_refinement(struct node *ptr_node, const int *links_old_ord_old, const int *links_new_ord_old)
{

    struct node *ptr_ch;

    int box_idx_x_ch; // Box index in X direcction of the child cell
    int box_idx_y_ch; // Box index in Y direcction of the child cell
    int box_idx_z_ch; // Box index in Z direcction of the child cell
    int box_idx_ch;   // Box index of the child cell

    int box_idx_x_node; // Box index in X direcction of the node cell
    int box_idx_y_node; // Box index in Y direcction of the node cell
    int box_idx_z_node; // Box index in Z direcction of the node cell
    int box_idx_node;   // Box index of the node cell

    int no_cells; // Total number of cells in the node

    // Cycle over new refinement zones
    for (int zone_idx = 0; zone_idx < ptr_node->zones_size; zone_idx++)
    {
        if (links_old_ord_old[zone_idx] < ptr_node->chn_size)
        {
            ptr_ch = ptr_node->pptr_chn[links_old_ord_old[zone_idx]];

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
                    box_idx_x_ch = ptr_ch->ptr_cell_idx_x[cell_idx] - ptr_ch->box_ts_x;
                    box_idx_y_ch = ptr_ch->ptr_cell_idx_y[cell_idx] - ptr_ch->box_ts_y;
                    box_idx_z_ch = ptr_ch->ptr_cell_idx_z[cell_idx] - ptr_ch->box_ts_z;

                    //** >> Removing particles in the cell from the child node **/
                    // Updating local mass, cell struct and box
                    for (int kk = 0; kk < 2; kk++)
                    {
                        for (int jj = 0; jj < 2; jj++)
                        {
                            for (int ii = 0; ii < 2; ii++)
                            {
                                box_idx_ch = (box_idx_x_ch + ii) + (box_idx_y_ch + jj) * ptr_ch->box_real_dim_x + (box_idx_z_ch + kk) * ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;
                                ptr_ch->local_mass -= ptr_ch->ptr_cell_struct[box_idx_ch].cell_mass;
                                ptr_ch->ptr_cell_struct[box_idx_ch].cell_mass = 0;
                                ptr_ch->ptr_cell_struct[box_idx_ch].ptcl_size = 0;
                                ptr_ch->ptr_box[box_idx_ch] = -4;
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
            ptr_ch->cell_size = no_cells;
        }
    }
}

static int adapt_child_box_and_cells(struct node *ptr_node, const int *links_old_ord_old, const int *links_new_ord_old)
{
    struct node *ptr_ch;

    int box_idx_x_ch_old; // Box index in X direcction of the child cell
    int box_idx_y_ch_old; // Box index in Y direcction of the child cell
    int box_idx_z_ch_old; // Box index in Z direcction of the child cell
    int box_idx_ch_old;   // Box index of the child cell

    int box_idx_x_ch_new; // Box index in X direcction of the child cell
    int box_idx_y_ch_new; // Box index in Y direcction of the child cell
    int box_idx_z_ch_new; // Box index in Z direcction of the child cell
    int box_idx_ch_new;   // Box index of the child cell

    int new_box_min_x;
    int new_box_min_y;
    int new_box_min_z;
    int new_box_max_x;
    int new_box_max_y;
    int new_box_max_z;

    int old_box_real_dim_x_ch;
    int old_box_real_dim_y_ch;
    //int old_box_real_dim_z_ch;

    int old_box_ts_x_ch;
    int old_box_ts_y_ch;
    int old_box_ts_z_ch;

    int cell_idx; // The cell index is simply i of the for loop

    int new_zone_idx; // ID of the new zone of refinement

    int size;

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
            size = 8 * ptr_node->ptr_zone_size[new_zone_idx] + ptr_ch->cell_size;

            //** >> Space checking of cells indexes of the child node
            if (space_check(&(ptr_ch->cell_cap), size, 2.0f, "p3i1i1i1", &(ptr_ch->ptr_cell_idx_x), &(ptr_ch->ptr_cell_idx_y), &(ptr_ch->ptr_cell_idx_z)) == _FAILURE_)
            {
                printf("Error, in space_check function\n");
                return _FAILURE_;
            }

            //** BOXES **/
            ptr_ch->box_check_fit = true;

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

            // Size of the "minimal new box"
            ptr_ch->box_dim_x = new_box_max_x - new_box_min_x + 1;
            ptr_ch->box_dim_y = new_box_max_y - new_box_min_y + 1;
            ptr_ch->box_dim_z = new_box_max_z - new_box_min_z + 1;

            //** >> Checking if the new dimensions fit in the old box **/
            aux_int = new_box_min_x - ptr_ch->box_ts_x;
            aux_int = aux_int < (new_box_min_y - ptr_ch->box_ts_y) ? aux_int : (new_box_min_y - ptr_ch->box_ts_y);
            aux_int = aux_int < (new_box_min_z - ptr_ch->box_ts_z) ? aux_int : (new_box_min_z - ptr_ch->box_ts_z);
            if (aux_int < bder_box)
            {
                ptr_ch->box_check_fit = false;
            }
            aux_int = new_box_max_x - ptr_ch->box_ts_x;
            if (aux_int > ptr_ch->box_real_dim_x - bder_box)
            {
                ptr_ch->box_check_fit = false;
            }
            aux_int = new_box_max_y - ptr_ch->box_ts_y;
            if (aux_int > ptr_ch->box_real_dim_y - bder_box)
            {
                ptr_ch->box_check_fit = false;
            }
            aux_int = new_box_max_z - ptr_ch->box_ts_z;
            if (aux_int > ptr_ch->box_real_dim_z - bder_box)
            {
                ptr_ch->box_check_fit = false;
            }

            //** >> The new box does not fit in the old box **/
            if (ptr_ch->box_check_fit == false)
            {
                // Real dimensions of the old box
                old_box_real_dim_x_ch = ptr_ch->box_real_dim_x;
                old_box_real_dim_y_ch = ptr_ch->box_real_dim_y;
                //old_box_real_dim_z_ch = ptr_ch->box_real_dim_z;

                // Translations between cell array and old box (Used to access only to the mass box)
                old_box_ts_x_ch = ptr_ch->box_ts_x; // Every cell in the level l in the box must be subtracted this value to obtain the box index
                old_box_ts_y_ch = ptr_ch->box_ts_y;
                old_box_ts_z_ch = ptr_ch->box_ts_z;

                // Real dimensions of the new box
                ptr_ch->box_real_dim_x = 5 > (n_exp - 1) ? (ptr_ch->box_dim_x + 10) : (ptr_ch->box_dim_x + 2 * n_exp - 2);
                ptr_ch->box_real_dim_y = 5 > (n_exp - 1) ? (ptr_ch->box_dim_y + 10) : (ptr_ch->box_dim_y + 2 * n_exp - 2);
                ptr_ch->box_real_dim_z = 5 > (n_exp - 1) ? (ptr_ch->box_dim_z + 10) : (ptr_ch->box_dim_z + 2 * n_exp - 2);

                // Translations between cell array and new box
                pos_x = (ptr_ch->box_real_dim_x - ptr_ch->box_dim_x) / 2; // Half of the distance of the box side less the "minimal box" side
                ptr_ch->box_ts_x = new_box_min_x - pos_x;             // Every cell in the level l in the box must be subtracted this value to obtain the box index
                pos_y = (ptr_ch->box_real_dim_y - ptr_ch->box_dim_y) / 2;
                ptr_ch->box_ts_y = new_box_min_y - pos_y;
                pos_z = (ptr_ch->box_real_dim_z - ptr_ch->box_dim_z) / 2;
                ptr_ch->box_ts_z = new_box_min_z - pos_z;

                //** >> Transfer from cell structure to auxiliary cell structure **/
                for (int cell_idx = 0; cell_idx < ptr_ch->cell_size; cell_idx++)
                {
                    box_idx_x_ch_old = ptr_ch->ptr_cell_idx_x[cell_idx] - old_box_ts_x_ch;
                    box_idx_y_ch_old = ptr_ch->ptr_cell_idx_y[cell_idx] - old_box_ts_y_ch;
                    box_idx_z_ch_old = ptr_ch->ptr_cell_idx_z[cell_idx] - old_box_ts_z_ch;
                    box_idx_ch_old = box_idx_x_ch_old + box_idx_y_ch_old * old_box_real_dim_x_ch + box_idx_z_ch_old * old_box_real_dim_x_ch * old_box_real_dim_y_ch;

                    //** >> Transfer from cell structure to auxiliary cell structure
                    ptr_ch->ptr_cell_struct_aux[box_idx_ch_old].cell_mass = ptr_ch->ptr_cell_struct[box_idx_ch_old].cell_mass;
                    ptr_ch->ptr_cell_struct_aux[box_idx_ch_old].ptcl_size = ptr_ch->ptr_cell_struct[box_idx_ch_old].ptcl_size;
                    
                    //** >> Space checking of cell structure particles
                    if (space_check(&(ptr_ch->ptr_cell_struct_aux[box_idx_ch_old].ptcl_cap), ptr_ch->ptr_cell_struct_aux[box_idx_ch_old].ptcl_size, 2.0f, "p1i1", &(ptr_ch->ptr_cell_struct_aux[box_idx_ch_old].ptr_ptcl)) == _FAILURE_)
                    {
                        printf("Error, in space_check function\n");
                        return _FAILURE_;
                    }

                    for (int j = 0; j < ptr_ch->ptr_cell_struct_aux[box_idx_ch_old].ptcl_size; j++)
                    {
                        ptr_ch->ptr_cell_struct_aux[box_idx_ch_old].ptr_ptcl[j] = ptr_ch->ptr_cell_struct[box_idx_ch_old].ptr_ptcl[j];
                    }
                }

                size = ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y * ptr_ch->box_real_dim_z;

                //** >> Space checking of boxes and cell structures
                if (space_check(&(ptr_ch->box_cap), size, 1.0f, "p4i1i1c1c1", &(ptr_ch->ptr_box), &(ptr_ch->ptr_box_aux), &(ptr_ch->ptr_cell_struct), &(ptr_ch->ptr_cell_struct_aux)) == _FAILURE_)
                {
                    printf("Error, in space_check function\n");
                    return _FAILURE_;
                }

                //** Transfer old information to new cell structures **/
                for (int j = 0; j < ptr_ch->box_cap; j++)
                {
                    ptr_ch->ptr_cell_struct[j].ptcl_size = 0;
                    ptr_ch->ptr_cell_struct[j].cell_mass = 0;
                }
                for (int cell_idx = 0; cell_idx < ptr_ch->cell_size; cell_idx++)
                {
                    box_idx_x_ch_old = ptr_ch->ptr_cell_idx_x[cell_idx] - old_box_ts_x_ch;
                    box_idx_y_ch_old = ptr_ch->ptr_cell_idx_y[cell_idx] - old_box_ts_y_ch;
                    box_idx_z_ch_old = ptr_ch->ptr_cell_idx_z[cell_idx] - old_box_ts_z_ch;
                    box_idx_ch_old = box_idx_x_ch_old + box_idx_y_ch_old * old_box_real_dim_x_ch + box_idx_z_ch_old * old_box_real_dim_x_ch * old_box_real_dim_y_ch;

                    if (ptr_ch->ptr_cell_struct_aux[box_idx_ch_old].ptcl_size > 0)
                    {
                        box_idx_x_ch_new = ptr_ch->ptr_cell_idx_x[cell_idx] - ptr_ch->box_ts_x;
                        box_idx_y_ch_new = ptr_ch->ptr_cell_idx_y[cell_idx] - ptr_ch->box_ts_y;
                        box_idx_z_ch_new = ptr_ch->ptr_cell_idx_z[cell_idx] - ptr_ch->box_ts_z;
                        box_idx_ch_new = box_idx_x_ch_new + box_idx_y_ch_new * ptr_ch->box_real_dim_x + box_idx_z_ch_new * ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;

                        //** >> Transfer from cell structure to auxiliary cell structure
                        ptr_ch->ptr_cell_struct[box_idx_ch_new].cell_mass = ptr_ch->ptr_cell_struct_aux[box_idx_ch_old].cell_mass;
                        ptr_ch->ptr_cell_struct[box_idx_ch_new].ptcl_size = ptr_ch->ptr_cell_struct_aux[box_idx_ch_old].ptcl_size;

                        //** >> Space checking of cell structure particles
                        if (space_check(&(ptr_ch->ptr_cell_struct[box_idx_ch_new].ptcl_cap), ptr_ch->ptr_cell_struct[box_idx_ch_new].ptcl_size, 2.0f, "p1i1", &(ptr_ch->ptr_cell_struct[box_idx_ch_new].ptr_ptcl)) == _FAILURE_)
                        {
                            printf("Error, in space_check function\n");
                            return _FAILURE_;
                        }

                        for (int j = 0; j < ptr_ch->ptr_cell_struct[box_idx_ch_new].ptcl_size; j++)
                        {
                            ptr_ch->ptr_cell_struct[box_idx_ch_new].ptr_ptcl[j] = ptr_ch->ptr_cell_struct_aux[box_idx_ch_old].ptr_ptcl[j];
                        }
                    }
                }

                // Initializing box at NO-EXIST (-4)
                for (int j = 0; j < size; j++)
                {
                    ptr_ch->ptr_box[j] = -4;
                }

                //* Potential, Acceleration and density capacity of the grid **/
                size = (ptr_ch->box_real_dim_x + 1) * (ptr_ch->box_real_dim_y + 1) * (ptr_ch->box_real_dim_z + 1);
                //** >> Space checking
                if (space_check(&(ptr_ch->grid_properties_cap), size, 1.0f, "p5v1v1v1v1v1", &(ptr_ch->ptr_pot), &(ptr_ch->ptr_ax), &(ptr_ch->ptr_ay), &(ptr_ch->ptr_az), &(ptr_ch->ptr_d)) == _FAILURE_)
                {
                    printf("Error, in space_check function\n");
                    return _FAILURE_;
                }
            }

            //** >> Grid points **/
            ptr_ch->grid_intr_size = 0;
            ptr_ch->grid_bder_size = 0;

            //** >> Reset grid properties: Density  **/
            size = (ptr_ch->box_real_dim_x + 1) * (ptr_ch->box_real_dim_y + 1) * (ptr_ch->box_real_dim_z + 1);
            for (int j = 0; j < size; j++)
            {
                //ptr_ch->ptr_pot[j] = 0;
                //ptr_ch->ptr_ax[j] = 0;
                //ptr_ch->ptr_ay[j] = 0;
                //ptr_ch->ptr_az[j] = 0;
                ptr_ch->ptr_d[j] = 0;
            }
        }
    }
    return _SUCCESS_;
}

static int create_new_child_nodes(struct node *ptr_node, const int *links_old_ord_old, const int *links_new_ord_old)
{

    struct node *ptr_ch; // child node
    int cell_idx;        // The cell index is simply i of the for loop

    int pos_x;  // Distance between the real box and the minimal box when the last is localized in the middle of the real one
    int pos_y;
    int pos_z;

    int size; // Size or number of elements in some array in child nodes

    //** >> Space checking in the number of child nodes of ptr_node
    if (space_check(&(ptr_node->chn_cap), ptr_node->zones_size, 2.0f, "p1n2", &(ptr_node->pptr_chn)) == _FAILURE_)
    {
        printf("Error, in space_check function\n");
        return _FAILURE_;
    }

    for (int i = ptr_node->chn_size; i < ptr_node->zones_size; i++)
    {
        ptr_ch = new_node();

        //** >> Global node properties **/
        ptr_ch->ID = i;
        ptr_ch->lv = ptr_node->lv + 1;
        //** >> Cells in the node **/
        size = 8 * ptr_node->ptr_zone_size[links_new_ord_old[i]];
        ptr_ch->cell_size = 0;

        //** >> Space checking of cells indexes of the child node
        if (space_check(&(ptr_ch->cell_cap), size, 2.0f, "p3i1i1i1", &(ptr_ch->ptr_cell_idx_x), &(ptr_ch->ptr_cell_idx_y), &(ptr_ch->ptr_cell_idx_z)) == _FAILURE_)
        {
            printf("Error, in space_check function\n");
            return _FAILURE_;
        }

        //** >> Boxes **/
        // MIN and MAX cell indexes values of the node.
        ptr_ch->box_min_x = INT_MAX;
        ptr_ch->box_min_y = INT_MAX;
        ptr_ch->box_min_z = INT_MAX;
        ptr_ch->box_max_x = 0;
        ptr_ch->box_max_y = 0;
        ptr_ch->box_max_z = 0;

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
        ptr_ch->box_real_dim_x = (ptr_ch->box_dim_x + 10) > (ptr_ch->box_dim_x + 2 * n_exp - 2) ? (ptr_ch->box_dim_x + 10) : (ptr_ch->box_dim_x + 2 * n_exp - 2);
        ptr_ch->box_real_dim_y = (ptr_ch->box_dim_y + 10) > (ptr_ch->box_dim_y + 2 * n_exp - 2) ? (ptr_ch->box_dim_y + 10) : (ptr_ch->box_dim_y + 2 * n_exp - 2);
        ptr_ch->box_real_dim_z = (ptr_ch->box_dim_z + 10) > (ptr_ch->box_dim_z + 2 * n_exp - 2) ? (ptr_ch->box_dim_z + 10) : (ptr_ch->box_dim_z + 2 * n_exp - 2);

        // Translations between cell array and box
        pos_x = (ptr_ch->box_real_dim_x - ptr_ch->box_dim_x) / 2; // Half of the distance of the box side less the "minimal box" side
        ptr_ch->box_ts_x = ptr_ch->box_min_x - pos_x;             // Every cell in the level l in the box must be subtracted this value to obtain the box index
        pos_y = (ptr_ch->box_real_dim_y - ptr_ch->box_dim_y) / 2;
        ptr_ch->box_ts_y = ptr_ch->box_min_y - pos_y;
        pos_z = (ptr_ch->box_real_dim_z - ptr_ch->box_dim_z) / 2;
        ptr_ch->box_ts_z = ptr_ch->box_min_z - pos_z;

        ptr_ch->box_check_fit = true;

        // Filling the boxes
        size = ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y * ptr_ch->box_real_dim_z;
        //** >> Space checking of boxes: box, box_aux and box_mass_aux
        if (space_check(&(ptr_ch->box_cap), size, 1.0f, "p4i1i1c1c1", &(ptr_ch->ptr_box), &(ptr_ch->ptr_box_aux), &(ptr_ch->ptr_cell_struct), &(ptr_ch->ptr_cell_struct_aux)) == _FAILURE_)
        {
            printf("Error, in space_check function\n");
            return _FAILURE_;
        }

        // Putting the value of NO-EXIST (-4) in every box index
        for (int j = 0; j < size; j++)
        {
            ptr_ch->ptr_box[j] = -4;
        }

        //** >> Particles in the node **/
        //** >> Cell structure **/
        for(int j = 0; j< ptr_ch->box_cap; j++)
        {
            ptr_ch->ptr_cell_struct[j].cell_mass = 0;
            ptr_ch->ptr_cell_struct[j].ptcl_size = 0;
        }
        ptr_ch->local_mass = 0;

        //** >> Grid points **/
        ptr_ch->grid_intr_size = 0;
        ptr_ch->grid_bder_size = 0;

        //** >> Potential, acceleration and density of the grid **/
        size = (ptr_ch->box_real_dim_x + 1) * (ptr_ch->box_real_dim_y + 1) * (ptr_ch->box_real_dim_z + 1);
        //** >> Space checking
        if (space_check(&(ptr_ch->grid_properties_cap), size, 1.0f, "p5v1v1v1v1v1", &(ptr_ch->ptr_pot), &(ptr_ch->ptr_ax), &(ptr_ch->ptr_ay), &(ptr_ch->ptr_az), &(ptr_ch->ptr_d)) == _FAILURE_)
        {
            printf("Error, in space_check function\n");
            return _FAILURE_;
        }

        //** >> Reset grid properties **/
        for (int j = 0; j < size; j++)
        {
            //ptr_ch->ptr_pot[j] = 0;
            //ptr_ch->ptr_ax[j] = 0;
            //ptr_ch->ptr_ay[j] = 0;
            //ptr_ch->ptr_az[j] = 0;
            ptr_ch->ptr_d[j] = 0;
        }

        //** >> Tree structure **/
        ptr_ch->chn_size = 0;
        ptr_ch->ptr_pt = ptr_node;
        ptr_node->pptr_chn[i] = ptr_ch;

        //** >> Auxililary arrays to go from old box to new box **/
        ptr_ch->cell_ref_size = 0;
        for (int j = 0; j < ptr_ch->zones_size; j++)
        {
            ptr_ch->ptr_zone_size[j] = 0;
        }
        ptr_ch->zones_size = 0;
    }
    return _SUCCESS_;
}

static int moving_old_child_to_new_child(struct node *ptr_node, const int *links_old_ord_old, const int *links_new_ord_old, const int *links_old_ord_new, const int *links_new_ord_new)
{
    struct node *ptr_ch_A;
    struct node *ptr_ch_B;
    int no_cells_ch_A; // Total number of cells in the child node A
    int no_cells_ch_B; // Total number of cells in the child node A

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

    int box_idxNbr_ch_A;
    int box_idxNbr_ch_B;

    int aux_const;

    // Cycle over new refinement zones,
    // Case of old child nodes to be reused
    printf("A\n");
    for (int zone_idx = 0; zone_idx < ptr_node->zones_size; zone_idx++)
    {
        if (links_old_ord_old[zone_idx] < ptr_node->chn_size)
        {
            new_zone_idx = links_new_ord_old[zone_idx];

            ptr_ch_A = ptr_node->pptr_chn[links_old_ord_old[zone_idx]];
            no_cells_ch_A = ptr_ch_A->cell_size;

           
            if(ptr_ch_A->box_check_fit == true)
            {
                 printf("A1\n");
                for (int cell_idx = 0; cell_idx < no_cells_ch_A; cell_idx += 8)
                {
                    box_idx_x_node = (ptr_ch_A->ptr_cell_idx_x[cell_idx] >> 1) - ptr_node->box_ts_x;
                    box_idx_y_node = (ptr_ch_A->ptr_cell_idx_y[cell_idx] >> 1) - ptr_node->box_ts_y;
                    box_idx_z_node = (ptr_ch_A->ptr_cell_idx_z[cell_idx] >> 1) - ptr_node->box_ts_z;
                    box_idx_node = box_idx_x_node + box_idx_y_node * ptr_node->box_real_dim_x + box_idx_z_node * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;

                    if (ptr_node->ptr_box_aux[box_idx_node] != new_zone_idx)
                    {
printf("A11\n");
                        box_idx_x_ch_A = ptr_ch_A->ptr_cell_idx_x[cell_idx] - ptr_ch_A->box_ts_x;
                        box_idx_y_ch_A = ptr_ch_A->ptr_cell_idx_y[cell_idx] - ptr_ch_A->box_ts_y;
                        box_idx_z_ch_A = ptr_ch_A->ptr_cell_idx_z[cell_idx] - ptr_ch_A->box_ts_z;
                        box_idx_ch_A = box_idx_x_ch_A + box_idx_y_ch_A * ptr_ch_A->box_real_dim_x + box_idx_z_ch_A * ptr_ch_A->box_real_dim_x * ptr_ch_A->box_real_dim_y;

                        //  ptr_ch_B: Is the child node where the information from the child A need to be sent
                        ptr_ch_B = ptr_node->pptr_chn[links_old_ord_new[ptr_node->ptr_box_aux[box_idx_node]]];

                        //** >> CELLS **/
                        no_cells_ch_B = ptr_ch_B->cell_size;

                        box_idx_x_ch_B = ptr_ch_A->ptr_cell_idx_x[cell_idx] - ptr_ch_B->box_ts_x;
                        box_idx_y_ch_B = ptr_ch_A->ptr_cell_idx_y[cell_idx] - ptr_ch_B->box_ts_y;
                        box_idx_z_ch_B = ptr_ch_A->ptr_cell_idx_z[cell_idx] - ptr_ch_B->box_ts_z;
                        box_idx_ch_B = box_idx_x_ch_B + box_idx_y_ch_B * ptr_ch_B->box_real_dim_x + box_idx_z_ch_B * ptr_ch_B->box_real_dim_x * ptr_ch_B->box_real_dim_y;

                        printf("A12\n");
                        for (int kk = 0; kk < 2; kk++)
                        {
                            for (int jj = 0; jj < 2; jj++)
                            {
                                for (int ii = 0; ii < 2; ii++)
                                {
                                    aux_const = ii + jj * 2 + kk * 4;
                                    // Moving from child A to child B
                                    ptr_ch_B->ptr_cell_idx_x[no_cells_ch_B + aux_const] = ptr_ch_A->ptr_cell_idx_x[cell_idx + aux_const];
                                    ptr_ch_B->ptr_cell_idx_y[no_cells_ch_B + aux_const] = ptr_ch_A->ptr_cell_idx_y[cell_idx + aux_const];
                                    ptr_ch_B->ptr_cell_idx_z[no_cells_ch_B + aux_const] = ptr_ch_A->ptr_cell_idx_z[cell_idx + aux_const];

                                    // Moving from end of child A to removed cells of child A
                                    ptr_ch_A->ptr_cell_idx_x[cell_idx + aux_const] = ptr_ch_A->ptr_cell_idx_x[no_cells_ch_A - 8 + aux_const];
                                    ptr_ch_A->ptr_cell_idx_y[cell_idx + aux_const] = ptr_ch_A->ptr_cell_idx_y[no_cells_ch_A - 8 + aux_const];
                                    ptr_ch_A->ptr_cell_idx_z[cell_idx + aux_const] = ptr_ch_A->ptr_cell_idx_z[no_cells_ch_A - 8 + aux_const];

                                    // Updating the box status of child node A and B
                                    box_idxNbr_ch_A = box_idx_ch_A + ii + jj * ptr_ch_A->box_real_dim_x + kk * ptr_ch_A->box_real_dim_x * ptr_ch_A->box_real_dim_y;
                                    box_idxNbr_ch_B = box_idx_ch_B + ii + jj * ptr_ch_B->box_real_dim_x + kk * ptr_ch_B->box_real_dim_x * ptr_ch_B->box_real_dim_y;
                                    ptr_ch_A->ptr_box[box_idxNbr_ch_A] = -4; // Putting the status of NO EXIST (-4) in the child node cell
                                    ptr_ch_B->ptr_box[box_idxNbr_ch_B] = -3; // Putting the status of EXIST (-3) in the child node cell
                                }
                            }
                        }

                        printf("A13\n");
                        //** >> Cell structure and local mass **/

                        for (int kk = 0; kk < 2; kk++)
                        {
                            for (int jj = 0; jj < 2; jj++)
                            {
                                for (int ii = 0; ii < 2; ii++)
                                {
                                    printf("D\n");
                                    box_idxNbr_ch_A = box_idx_ch_A + ii + jj * ptr_ch_A->box_real_dim_x + kk * ptr_ch_A->box_real_dim_x * ptr_ch_A->box_real_dim_y;
                                    box_idxNbr_ch_B = box_idx_ch_B + ii + jj * ptr_ch_B->box_real_dim_x + kk * ptr_ch_B->box_real_dim_x * ptr_ch_B->box_real_dim_y;
                                            printf("E\n"); 
                                            printf("box_idxNbr_ch_A = %d, box_idxNbr_ch_B = %d\n",box_idxNbr_ch_A,box_idxNbr_ch_B); 
                                            printf("ptr_ch_A:\n");
                                            printf("Box_idx_A = %d\n",box_idx_ch_A);
                                            printf("box_real_dim_x = %d\n", ptr_ch_A->box_real_dim_x);
                                            printf("box_real_dim_y = %d\n", ptr_ch_A->box_real_dim_y);
                                            printf("box_real_dim_z = %d\n", ptr_ch_A->box_real_dim_z);
                                            printf("\n");
                                            printf("ptr_ch_B:\n");
                                            printf("Box_idx_B = %d\n",box_idx_ch_B);
                                            printf("box_real_dim_x = %d\n", ptr_ch_B->box_real_dim_x);
                                            printf("box_real_dim_y = %d\n", ptr_ch_B->box_real_dim_y);
                                            printf("box_real_dim_z = %d\n", ptr_ch_B->box_real_dim_z);
                                    ptr_ch_B->ptr_cell_struct[box_idxNbr_ch_B].cell_mass = ptr_ch_A->ptr_cell_struct[box_idxNbr_ch_A].cell_mass;
                                    ptr_ch_B->ptr_cell_struct[box_idxNbr_ch_B].ptcl_size = ptr_ch_A->ptr_cell_struct[box_idxNbr_ch_A].ptcl_size;
                                            printf("F\n");
                                    //** >> Space checking of cell structure particles
                                    if (space_check(&(ptr_ch_B->ptr_cell_struct[box_idxNbr_ch_B].ptcl_cap), ptr_ch_B->ptr_cell_struct[box_idxNbr_ch_B].ptcl_size , 2.0f, "p1i1", &(ptr_ch_B->ptr_cell_struct[box_idxNbr_ch_B].ptr_ptcl)) == _FAILURE_)
                                    {
                                        printf("Error, in space_check function\n");
                                        return _FAILURE_;
                                    }
                                    printf("G\n");
                                    for (int j = 0; j < ptr_ch_B->ptr_cell_struct[box_idxNbr_ch_B].ptcl_size; j++)
                                    {
                                        ptr_ch_B->ptr_cell_struct[box_idxNbr_ch_B].ptr_ptcl[j] = ptr_ch_A->ptr_cell_struct[box_idxNbr_ch_A].ptr_ptcl[j];
                                    }
                                    printf("H\n");
                                    ptr_ch_B->local_mass += ptr_ch_A->ptr_cell_struct[box_idxNbr_ch_A].cell_mass;
                                    ptr_ch_A->local_mass -= ptr_ch_A->ptr_cell_struct[box_idxNbr_ch_A].cell_mass;
                                    ptr_ch_A->ptr_cell_struct[box_idxNbr_ch_A].cell_mass = 0;
                                    ptr_ch_A->ptr_cell_struct[box_idxNbr_ch_A].ptcl_size = 0;
                                    printf("I\n");
                                }
                            }
                        }
                        printf("A14\n");
                        ptr_ch_B->cell_size += 8;
                        no_cells_ch_A -= 8;
                        cell_idx -= 8;
                    }

                }
            }
            
            else
            {
                printf("B\n");
                for (int cell_idx = 0; cell_idx < no_cells_ch_A; cell_idx += 8)
                {
                    box_idx_x_node = (ptr_ch_A->ptr_cell_idx_x[cell_idx] >> 1) - ptr_node->box_ts_x;
                    box_idx_y_node = (ptr_ch_A->ptr_cell_idx_y[cell_idx] >> 1) - ptr_node->box_ts_y;
                    box_idx_z_node = (ptr_ch_A->ptr_cell_idx_z[cell_idx] >> 1) - ptr_node->box_ts_z;
                    box_idx_node = box_idx_x_node + box_idx_y_node * ptr_node->box_real_dim_x + box_idx_z_node * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;
                    
                    box_idx_x_ch_A = ptr_ch_A->ptr_cell_idx_x[cell_idx] - ptr_ch_A->box_ts_x;
                    box_idx_y_ch_A = ptr_ch_A->ptr_cell_idx_y[cell_idx] - ptr_ch_A->box_ts_y;
                    box_idx_z_ch_A = ptr_ch_A->ptr_cell_idx_z[cell_idx] - ptr_ch_A->box_ts_z;
                    box_idx_ch_A = box_idx_x_ch_A + box_idx_y_ch_A * ptr_ch_A->box_real_dim_x + box_idx_z_ch_A * ptr_ch_A->box_real_dim_x * ptr_ch_A->box_real_dim_y;

                    if (ptr_node->ptr_box_aux[box_idx_node] == new_zone_idx)
                    {
                        //** changing the box and transfer box_mass_aux to box_mass of the cells 
                        if (ptr_ch_A->ptr_box[box_idx_ch_A] == -4)
                        {
                            for (int kk = 0; kk < 2; kk++)
                            {
                                for (int jj = 0; jj < 2; jj++)
                                {
                                    for (int ii = 0; ii < 2; ii++)
                                    {
                                        box_idxNbr_ch_A = box_idx_ch_A + ii + jj * ptr_ch_A->box_real_dim_x + kk * ptr_ch_A->box_real_dim_x * ptr_ch_A->box_real_dim_y;
                                        ptr_ch_A->ptr_box[box_idxNbr_ch_A] = -3; // Putting the status of EXIST (-3) in the child node cell
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        ptr_ch_B = ptr_node->pptr_chn[links_old_ord_new[ptr_node->ptr_box_aux[box_idx_node]]];
                        
                        //** >> CELLS **/
                        no_cells_ch_B = ptr_ch_B->cell_size;

                        box_idx_x_ch_B = ptr_ch_A->ptr_cell_idx_x[cell_idx] - ptr_ch_B->box_ts_x;
                        box_idx_y_ch_B = ptr_ch_A->ptr_cell_idx_y[cell_idx] - ptr_ch_B->box_ts_y;
                        box_idx_z_ch_B = ptr_ch_A->ptr_cell_idx_z[cell_idx] - ptr_ch_B->box_ts_z;
                        box_idx_ch_B = box_idx_x_ch_B + box_idx_y_ch_B * ptr_ch_B->box_real_dim_x + box_idx_z_ch_B * ptr_ch_B->box_real_dim_x * ptr_ch_B->box_real_dim_y;

                        for (int kk = 0; kk < 2; kk++)
                        {
                            for (int jj = 0; jj < 2; jj++)
                            {
                                for (int ii = 0; ii < 2; ii++)
                                {
                                    aux_const = ii + jj * 2 + kk * 4;

                                    // Moving from child A to child B
                                    ptr_ch_B->ptr_cell_idx_x[no_cells_ch_B + aux_const] = ptr_ch_A->ptr_cell_idx_x[cell_idx + aux_const];
                                    ptr_ch_B->ptr_cell_idx_y[no_cells_ch_B + aux_const] = ptr_ch_A->ptr_cell_idx_y[cell_idx + aux_const];
                                    ptr_ch_B->ptr_cell_idx_z[no_cells_ch_B + aux_const] = ptr_ch_A->ptr_cell_idx_z[cell_idx + aux_const];
                                    // Moving from end of child A to removed cells of child A
                                    ptr_ch_A->ptr_cell_idx_x[cell_idx + aux_const] = ptr_ch_A->ptr_cell_idx_x[no_cells_ch_A - 8 + aux_const];
                                    ptr_ch_A->ptr_cell_idx_y[cell_idx + aux_const] = ptr_ch_A->ptr_cell_idx_y[no_cells_ch_A - 8 + aux_const];
                                    ptr_ch_A->ptr_cell_idx_z[cell_idx + aux_const] = ptr_ch_A->ptr_cell_idx_z[no_cells_ch_A - 8 + aux_const];

                                    // Updating the box status of child node B
                                    box_idxNbr_ch_B = box_idx_ch_B + ii + jj * ptr_ch_B->box_real_dim_x + kk * ptr_ch_B->box_real_dim_x * ptr_ch_B->box_real_dim_y;
                                    ptr_ch_B->ptr_box[box_idx_ch_B] = -3; // Putting the status of EXIST (-3) in the child node cell
                                }
                            }
                        }
                        ptr_ch_B->cell_size += 8;
                        no_cells_ch_A -= 8;
                        cell_idx -= 8;

                        //** >> Cell structure and local mass **/
                        for (int kk = 0; kk < 2; kk++)
                        {
                            for (int jj = 0; jj < 2; jj++)
                            {
                                for (int ii = 0; ii < 2; ii++)
                                {
                                    box_idxNbr_ch_A = box_idx_ch_A + ii + jj * ptr_ch_A->box_real_dim_x + kk * ptr_ch_A->box_real_dim_x * ptr_ch_A->box_real_dim_y;
                                    box_idxNbr_ch_B = box_idx_ch_B + ii + jj * ptr_ch_B->box_real_dim_x + kk * ptr_ch_B->box_real_dim_x * ptr_ch_B->box_real_dim_y;

                                    ptr_ch_B->ptr_cell_struct[box_idxNbr_ch_B].cell_mass = ptr_ch_A->ptr_cell_struct[box_idxNbr_ch_A].cell_mass;
                                    ptr_ch_B->ptr_cell_struct[box_idxNbr_ch_B].ptcl_size = ptr_ch_A->ptr_cell_struct[box_idxNbr_ch_A].ptcl_size;

                                    //** >> Space checking of cell structure particles
                                    if (space_check(&(ptr_ch_B->ptr_cell_struct[box_idxNbr_ch_B].ptcl_cap), ptr_ch_B->ptr_cell_struct[box_idxNbr_ch_B].ptcl_size, 2.0f, "p1i1", &(ptr_ch_B->ptr_cell_struct[box_idxNbr_ch_B].ptr_ptcl)) == _FAILURE_)
                                    {
                                        printf("Error, in space_check function\n");
                                        return _FAILURE_;
                                    }

                                    for (int j = 0; j < ptr_ch_B->ptr_cell_struct[box_idxNbr_ch_B].ptcl_size; j++)
                                    {
                                        ptr_ch_B->ptr_cell_struct[box_idxNbr_ch_B].ptr_ptcl[j] = ptr_ch_A->ptr_cell_struct[box_idxNbr_ch_A].ptr_ptcl[j];
                                    }

                                    ptr_ch_B->local_mass += ptr_ch_A->ptr_cell_struct[box_idxNbr_ch_A].cell_mass;
                                    ptr_ch_A->local_mass -= ptr_ch_A->ptr_cell_struct[box_idxNbr_ch_A].cell_mass;
                                    ptr_ch_A->ptr_cell_struct[box_idxNbr_ch_A].cell_mass = 0;
                                    ptr_ch_A->ptr_cell_struct[box_idxNbr_ch_A].ptcl_size = 0;
                                }
                            }
                        }
                    }
                }
            }
            ptr_ch_A->cell_size = no_cells_ch_A;
        }
    }


    printf("C\n");
    //Case of old child nodes that will not be used
    cntr_ch = 0;
    for (int ch = 0; ch < ptr_node->chn_size - ptr_node->zones_size; ch++)
    {
        while (links_old_ord_old[cntr_ch] == (cntr_ch + ch) && cntr_ch < ptr_node->zones_size)
        {
            cntr_ch++;
        }

        ptr_ch_A = ptr_node->pptr_chn[cntr_ch + ch];

        no_cells_ch_A = ptr_ch_A->cell_size;
        for (int cell_idx = 0; cell_idx < no_cells_ch_A; cell_idx += 8)
        {
            box_idx_x_node = (ptr_ch_A->ptr_cell_idx_x[cell_idx] >> 1) - ptr_node->box_ts_x;
            box_idx_y_node = (ptr_ch_A->ptr_cell_idx_y[cell_idx] >> 1) - ptr_node->box_ts_y;
            box_idx_z_node = (ptr_ch_A->ptr_cell_idx_z[cell_idx] >> 1) - ptr_node->box_ts_z;
            box_idx_node = box_idx_x_node + box_idx_y_node * ptr_node->box_real_dim_x + box_idx_z_node * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;
            if (ptr_node->ptr_box_aux[box_idx_node] >= 0)
            {
                ptr_ch_B = ptr_node->pptr_chn[links_old_ord_new[ptr_node->ptr_box_aux[box_idx_node]]];

                //** >> CELLS **/
                no_cells_ch_B = ptr_ch_B->cell_size;

                box_idx_x_ch_A = ptr_ch_A->ptr_cell_idx_x[cell_idx] - ptr_ch_A->box_ts_x;
                box_idx_y_ch_A = ptr_ch_A->ptr_cell_idx_y[cell_idx] - ptr_ch_A->box_ts_y;
                box_idx_z_ch_A = ptr_ch_A->ptr_cell_idx_z[cell_idx] - ptr_ch_A->box_ts_z;
                box_idx_ch_A = box_idx_x_ch_A + box_idx_y_ch_A * ptr_ch_A->box_real_dim_x + box_idx_z_ch_A * ptr_ch_A->box_real_dim_x * ptr_ch_A->box_real_dim_y;

                box_idx_x_ch_B = ptr_ch_A->ptr_cell_idx_x[cell_idx] - ptr_ch_B->box_ts_x;
                box_idx_y_ch_B = ptr_ch_A->ptr_cell_idx_y[cell_idx] - ptr_ch_B->box_ts_y;
                box_idx_z_ch_B = ptr_ch_A->ptr_cell_idx_z[cell_idx] - ptr_ch_B->box_ts_z;
                box_idx_ch_B = box_idx_x_ch_B + box_idx_y_ch_B * ptr_ch_B->box_real_dim_x + box_idx_z_ch_B * ptr_ch_B->box_real_dim_x * ptr_ch_B->box_real_dim_y;


                for (int kk = 0; kk < 2; kk++)
                {
                    for (int jj = 0; jj < 2; jj++)
                    {
                        for (int ii = 0; ii < 2; ii++)
                        {
                            aux_const = ii + jj * 2 + kk * 4;
                            ptr_ch_B->ptr_cell_idx_x[no_cells_ch_B + aux_const] = ptr_ch_A->ptr_cell_idx_x[cell_idx + aux_const];
                            ptr_ch_B->ptr_cell_idx_y[no_cells_ch_B + aux_const] = ptr_ch_A->ptr_cell_idx_y[cell_idx + aux_const];
                            ptr_ch_B->ptr_cell_idx_z[no_cells_ch_B + aux_const] = ptr_ch_A->ptr_cell_idx_z[cell_idx + aux_const];

                            // Updating the box status of child node B
                            box_idxNbr_ch_B = box_idx_ch_B + ii + jj * ptr_ch_B->box_real_dim_x + kk * ptr_ch_B->box_real_dim_x * ptr_ch_B->box_real_dim_y;
                            ptr_ch_B->ptr_box[box_idx_ch_B] = -3; // Putting the status of EXIST (-3) in the child node cell
                        }
                    }
                }
                ptr_ch_B->cell_size += 8;

                //** >> Cell structure and local mass **/
                for (int kk = 0; kk < 2; kk++)
                {
                    for (int jj = 0; jj < 2; jj++)
                    {
                        for (int ii = 0; ii < 2; ii++)
                        {
                            box_idxNbr_ch_A = box_idx_ch_A + ii + jj * ptr_ch_A->box_real_dim_x + kk * ptr_ch_A->box_real_dim_x * ptr_ch_A->box_real_dim_y;
                            box_idxNbr_ch_B = box_idx_ch_B + ii + jj * ptr_ch_B->box_real_dim_x + kk * ptr_ch_B->box_real_dim_x * ptr_ch_B->box_real_dim_y;

                            ptr_ch_B->ptr_cell_struct[box_idxNbr_ch_B].cell_mass = ptr_ch_A->ptr_cell_struct[box_idxNbr_ch_A].cell_mass;
                            ptr_ch_B->ptr_cell_struct[box_idxNbr_ch_B].ptcl_size = ptr_ch_A->ptr_cell_struct[box_idxNbr_ch_A].ptcl_size;

                            //** >> Space checking of cell structure particles
                            if (space_check(&(ptr_ch_B->ptr_cell_struct[box_idxNbr_ch_B].ptcl_cap), ptr_ch_B->ptr_cell_struct[box_idxNbr_ch_B].ptcl_size, 2.0f, "p1i1", &(ptr_ch_B->ptr_cell_struct[box_idxNbr_ch_B].ptr_ptcl)) == _FAILURE_)
                            {
                                printf("Error, in space_check function\n");
                                return _FAILURE_;
                            }

                            for (int j = 0; j < ptr_ch_B->ptr_cell_struct[box_idxNbr_ch_B].ptcl_size; j++)
                            {
                                ptr_ch_B->ptr_cell_struct[box_idxNbr_ch_B].ptr_ptcl[j] = ptr_ch_A->ptr_cell_struct[box_idxNbr_ch_A].ptr_ptcl[j];
                            }

                            ptr_ch_B->local_mass += ptr_ch_A->ptr_cell_struct[box_idxNbr_ch_A].cell_mass;
                        }
                    }
                }
                
            }
        }
    }
    return _SUCCESS_;
}

static int moving_new_zones_to_new_child(struct node *ptr_node, int *links_old_ord_new)
{
    struct node *ptr_ch; // child node

    int box_idx_x_node; // Box index in X direcction of the node cell
    int box_idx_y_node; // Box index in Y direcction of the node cell
    int box_idx_z_node; // Box index in Z direcction of the node cell
    int box_idx_node;   // Box index of the node cell

    int box_idx_x_ch; // Box index in X direcction of the node cell
    int box_idx_y_ch; // Box index in Y direcction of the node cell
    int box_idx_z_ch; // Box index in Z direcction of the node cell
    int box_idx_ch;   // Box index of the node cell

    int box_idx_ptcl_ch;   // Box index of particles in the child node

    int cell_idx_x_aux;
    int cell_idx_y_aux;
    int cell_idx_z_aux;

    int cell_idx_node;
    int cell_idx_ch;

    int ptcl_idx;

    int no_cells_ch;    //Number of cell in the child node



    for (int zone_idx = 0; zone_idx < ptr_node->zones_size; zone_idx++)
    {
        ptr_ch = ptr_node->pptr_chn[links_old_ord_new[zone_idx]];

        no_cells_ch = ptr_ch->cell_size;
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
                //** >> CELLS, BOX STATUS AND CAPACITY OF CELL STRUCTURE **/
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

                            // Cell structure capacity
                            //** >> Space checking of the capacity of the particles in the child cells **/
                            if (space_check(&(ptr_ch->ptr_cell_struct[box_idx_ch].ptcl_cap), ptr_node->ptr_cell_struct[box_idx_node].ptcl_size + 1, 2.0f, "p1i1", &(ptr_ch->ptr_cell_struct[box_idx_ch].ptr_ptcl)) == _FAILURE_)
                            {
                                printf("Error, in space_check function\n");
                                return _FAILURE_;
                            }
                        }
                    }
                }
                no_cells_ch += 8;

                //Transfer particles from node to child
                for (int j = 0; j < ptr_node->ptr_cell_struct[box_idx_node].ptcl_size; j++)
                {
                    ptcl_idx = ptr_node->ptr_cell_struct[box_idx_node].ptr_ptcl[j];
                    box_idx_ptcl_ch = ptcl_idx_to_box_idx(ptr_ch, ptcl_idx);
                    ptr_ch->ptr_cell_struct[box_idx_ptcl_ch].cell_mass += GL_ptcl_mass[ptcl_idx];
                    ptr_ch->local_mass += GL_ptcl_mass[ptcl_idx];
                    ptr_ch->ptr_cell_struct[box_idx_ptcl_ch].ptr_ptcl[ptr_ch->ptr_cell_struct[box_idx_ptcl_ch].ptcl_size] = ptcl_idx;
                    ptr_ch->ptr_cell_struct[box_idx_ptcl_ch].ptcl_size += 1;
                }
            }
        }
        ptr_ch->cell_size = no_cells_ch;
    }
    return _SUCCESS_;
}

static void update_border_child_boxes(struct node *ptr_node, const int *links_old_ord_old)
{

    struct node *ptr_ch;

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
}

static void reorganization_child_node(struct node *ptr_node, const int *links_old_ord_old, const int *links_new_ord_old, const int *links_old_ord_new)
{
    struct node *ptr_ch;
    struct node **pptr_aux;

    int cntr_ch;

    pptr_aux = (struct node **)malloc((links_old_ord_old[ptr_node->zones_size - 1] + 1) * sizeof(struct node *));

    // Updating the ID of the child nodes
    for (int zone_idx = 0; zone_idx < ptr_node->zones_size; zone_idx++)
    {
        ptr_ch = ptr_node->pptr_chn[links_old_ord_old[zone_idx]];
        ptr_ch->ID = links_new_ord_old[zone_idx];
    }

    // Reorganization of child nodes
    for (int i = 0; i < ptr_node->zones_size; i++)
    {
        pptr_aux[i] = ptr_node->pptr_chn[links_old_ord_new[i]];
    }

    cntr_ch = 0;
    for (int i = 0; i < links_old_ord_old[ptr_node->zones_size - 1] + 1 - ptr_node->zones_size; i++)
    {
        while (links_old_ord_old[cntr_ch] == cntr_ch + i)
        {
            cntr_ch++;
        }
        pptr_aux[ptr_node->zones_size + i] = ptr_node->pptr_chn[cntr_ch + i];
    }

    for (int i = 0; i < links_old_ord_old[ptr_node->zones_size - 1] + 1; i++)
    {
        ptr_node->pptr_chn[i] = pptr_aux[i];
    }

    int no_grandchildren2 = 0;
    for (int ch = 0; ch < ptr_node->zones_size; ch++)
    {
        for (int grandch = 0; grandch < ptr_node->pptr_chn[ch]->chn_size; grandch++)
        {
            no_grandchildren2++;
        }
    }
    for (int ch = ptr_node->zones_size; ch < ptr_node->chn_size; ch++)
    {
        for (int grandch = 0; grandch < ptr_node->pptr_chn[ch]->chn_size; grandch++)
        {
            no_grandchildren2++;
        }
    }

    free(pptr_aux);
}

static int reorganization_grandchild_node(struct node *ptr_node)
{
    int no_grandchildren = 0;

    for (int ch = 0; ch < ptr_node->zones_size; ch++)
    {
        for (int grandch = 0; grandch < ptr_node->pptr_chn[ch]->chn_size; grandch++)
        {
            no_grandchildren++;
        }
    }
    for (int ch = ptr_node->zones_size; ch < ptr_node->chn_size; ch++)
    {
        for (int grandch = 0; grandch < ptr_node->pptr_chn[ch]->chn_size; grandch++)
        {
            no_grandchildren++;
        }
    }

    if(no_grandchildren > 0)
    {
        int box_idx_x_node; // Box index in X direcction of the node cell
        int box_idx_y_node; // Box index in Y direcction of the node cell
        int box_idx_z_node; // Box index in Z direcction of the node cell
        int box_idx_node;   // Box index of the node cell

        int child_ID;

        int size;

        struct node **pptr_aux;
        int cntr_ch = 0; //Counter the number of grandchildren
        pptr_aux = (struct node **)malloc(no_grandchildren * sizeof(struct node *));

        for (int ch = 0; ch < ptr_node->zones_size; ch++)
        {
            for (int grandch = 0; grandch < ptr_node->pptr_chn[ch]->chn_size; grandch++)
            {
                pptr_aux[cntr_ch] = ptr_node->pptr_chn[ch]->pptr_chn[grandch];
                ptr_node->pptr_chn[ch]->pptr_chn[grandch] = NULL;
                cntr_ch++;
            }
            ptr_node->pptr_chn[ch]->chn_size = 0; // Reseting the number of grandchildren to 0
        }
        for (int ch = ptr_node->zones_size; ch < ptr_node->chn_size; ch++)
        {
            for (int grandch = 0; grandch < ptr_node->pptr_chn[ch]->chn_size; grandch++)
            {
                pptr_aux[cntr_ch] = ptr_node->pptr_chn[ch]->pptr_chn[grandch];
                ptr_node->pptr_chn[ch]->pptr_chn[grandch] = NULL;
                cntr_ch++;
            }
            ptr_node->pptr_chn[ch]->chn_size = 0; // Reseting the number of grandchildren to 0
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

            if (space_check(&(ptr_node->pptr_chn[child_ID]->chn_cap), size + 1, 2.0f, "p1n2", &(ptr_node->pptr_chn[child_ID]->pptr_chn)) == _FAILURE_)
            {
                printf("Error, in space_check function\n");
                return _FAILURE_;
            }
            ptr_node->pptr_chn[child_ID]->pptr_chn[size] = pptr_aux[grandch];
            pptr_aux[grandch]->ptr_pt = ptr_node->pptr_chn[child_ID];
            pptr_aux[grandch]->ID = size;
            ptr_node->pptr_chn[child_ID]->chn_size = size + 1;
        }
        free(pptr_aux);
    }

    return _SUCCESS_;
}

static void moved_unused_child_node_to_memory_pool(struct node *ptr_node)
{
    //** >> Putting unused child nodes to the stack of memory pool **/
    for (int i = ptr_node->zones_size; i < ptr_node->chn_size; i++)
    {
        add_node_to_stack(ptr_node->pptr_chn[i]);
    }
}

static void updating_ref_zones_grandchildren(struct node *ptr_node)
{
    struct node *ptr_ch;
    struct node *ptr_grandch;

    int box_idx_x_ch; // Box index in X direcction of the child cell
    int box_idx_y_ch; // Box index in Y direcction of the child cell
    int box_idx_z_ch; // Box index in Z direcction of the child cell
    int box_idx_ch;   // Box index of the child cell

    for (int ch = 0; ch < ptr_node->zones_size; ch++)
    {
        ptr_ch = ptr_node->pptr_chn[ch];
        for (int grandch = 0; grandch < ptr_ch->chn_size; grandch++)
        {
            ptr_grandch = ptr_ch->pptr_chn[grandch];
            for (int cell_idx = 0; cell_idx < ptr_grandch->cell_size; cell_idx +=8)
            {
                box_idx_x_ch = (ptr_grandch->ptr_cell_idx_x[cell_idx] >> 1) - ptr_ch->box_ts_x;
                box_idx_y_ch = (ptr_grandch->ptr_cell_idx_y[cell_idx] >> 1) - ptr_ch->box_ts_y;
                box_idx_z_ch = (ptr_grandch->ptr_cell_idx_z[cell_idx] >> 1) - ptr_ch->box_ts_z;
                box_idx_ch = box_idx_x_ch + box_idx_y_ch * ptr_ch->box_real_dim_x + box_idx_z_ch * ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;
                ptr_ch->ptr_box[box_idx_ch] = grandch;
            }
        }
    }
}

static int update_child_grid_points (struct node *ptr_node)
{

    struct node *ptr_ch;

    int box_idxNbr_i0_j0_k0_ch;    // Box index of the child node at x=0,y=0,z=0
    int box_idxNbr_im1_j0_k0_ch;   // Box index of the child node of the neighbor at x=-1,y=0,z=0
    int box_idxNbr_i0_jm1_k0_ch;   // Box index of the child node of the neighbor at x=0,y=-1,z=0
    int box_idxNbr_im1_jm1_k0_ch;  // Box index of the child node of the neighbor at x=-1,y=-1,z=0
    int box_idxNbr_i0_j0_km1_ch;   // Box index of the child node of the neighbor at x=0,y=0,z=-1
    int box_idxNbr_im1_j0_km1_ch;  // Box index of the child node of the neighbor at x=-1,y=0,z=-1
    int box_idxNbr_i0_jm1_km1_ch;  // Box index of the child node of the neighbor at x=0,y=-1,z=-1
    int box_idxNbr_im1_jm1_km1_ch; // Box index of the child node of the neighbor at x=-1,y=-1,z=-1

    int box_idx_ch; // Box index of the child node

    int box_grid_idx_ch; // Box grid index of the child node

    bool is_bder_grid_point; // Ask if the grid point is interior

    for(int ch = 0; ch < ptr_node->zones_size; ch++)
    {
        ptr_ch = ptr_node->pptr_chn[ch];
        //** >> Grid points **/
        for (int kk = ptr_ch->box_min_z - ptr_ch->box_ts_z; kk < ptr_ch->box_max_z - ptr_ch->box_ts_z + 2; kk++)
        {
            for (int jj = ptr_ch->box_min_y - ptr_ch->box_ts_y; jj < ptr_ch->box_max_y - ptr_ch->box_ts_y + 2; jj++)
            {
                for (int ii = ptr_ch->box_min_x - ptr_ch->box_ts_x; ii < ptr_ch->box_max_x - ptr_ch->box_ts_x + 2; ii++)
                {
                    box_grid_idx_ch = ii + jj * (ptr_ch->box_real_dim_x + 1) + kk * (ptr_ch->box_real_dim_x + 1) * (ptr_ch->box_real_dim_y + 1);

                    box_idx_ch = ii + jj * ptr_ch->box_real_dim_x + kk * ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;
                    box_idxNbr_i0_j0_k0_ch = box_idx_ch;
                    box_idxNbr_im1_j0_k0_ch = box_idx_ch - 1;
                    box_idxNbr_i0_jm1_k0_ch = box_idx_ch - ptr_ch->box_real_dim_x;
                    box_idxNbr_im1_jm1_k0_ch = box_idx_ch - 1 - ptr_ch->box_real_dim_x;
                    box_idxNbr_i0_j0_km1_ch = box_idx_ch - ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;
                    box_idxNbr_im1_j0_km1_ch = box_idx_ch - 1 - ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;
                    box_idxNbr_i0_jm1_km1_ch = box_idx_ch - ptr_ch->box_real_dim_x - ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;
                    box_idxNbr_im1_jm1_km1_ch = box_idx_ch - 1 - ptr_ch->box_real_dim_x - ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;

                    //** >> The grid point exist **/
                    if (ptr_ch->ptr_box[box_idxNbr_i0_j0_k0_ch] > -4 ||
                        ptr_ch->ptr_box[box_idxNbr_im1_j0_k0_ch] > -4 ||
                        ptr_ch->ptr_box[box_idxNbr_i0_jm1_k0_ch] > -4 ||
                        ptr_ch->ptr_box[box_idxNbr_im1_jm1_k0_ch] > -4 ||
                        ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] > -4 ||
                        ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] > -4 ||
                        ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] > -4 ||
                        ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] > -4)
                    {

                        is_bder_grid_point = false;
                        //** Connection to the right  **/
                        if (ptr_ch->ptr_box[box_idxNbr_i0_j0_k0_ch] < -3 &&
                            ptr_ch->ptr_box[box_idxNbr_i0_jm1_k0_ch] < -3 &&
                            ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] < -3 &&
                            ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] < -3)
                        {
                            is_bder_grid_point = true;
                        }
                        //** Connection to the left  **/
                        else if (ptr_ch->ptr_box[box_idxNbr_im1_j0_k0_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_jm1_k0_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] < -3)
                        {
                            is_bder_grid_point = true;
                        }
                        //** Backward connection   **/
                        else if (ptr_ch->ptr_box[box_idxNbr_i0_j0_k0_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_j0_k0_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] < -3)
                        {
                            is_bder_grid_point = true;
                        }
                        //** Forward connection   **/
                        else if (ptr_ch->ptr_box[box_idxNbr_i0_jm1_k0_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_jm1_k0_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] < -3)
                        {
                            is_bder_grid_point = true;
                        }
                        //** Upward connection **/
                        else if (ptr_ch->ptr_box[box_idxNbr_i0_j0_k0_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_j0_k0_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_i0_jm1_k0_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_jm1_k0_ch] < -3)
                        {
                            is_bder_grid_point = true;
                        }
                        //** Down connection **/
                        else if (ptr_ch->ptr_box[box_idxNbr_i0_j0_km1_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_j0_km1_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_i0_jm1_km1_ch] < -3 &&
                                 ptr_ch->ptr_box[box_idxNbr_im1_jm1_km1_ch] < -3)
                        {
                            is_bder_grid_point = true;
                        }

                        //** >> Adding the grid point **/

                        //** >> Border grid point**/
                        if (is_bder_grid_point == true)
                        {
                            //** >> Space checking of border grid points array**/
                            if (space_check(&(ptr_ch->grid_bder_cap), ptr_ch->grid_bder_size + 1, 2.0f, "p1i1", &(ptr_ch->ptr_grid_bder)) == _FAILURE_)
                            {
                                printf("Error, in space_check function\n");
                                return _FAILURE_;
                            }

                            //** >> Adding the grid point to the border array **/
                            ptr_ch->ptr_grid_bder[ptr_ch->grid_bder_size] = box_grid_idx_ch;
                            ptr_ch->grid_bder_size += 1; // Increasing the number of border grid points in the array
                        }
                        //** Interior grid point **/
                        else
                        {
                            //** >> Space checking of interior grid points array**/
                            if (space_check(&(ptr_ch->grid_intr_cap), ptr_ch->grid_intr_size + 1, 2.0f, "p1i1", &(ptr_ch->ptr_grid_intr)) == _FAILURE_)
                            {
                                printf("Error, in space_check function\n");
                                return _FAILURE_;
                            }

                            //** >> Adding the grid point to the interior array **/
                            ptr_ch->ptr_grid_intr[ptr_ch->grid_intr_size] = box_grid_idx_ch;
                            ptr_ch->grid_intr_size += 1; // Increasing the number of interior grid points in the array
                        }
                    }
                }
            }
        }
    }

    return _SUCCESS_;
}

static void exchange_box_aux_to_box(struct node *ptr_node)
{
    int *ptr_aux;

    ptr_aux = ptr_node->ptr_box;
    ptr_node->ptr_box = ptr_node->ptr_box_aux;
    ptr_node->ptr_box_aux = ptr_aux;

}

static int tentacles_updating(struct node *ptr_node, int tentacle_lv)
{
    int no_tentacles = GL_tentacles_size[tentacle_lv];
    int size = no_tentacles + ptr_node->zones_size;

    //** >> Space checking of the capacity of the refined cells **/
    if (space_check(&(GL_tentacles_cap[tentacle_lv]), size, 4.0f, "p1n2", &(GL_tentacles[tentacle_lv])) == _FAILURE_)
    {
        printf("Error, in space_check function\n");
        return _FAILURE_;
    }

    for (int i = 0; i < ptr_node->zones_size; i++)
    {
        GL_tentacles[tentacle_lv][no_tentacles + i] = ptr_node->pptr_chn[i];
    }

    GL_tentacles_size[tentacle_lv] = size;
    return _SUCCESS_;
}

static void update_chn_size(struct node *ptr_node)
{
    ptr_node->chn_size = ptr_node->zones_size;
}

static void updating_tentacles_max_lv()
{
    int no_lvs = GL_tentacles_level_max < (lmax - lmin) ? GL_tentacles_level_max : (GL_tentacles_level_max - 1);

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

        clock_t aux_clock;

        struct node *ptr_node;
        int no_pts; // Number of parents in the cycle
        int no_lvs; // Number of level of refinement to adapt

        int *links_old_ord_old; // Storing the old zone of refinement id using old order
        int *links_new_ord_old; // Storing the new zone of refinement id using old order
        int *links_old_ord_new; // Storing the old zone of refinement id using new order
        int *links_new_ord_new; // Storing the new zone of refinement id using new order
        int links_cap;

        links_cap = 256;
        links_old_ord_old = (int *) malloc (links_cap * sizeof(int));
        links_new_ord_old = (int *) malloc (links_cap * sizeof(int));
        links_old_ord_new = (int *) malloc (links_cap * sizeof(int));
        links_new_ord_new = (int *) malloc (links_cap * sizeof(int));

        no_lvs = GL_tentacles_level_max < (lmax - lmin) ? GL_tentacles_level_max : (GL_tentacles_level_max - 1);

        for (int lv = no_lvs; lv > -1; lv--)
        {
            GL_tentacles_size[lv + 1] = 0;
            no_pts = GL_tentacles_size[lv];
            //** >> For cycle over parent nodes **/
            for (int i = 0; i < no_pts; i++)
            {
                ptr_node = GL_tentacles[lv][i];

                //** Updating the box mass information **/
                printf("\n\nupdating_cell_struct\n\n");
                aux_clock = clock();
                if (updating_cell_struct(ptr_node) == _FAILURE_)
                {
                    printf("Error at function updating_cell_struct()\n");
                    return _FAILURE_;
                }
                GL_times[30] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

                //** Initialization of the box_aux **/
                printf("\n\nInitialization box aux\n\n");
                aux_clock = clock();
                initialization_box_aux(ptr_node);
                GL_times[31] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

                //** Initialization of the auiliary refinement arrays**/
                printf("\n\nInitialization ref aux\n\n\n");
                aux_clock = clock();
                initialization_ref_aux(ptr_node);
                GL_times[32] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;
                
                //** >> Filling the refinement cells array **/
                printf("\n\nFill cell ref\n\n\n");
                aux_clock = clock();
                if (fill_cell_ref(ptr_node) == _FAILURE_)
                {
                    printf("Error at function fill_cell_ref()\n");
                    return _FAILURE_;
                }
                GL_times[33] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

                //** >> Filling the different zones of refinement **/
                printf("\n\nFill zones ref\n\n");
                aux_clock = clock();
                if (fill_zones_ref(ptr_node) == _FAILURE_)
                {
                    printf("Error at function fill_zones_ref()\n");
                    return _FAILURE_;
                }
                GL_times[34] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

                //** >> Create links **/
                printf("\n\nCreate links\n\n");
                aux_clock = clock();
                if(ptr_node->zones_size > 0)
                {
                    if (create_links(ptr_node, &links_old_ord_old, &links_new_ord_old, &links_old_ord_new, &links_new_ord_new, &links_cap) == _FAILURE_)
                    {
                        printf("Error at function create_links()\n");
                        return _FAILURE_;
                    }
                }
                GL_times[36] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

                //** >> Removing cells that no longer require refinement **/
                printf("\n\nRemoving cells nolonger requiere refinement\n\n");
                aux_clock = clock();
                if (0 < ptr_node->zones_size && 0 < ptr_node->chn_size)
                {
                    remov_cells_nolonger_require_refinement(ptr_node, links_old_ord_old, links_new_ord_old);
                }
                GL_times[37] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

                //** >> Adapting child boxes to the new space **/
                printf("\n\nAdapt child box\n\n");
                aux_clock = clock();
                if (0 < ptr_node->zones_size && 0 < ptr_node->chn_size)
                {
                    if (adapt_child_box_and_cells(ptr_node, links_old_ord_old, links_new_ord_old) == _FAILURE_)
                    {
                        printf("Error at function adapt_child_box()\n");
                        return _FAILURE_;
                    }
                }
                GL_times[38] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;
                
                //** >> Creating and defining new children nodes for excess in refinement zones and and linking them to the parent node ptr_node **/
                printf("\n\nCreate new child nodes\n\n");
                aux_clock = clock();
                if (ptr_node->zones_size > ptr_node->chn_size)
                {

                    if (create_new_child_nodes(ptr_node, links_old_ord_old, links_new_ord_old) == _FAILURE_)
                    {
                        printf("Error, in new_child_nodes function\n");
                        return _FAILURE_;
                    }
                }
                GL_times[39] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

                //** >> Adapting the information from old child nodes to new child nodes **/
                printf("\n\nMoving old child to new child\n\n");
                aux_clock = clock();
                if (0 < ptr_node->zones_size && 0 < ptr_node->chn_size)
                {
                    if (moving_old_child_to_new_child(ptr_node, links_old_ord_old, links_new_ord_old, links_old_ord_new, links_new_ord_new) == _FAILURE_)
                    {
                        printf("Error at function moving_old_child_to_new_child()\n");
                        return _FAILURE_;
                    }
                }
                GL_times[40] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

                //** >> Moving new zones of refienemnt information to all child nodes **/
                printf("\n\nMoving new zones to new child\n\n");
                aux_clock = clock();
                if (0 < ptr_node->zones_size)
                {
                    if (moving_new_zones_to_new_child(ptr_node, links_old_ord_new) == _FAILURE_)
                    {
                        printf("Error at function moving_new_zones_to_new_child()\n");
                        return _FAILURE_;
                    }
                }
                GL_times[41] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

                //** >> Update border of the child boxes **//
                printf("\n\nUpdate border child boxes\n\n");
                aux_clock = clock();
                if (ptr_node->zones_size > 0)
                {
                    update_border_child_boxes(ptr_node, links_old_ord_old);
                }
                GL_times[42] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

                //** >> Reorganization child nodes **/
                printf("\n\nReorganization child nodes\n\n");
                aux_clock = clock();
                if (ptr_node->zones_size > 0)
                {
                    reorganization_child_node(ptr_node, links_old_ord_old, links_new_ord_old, links_old_ord_new);
                }
                GL_times[43] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

                //** >> Reorganization grandchild nodes **/
                printf("\n\nReorganization grandchild nodes\n\n");
                aux_clock = clock();
                if (ptr_node->zones_size > 0 && ptr_node->chn_size > 0 && lv < no_lvs)
                {
                    if (reorganization_grandchild_node(ptr_node) == _FAILURE_)
                    {
                        printf("Error at function reorganization_grandchild_node()\n");
                        return _FAILURE_;
                    }
                }
                GL_times[44] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

                //** >> Moved Unused child node to the stack of memory pool **/
                printf("\n\nMoved Unused child node to the stack of memory pool\n\n");
                aux_clock = clock();
                if (ptr_node->chn_size > ptr_node->zones_size)
                {
                    moved_unused_child_node_to_memory_pool(ptr_node);
                }
                GL_times[45] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

                //** >> Updating refinement zones of the grandchildren **/
                printf("\n\nUpdating refinement zones of the grandchildren\n\n");
                aux_clock = clock();
                if (ptr_node->zones_size > 0 && ptr_node->chn_size > 0 && lv < no_lvs)
                {
                    updating_ref_zones_grandchildren(ptr_node);
                }
                GL_times[46] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

                //** >> Updating children grid points **/
                printf("\n\nUpdating children grid points\n\n");
                aux_clock = clock();
                if (ptr_node->zones_size > 0)
                {
                    if (update_child_grid_points(ptr_node) == _FAILURE_)
                    {
                        printf("Error at function update_child_grid_points()\n");
                        return _FAILURE_;
                    }
                }
                GL_times[47] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

                //** >> Exchange between auiliary box and box **/
                printf("\n\nExchange between auiliary box and box\n\n");
                aux_clock = clock();
                if (ptr_node->zones_size > 0 || ptr_node->chn_size > 0)
                {
                    exchange_box_aux_to_box(ptr_node);
                }
                GL_times[48] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

                //** >> Tentacles updating **/
                printf("\n\nTentacles Updating\n\n");
                aux_clock = clock();
                if (0 < ptr_node->zones_size)
                {
                    if (tentacles_updating(ptr_node, lv + 1) == _FAILURE_)
                    {
                        printf("Error at function tentacles_updating()\n");
                        return _FAILURE_;
                    }
                }
                GL_times[49] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

                //** >> Updating children size
                printf("\n\nUpdating chn size\n\n");
                aux_clock = clock();
                update_chn_size(ptr_node);
                GL_times[50] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;
            }
        }

        //** >> Tentacles Updating lv max **/
        //printf("\n\nTentacles updating lv max\n\n");
        aux_clock = clock();
        updating_tentacles_max_lv();
        GL_times[51] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

        aux_clock = clock();
        free(links_old_ord_old);
        free(links_new_ord_old);
        free(links_old_ord_new);
        free(links_new_ord_new);
        GL_times[52] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;
    }

        return _SUCCESS_;
}