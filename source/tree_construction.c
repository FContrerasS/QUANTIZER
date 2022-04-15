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

static int particle_distribution_in_cell(int **pptr_cell_ptcl, int *ptr_cell_ptcl_cap, int *ptr_cell_ptcl_size, const struct node *ptr_node)
{
    //** >> Filling the cells with the particles: **/

    int box_idx_x;         // Box index in X direcction
    int box_idx_y;         // Box index in Y direcction
    int box_idx_z;         // Box index in Z direcction
    int box_idx;           // Box index
    int ptcl_idx;          // Particle index

    int lv; // Level of refinement

    int size;

    lv = ptr_node->lv; 
    size = ptr_node->box_real_dim_x * ptr_node->box_real_dim_y * ptr_node->box_real_dim_z;

    //** >> Initializing pointers of cell particles **/
    for (int i = 0; i < size; i++)
    {
        ptr_cell_ptcl_size[i] = 0;
    }


    //** >> Filling cell particles arrays **/
    for (int i = 0; i < ptr_node->ptcl_size; i++)
    {
        ptcl_idx = ptr_node->ptr_ptcl[i];
        box_idx_x = GL_ptcl_x[ptcl_idx] * (1 << lv) - ptr_node->box_ts_x; // Particle indexes in the level
        box_idx_y = GL_ptcl_y[ptcl_idx] * (1 << lv) - ptr_node->box_ts_y;
        box_idx_z = GL_ptcl_z[ptcl_idx] * (1 << lv) - ptr_node->box_ts_z;
        box_idx = box_idx_x + box_idx_y * ptr_node->box_real_dim_x + box_idx_z * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;

        // //** >> Space checking of the capacity of particles in the cells **/
        if (space_check(&(ptr_cell_ptcl_cap[box_idx]), ptr_cell_ptcl_size[box_idx] + 1, "p1i1", &(pptr_cell_ptcl[box_idx])) == _FAILURE_)
        {
            printf("Error, in space_check function\n");
            return _FAILURE_;
        }

        //** >> Placing the particle in the corresponding cell **/
        pptr_cell_ptcl[box_idx][ptr_cell_ptcl_size[box_idx]] = ptcl_idx; // putting the particle index ptcl_idx in the cell at position box_idx
        ptr_cell_ptcl_size[box_idx] += 1;                                // Increasing the already number of particles in the cell
    }

    return _SUCCESS_;

} // end function particle_distribution_in_cell

static void filling_mass_box(int **pptr_cell_ptcl, const int *ptr_cell_ptcl_size, struct node *ptr_node)
{

    //** >> Filling the mass box of the main node structure **/
    // int cell_idx; // The cell index is simply i of the for loop
    int ptcl_idx; // Particle index

    int size = ptr_node->box_real_dim_x * ptr_node->box_real_dim_y * ptr_node->box_real_dim_z;

    vtype total_mass = 0;

    //** >> Check all elements in the box **/
    for (int i = 0; i < size; i++)
    {
        //** >> Analizing BORDER cells (-2) and EXIST cells (-3) **/
        if (ptr_node->ptr_box[i] > -4)
        {
            for (int j = 0; j < ptr_cell_ptcl_size[i]; j++)
            {
                ptcl_idx = pptr_cell_ptcl[i][j];
                total_mass += GL_ptcl_mass[ptcl_idx];
                ptr_node->ptr_box_mass[i] += GL_ptcl_mass[ptcl_idx]; // The mass was initialized to zero
            }
        }
    }

} // end function filling_density_box_main_node

static int fill_cell_ref(struct node *ptr_node)
{
    //** >> Adding cells which satisfy the refinement criteria to the array ptr_cell_ref and chaning the box the status of refinement -1 **/
       
    int size;                         // Size of the refinement cells array

    int box_idx_x;  // Box index in X direcction
    int box_idx_y;  // Box index in Y direcction
    int box_idx_z;  // Box index in Z direcction
    int box_idx;    // Box index
    int box_idxNbr; // Box index in the neigborhood
    // int cell_idx; // The cell index is simply i of the for loop


    int cntr;     // Counter used to add cell index to the ptr_cell_ref

    //** >> Changing the box status from EXIST (-3) to REFINEMENT REQUIRED (-1) **/
    size = 0; // Initial size of the cell refined array must be 0
    for (int i = 0; i < ptr_node->cell_size; i++)
    {
        box_idx_x = ptr_node->ptr_cell_idx_x[i] - ptr_node->box_ts_x;
        box_idx_y = ptr_node->ptr_cell_idx_y[i] - ptr_node->box_ts_y;
        box_idx_z = ptr_node->ptr_cell_idx_z[i] - ptr_node->box_ts_z;
        box_idx = box_idx_x + box_idx_y * ptr_node->box_real_dim_x + box_idx_z * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;

        // Refinement criterion in the box_mass in no border box points
        if (ptr_node->ptr_box[box_idx] != -2 && ptr_node->ptr_box_mass[box_idx] >= ref_criterion_mass) // No border (-2)
        {
            if (ptr_node->ptr_box[box_idx] == -3) // Cell has not been added yet
            {
                size++;
                //** >> Chaning the cell box status from EXIST (-3) to REFINEMENT REQUIRED (-1) **/
                ptr_node->ptr_box[box_idx] = -1;
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
                        if (ptr_node->ptr_box[box_idxNbr] == -3 || ptr_node->ptr_box[box_idxNbr] == -2) // Cell has not been added yet
                        {
                            size++;
                            //** >> Chaning the cell box status from EXIST (-3) to REFINEMENT REQUIRED (-1) **/
                            ptr_node->ptr_box[box_idxNbr] = -1;
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

        if (ptr_node->ptr_box[box_idx] == -1) // Cell require refinement
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
    int cntr;     // Counter Number of cells added to any refinement zone
    int cntr2;  // Counter number of cells added per inspection
    int insp;     // Numer of inspected cells in the zone

    int zone_size;    // Number of cells in the zone
    int zone_idx_max; // Maximum id of the zone. It is equal to to the capacity in the number of zones
    int zone_idx;     // Index of the zone

    int box_idx_x;    // Box index in X direcction
    int box_idx_y;    // Box index in Y direcction
    int box_idx_z;    // Box index in Z direcction
    int box_idx;      // Box index

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

        if (ptr_node->ptr_box[box_idx] == -1) // A cell without zone has been founded
        {
            zone_size = 0; // Initial number of element in the zone

            //** >> Including the first element of the box to the auxiliary array ptr_aux_idx **/
            ptr_node->ptr_aux_idx[0] = box_idx;

            //** >>  Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0) **/
            ptr_node->ptr_box[box_idx] = zone_idx;

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
                if (ptr_node->ptr_box[box_idxNbr_x_plus] == -1)
                {
                    ptr_node->ptr_aux_idx[zone_size + cntr2] = box_idxNbr_x_plus; // Including the neighboring element of the box to the auxiliary array
                    ptr_node->ptr_box[box_idxNbr_x_plus] = zone_idx;          // Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0)
                    cntr2++;                                       // +1 to the number of cells added in the current inspection
                }
                // Second neighbor
                if (ptr_node->ptr_box[box_idxNbr_x_minus] == -1)
                {
                    ptr_node->ptr_aux_idx[zone_size + cntr2] = box_idxNbr_x_minus; // Including the neighboring element of the box to the auxiliary array
                    ptr_node->ptr_box[box_idxNbr_x_minus] = zone_idx;          // Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0)
                    cntr2++;                                       // +1 to the number of cells added in the current inspection
                }
                // Third neighbor
                if (ptr_node->ptr_box[box_idxNbr_y_plus] == -1)
                {
                    ptr_node->ptr_aux_idx[zone_size + cntr2] = box_idxNbr_y_plus; // Including the neighboring element of the box to the auxiliary array
                    ptr_node->ptr_box[box_idxNbr_y_plus] = zone_idx;          // Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0)
                    cntr2++;                                       // +1 to the number of cells added in the current inspection
                }
                // Fourth neighbor
                if (ptr_node->ptr_box[box_idxNbr_y_minus] == -1)
                {
                    ptr_node->ptr_aux_idx[zone_size + cntr2] = box_idxNbr_y_minus; // Including the neighboring element of the box to the auxiliary array
                    ptr_node->ptr_box[box_idxNbr_y_minus] = zone_idx;          // Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0)
                    cntr2++;                                       // +1 to the number of cells added in the current inspection
                }
                // Fifth neighbor
                if (ptr_node->ptr_box[box_idxNbr_z_plus] == -1)
                {
                    ptr_node->ptr_aux_idx[zone_size + cntr2] = box_idxNbr_z_plus; // Including the neighboring element of the box to the auxiliary array
                    ptr_node->ptr_box[box_idxNbr_z_plus] = zone_idx;          // Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0)
                    cntr2++;                                       // +1 to the number of cells added in the current inspection
                }
                // Sixth neighbor
                if (ptr_node->ptr_box[box_idxNbr_z_minus] == -1)
                {
                    ptr_node->ptr_aux_idx[zone_size + cntr2] = box_idxNbr_z_minus; // Including the neighboring element of the box to the auxiliary array
                    ptr_node->ptr_box[box_idxNbr_z_minus] = zone_idx;          // Changing the box status from REFINEMENT REQUIRED (-1) to the refinement zone ID (>= 0)
                    cntr2++;                                       // +1 to the number of cells added in the current inspection
                }
                zone_size += cntr2; // Increasing the number of cells in the zone
                cntr += cntr2;      // Increasing the number of cells added to the zones
                insp++;             // Increasing the number of inspections
            }// End while cycle, now the box contains the the information about all cells of the zone "zone_idx"

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
    } // At this point the box contains the information of all refinement zones

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
        zone_idx = ptr_node->ptr_box[box_idx];
        cntr = ptr_node->ptr_aux_idx[zone_idx];          // Counter the element in the zone "zone_idx"
        ptr_node->pptr_zones[zone_idx][cntr] = cell_idx; // Adding the index of the cell array in the block to the zone
        ptr_node->ptr_aux_idx[zone_idx] += 1;            // Counter the number of elements added in the zone "zone_idx"
    }

    return _SUCCESS_;
} // End function fill_zones_ref

static int fill_child_nodes(int **pptr_cell_ptcl, const int *ptr_cell_ptcl_size, struct node *ptr_node_pt)
{

    struct node *ptr_node_ch = NULL; // Pointer to the child node
    int pt_cells_no;                 // Number of parent cells in the child
    int *ptr_pt_cells = NULL;        // Parent cells in the child

    int cell_idx;  // Cell index
    int aux_idx_x; // Auxiliary indexes
    int aux_idx_y;
    int aux_idx_z;

    int box_idx_x;  // Box index in X direcction
    int box_idx_y;  // Box index in Y direcction
    int box_idx_z;  // Box index in Z direcction
    int box_idx;    // Box index
    int box_idxNbr; // Box index in the neigborhood

    int box_idxNbr_i0_j0_k0;  // Box index at x=0,y=0,z=0
    int box_idxNbr_im1_j0_k0; // Box index of the neighbor at x=-1,y=0,z=0
    int box_idxNbr_i0_jm1_k0; // Box index of the neighbor at x=0,y=-1,z=0
    int box_idxNbr_im1_jm1_k0; // Box index of the neighbor at x=-1,y=-1,z=0
    int box_idxNbr_i0_j0_km1; // Box index of the neighbor at x=0,y=0,z=-1
    int box_idxNbr_im1_j0_km1;   // Box index of the neighbor at x=-1,y=0,z=-1
    int box_idxNbr_i0_jm1_km1;  // Box index of the neighbor at x=0,y=-1,z=-1
    int box_idxNbr_im1_jm1_km1;  // Box index of the neighbor at x=-1,y=-1,z=-1

    int box_grid_idx; // Box grid index

    bool is_bder_grid_point; // Ask if the grid point is interior

    int cntr; // Counter
    int pos_x;  // Distance between the real box and the minimal box when the last is localized in the middle of the real one
    int pos_y;
    int pos_z;

    int cap;  // Capacity of arrays in child nodes
    int size; // Size or number of elements in some array in child nodes

    //** >> Creating child nodes in parent node **/
    ptr_node_pt->chn_cap = ptr_node_pt->zones_cap; // Same amount than refinement zones
    ptr_node_pt->chn_size = ptr_node_pt->zones_size;
    ptr_node_pt->pptr_chn = (struct node **)malloc(ptr_node_pt->chn_cap * sizeof(struct node *)); // Allocating children pointers
    for (int i = 0; i < ptr_node_pt->chn_cap; i++)
    {
        ptr_node_pt->pptr_chn[i] = NULL;
    }

    //** >> Filling child nodes **/
    for (int i = 0; i < ptr_node_pt->chn_size; i++)
    {
        //** >> Create and Initialize child node **/
        ptr_node_ch = (struct node *)malloc(sizeof(struct node));
        initialize_node(ptr_node_ch);

        //** >> Other variables **/
        pt_cells_no = ptr_node_pt->ptr_zone_size[i]; // Size or number of cells in the parent node zone
        ptr_pt_cells = ptr_node_pt->pptr_zones[i];   // Cells in the parent node zone

        //** >> Global prperties **/
        ptr_node_ch->ID = i;
        ptr_node_ch->lv = ptr_node_pt->lv + 1;

        //** >> Cells in the node **/
        cap = 2 * 8 * pt_cells_no; // Capacity is the double of total child cells
        ptr_node_ch->cell_cap = cap;
        size = 8 * pt_cells_no;
        ptr_node_ch->cell_size = size;
        ptr_node_ch->ptr_cell_idx_x = (int *)malloc(ptr_node_ch->cell_cap * sizeof(int));
        ptr_node_ch->ptr_cell_idx_y = (int *)malloc(ptr_node_ch->cell_cap * sizeof(int));
        ptr_node_ch->ptr_cell_idx_z = (int *)malloc(ptr_node_ch->cell_cap * sizeof(int));

        for (int j = 0; j < pt_cells_no; j++)
        {
            aux_idx_x = ptr_node_pt->ptr_cell_idx_x[ptr_pt_cells[j]] * 2;
            aux_idx_y = ptr_node_pt->ptr_cell_idx_y[ptr_pt_cells[j]] * 2;
            aux_idx_z = ptr_node_pt->ptr_cell_idx_z[ptr_pt_cells[j]] * 2;
            for (int kk = 0; kk < 2; kk++)
            {
                for (int jj = 0; jj < 2; jj++)
                {
                    for (int ii = 0; ii < 2; ii++)
                    {
                        cell_idx = ii + jj * 2 + kk * 4 + j * 8; // Child cell index
                        ptr_node_ch->ptr_cell_idx_x[cell_idx] = aux_idx_x + ii;
                        ptr_node_ch->ptr_cell_idx_y[cell_idx] = aux_idx_y + jj;
                        ptr_node_ch->ptr_cell_idx_z[cell_idx] = aux_idx_z + kk;
                    }
                }
            }
        }

        //** >> Particles in the node **/
        // Size and Capacity
        size = 0; // Number of particles in the zone
        for (int j = 0; j < ptr_node_pt->ptr_zone_size[i]; j++)
        {
            box_idx_x = ptr_node_pt->ptr_cell_idx_x[ptr_node_pt->pptr_zones[i][j]] - ptr_node_pt->box_ts_x;
            box_idx_y = ptr_node_pt->ptr_cell_idx_y[ptr_node_pt->pptr_zones[i][j]] - ptr_node_pt->box_ts_y;
            box_idx_z = ptr_node_pt->ptr_cell_idx_z[ptr_node_pt->pptr_zones[i][j]] - ptr_node_pt->box_ts_z;
            box_idx = box_idx_x + box_idx_y * ptr_node_pt->box_real_dim_x + box_idx_z * ptr_node_pt->box_real_dim_x * ptr_node_pt->box_real_dim_y;
            size += ptr_cell_ptcl_size[box_idx];
        }
        cap =  size; // Capacity in the number of particles in the node
        ptr_node_ch->ptcl_cap = cap;
        ptr_node_ch->ptcl_size = size;
        // Putting particles in the child node
        ptr_node_ch->ptr_ptcl = (int *)malloc(cap * sizeof(int));
        cntr = 0;
        for (int j = 0; j < ptr_node_pt->ptr_zone_size[i]; j++)
        {
            box_idx_x = ptr_node_pt->ptr_cell_idx_x[ptr_node_pt->pptr_zones[i][j]] - ptr_node_pt->box_ts_x;
            box_idx_y = ptr_node_pt->ptr_cell_idx_y[ptr_node_pt->pptr_zones[i][j]] - ptr_node_pt->box_ts_y;
            box_idx_z = ptr_node_pt->ptr_cell_idx_z[ptr_node_pt->pptr_zones[i][j]] - ptr_node_pt->box_ts_z;
            box_idx = box_idx_x + box_idx_y * ptr_node_pt->box_real_dim_x + box_idx_z * ptr_node_pt->box_real_dim_x * ptr_node_pt->box_real_dim_y;
            for (int k = 0; k < ptr_cell_ptcl_size[box_idx]; k++)
            {
                ptr_node_ch->ptr_ptcl[cntr] = pptr_cell_ptcl[box_idx][k];
                cntr++;
            }
        }

        //** >> Boxes **/
        // MIN and MAX cell indexes values of the node.
        for (int j = 0; j < ptr_node_pt->ptr_zone_size[i]; j++)
        {
            cell_idx = ptr_node_pt->pptr_zones[i][j]; // Cell index in the parent zone
            if (ptr_node_ch->box_min_x > ptr_node_pt->ptr_cell_idx_x[cell_idx])
            {
                ptr_node_ch->box_min_x = ptr_node_pt->ptr_cell_idx_x[cell_idx];
            }
            if (ptr_node_ch->box_min_y > ptr_node_pt->ptr_cell_idx_y[cell_idx])
            {
                ptr_node_ch->box_min_y = ptr_node_pt->ptr_cell_idx_y[cell_idx];
            }
            if (ptr_node_ch->box_min_z > ptr_node_pt->ptr_cell_idx_z[cell_idx])
            {
                ptr_node_ch->box_min_z = ptr_node_pt->ptr_cell_idx_z[cell_idx];
            }
            if (ptr_node_ch->box_max_x < ptr_node_pt->ptr_cell_idx_x[cell_idx])
            {
                ptr_node_ch->box_max_x = ptr_node_pt->ptr_cell_idx_x[cell_idx];
            }
            if (ptr_node_ch->box_max_y < ptr_node_pt->ptr_cell_idx_y[cell_idx])
            {
                ptr_node_ch->box_max_y = ptr_node_pt->ptr_cell_idx_y[cell_idx];
            }
            if (ptr_node_ch->box_max_z < ptr_node_pt->ptr_cell_idx_z[cell_idx])
            {
                ptr_node_ch->box_max_z = ptr_node_pt->ptr_cell_idx_z[cell_idx];
            }
        }

        // Chaning the min and max of the "minimal box" from parent units to child units
        ptr_node_ch->box_min_x = 2 * ptr_node_ch->box_min_x;
        ptr_node_ch->box_min_y = 2 * ptr_node_ch->box_min_y;
        ptr_node_ch->box_min_z = 2 * ptr_node_ch->box_min_z;
        ptr_node_ch->box_max_x = 2 * ptr_node_ch->box_max_x + 1;
        ptr_node_ch->box_max_y = 2 * ptr_node_ch->box_max_y + 1;
        ptr_node_ch->box_max_z = 2 * ptr_node_ch->box_max_z + 1;

        // Size of the "minimal box"
        ptr_node_ch->box_dim_x = ptr_node_ch->box_max_x - ptr_node_ch->box_min_x + 1;
        ptr_node_ch->box_dim_y = ptr_node_ch->box_max_y - ptr_node_ch->box_min_y + 1;
        ptr_node_ch->box_dim_z = ptr_node_ch->box_max_z - ptr_node_ch->box_min_z + 1;

        // Real dimensions of the box, or real capacity
        // ptr_node_ch->box_real_dim_x = fmax(ptr_node_ch->box_dim_x * 3, ptr_node_ch->box_dim_x + 2 * n_exp - 2);
        // ptr_node_ch->box_real_dim_y = fmax(ptr_node_ch->box_dim_y * 3, ptr_node_ch->box_dim_y + 2 * n_exp - 2);
        // ptr_node_ch->box_real_dim_z = fmax(ptr_node_ch->box_dim_z * 3, ptr_node_ch->box_dim_z + 2 * n_exp - 2);
        // ptr_node_ch->box_real_dim_x = fmax(ptr_node_ch->box_dim_x + 10, ptr_node_ch->box_dim_x + 2 * n_exp - 2);
        // ptr_node_ch->box_real_dim_y = fmax(ptr_node_ch->box_dim_y + 10, ptr_node_ch->box_dim_y + 2 * n_exp - 2);
        // ptr_node_ch->box_real_dim_z = fmax(ptr_node_ch->box_dim_z + 10, ptr_node_ch->box_dim_z + 2 * n_exp - 2);
        ptr_node_ch->box_real_dim_x = ptr_node_ch->box_dim_x + 10 > ptr_node_ch->box_dim_x + 2 * n_exp - 2 ? ptr_node_ch->box_dim_x + 10 : ptr_node_ch->box_dim_x + 2 * n_exp - 2;
        ptr_node_ch->box_real_dim_y = ptr_node_ch->box_dim_y + 10 > ptr_node_ch->box_dim_y + 2 * n_exp - 2 ? ptr_node_ch->box_dim_y + 10 : ptr_node_ch->box_dim_y + 2 * n_exp - 2;
        ptr_node_ch->box_real_dim_z = ptr_node_ch->box_dim_z + 10 > ptr_node_ch->box_dim_z + 2 * n_exp - 2 ? ptr_node_ch->box_dim_z + 10 : ptr_node_ch->box_dim_z + 2 * n_exp - 2;

        // Translations between cell array and box
        pos_x = (ptr_node_ch->box_real_dim_x - ptr_node_ch->box_dim_x) / 2; // Half of the distance of the box side less the "minimal box" side
        ptr_node_ch->box_ts_x = ptr_node_ch->box_min_x - pos_x;             // Every cell in the level l in the box must be subtracted this value to obtain the box index
        pos_y = (ptr_node_ch->box_real_dim_y - ptr_node_ch->box_dim_y) / 2;
        ptr_node_ch->box_ts_y = ptr_node_ch->box_min_y - pos_y;
        pos_z = (ptr_node_ch->box_real_dim_z - ptr_node_ch->box_dim_z) / 2;
        ptr_node_ch->box_ts_z = ptr_node_ch->box_min_z - pos_z;

        // Filling the box status
        cap = ptr_node_ch->box_real_dim_x * ptr_node_ch->box_real_dim_y * ptr_node_ch->box_real_dim_z; // In general, the size of each side must be 3 times bigger than the same side of the "minimal box"
        size = ptr_node_ch->cell_size;                                                                 // Number of cells in the block
        ptr_node_ch->ptr_box = (int *)malloc(cap * sizeof(int));
        ptr_node_ch->ptr_box_aux = (int *)malloc(cap * sizeof(int));
        // Putting the value of NO-EXIST (-4) in every box index
        for (int j = 0; j < cap; j++)
        {
            ptr_node_ch->ptr_box[j] = -4;
        }
        // Changing from NO-EXIST (-4) to EXIST (-3) to all cells in the block of the child node
        for (int j = 0; j < size; j++)
        {
            box_idx_x = ptr_node_ch->ptr_cell_idx_x[j] - ptr_node_ch->box_ts_x;
            box_idx_y = ptr_node_ch->ptr_cell_idx_y[j] - ptr_node_ch->box_ts_y;
            box_idx_z = ptr_node_ch->ptr_cell_idx_z[j] - ptr_node_ch->box_ts_z;
            box_idx = box_idx_x + box_idx_y * ptr_node_ch->box_real_dim_x + box_idx_z * ptr_node_ch->box_real_dim_x * ptr_node_ch->box_real_dim_y;
            ptr_node_ch->ptr_box[box_idx] = -3;
        }
        // Changing the EXIST (-3) TO BORDER (-2) to all border cells in the block
        for (int j = 0; j < size; j++)
        {
            box_idx_x = ptr_node_ch->ptr_cell_idx_x[j] - ptr_node_ch->box_ts_x;
            box_idx_y = ptr_node_ch->ptr_cell_idx_y[j] - ptr_node_ch->box_ts_y;
            box_idx_z = ptr_node_ch->ptr_cell_idx_z[j] - ptr_node_ch->box_ts_z;
            box_idx = box_idx_x + box_idx_y * ptr_node_ch->box_real_dim_x + box_idx_z * ptr_node_ch->box_real_dim_x * ptr_node_ch->box_real_dim_y;

            for (int kk = -1; kk < 2; kk++)
            {
                for (int jj = -1; jj < 2; jj++)
                {
                    for (int ii = -1; ii < 2; ii++)
                    {
                        box_idxNbr = box_idx + ii + jj * ptr_node_ch->box_real_dim_x + kk * ptr_node_ch->box_real_dim_x * ptr_node_ch->box_real_dim_y;
                        if (ptr_node_ch->ptr_box[box_idxNbr] == -4) // Border cells are those such that at least on of their first neighbors are NO-EXIST cells.
                        {
                            ptr_node_ch->ptr_box[box_idx] = -2;
                            ii = 2;
                            jj = 2;
                            kk = 2;
                        }
                    }
                }
            }
        }
        
        //** >> Grid points **/
        for (int kk = ptr_node_ch->box_min_z - ptr_node_ch->box_ts_z; kk < ptr_node_ch->box_max_z - ptr_node_ch->box_ts_z + 2; kk++)
        {
            for (int jj = ptr_node_ch->box_min_y - ptr_node_ch->box_ts_y; jj < ptr_node_ch->box_max_y - ptr_node_ch->box_ts_y + 2; jj++)
            {
                for (int ii = ptr_node_ch->box_min_x - ptr_node_ch->box_ts_x; ii < ptr_node_ch->box_max_x - ptr_node_ch->box_ts_x + 2; ii++)
                {
                    box_grid_idx = ii + jj * (ptr_node_ch->box_real_dim_x + 1) + kk * (ptr_node_ch->box_real_dim_x + 1) * (ptr_node_ch->box_real_dim_y + 1);
    
                    box_idx = ii + jj * ptr_node_ch->box_real_dim_x + kk * ptr_node_ch->box_real_dim_x * ptr_node_ch->box_real_dim_y;
                    box_idxNbr_i0_j0_k0 = box_idx;
                    box_idxNbr_im1_j0_k0 = box_idx - 1;
                    box_idxNbr_i0_jm1_k0 = box_idx - ptr_node_ch->box_real_dim_x;
                    box_idxNbr_im1_jm1_k0 = box_idx - 1 - ptr_node_ch->box_real_dim_x;
                    box_idxNbr_i0_j0_km1 = box_idx - ptr_node_ch->box_real_dim_x * ptr_node_ch->box_real_dim_y;
                    box_idxNbr_im1_j0_km1 = box_idx - 1 - ptr_node_ch->box_real_dim_x * ptr_node_ch->box_real_dim_y;
                    box_idxNbr_i0_jm1_km1 = box_idx - ptr_node_ch->box_real_dim_x - ptr_node_ch->box_real_dim_x * ptr_node_ch->box_real_dim_y;
                    box_idxNbr_im1_jm1_km1 = box_idx - 1 - ptr_node_ch->box_real_dim_x - ptr_node_ch->box_real_dim_x * ptr_node_ch->box_real_dim_y;

                    //** >> The grid point exist **/
                    if (ptr_node_ch->ptr_box[box_idxNbr_i0_j0_k0] > -4 ||
                        ptr_node_ch->ptr_box[box_idxNbr_im1_j0_k0] > -4 ||
                        ptr_node_ch->ptr_box[box_idxNbr_i0_jm1_k0] > -4 ||
                        ptr_node_ch->ptr_box[box_idxNbr_im1_jm1_k0] > -4 ||
                        ptr_node_ch->ptr_box[box_idxNbr_i0_j0_km1] > -4 ||
                        ptr_node_ch->ptr_box[box_idxNbr_im1_j0_km1] > -4 ||
                        ptr_node_ch->ptr_box[box_idxNbr_i0_jm1_km1] > -4 ||
                        ptr_node_ch->ptr_box[box_idxNbr_im1_jm1_km1] > -4)
                    {
                        
                        is_bder_grid_point = false;
                        //** Connection to the right  **/
                        if (ptr_node_ch->ptr_box[box_idxNbr_i0_j0_k0] < -3 &&
                            ptr_node_ch->ptr_box[box_idxNbr_i0_jm1_k0] < -3 &&
                            ptr_node_ch->ptr_box[box_idxNbr_i0_j0_km1] < -3 &&
                            ptr_node_ch->ptr_box[box_idxNbr_i0_jm1_km1] < -3)
                        {
                            is_bder_grid_point = true;
                        }
                        //** Connection to the left  **/
                        else if (ptr_node_ch->ptr_box[box_idxNbr_im1_j0_k0] < -3 &&
                                 ptr_node_ch->ptr_box[box_idxNbr_im1_jm1_k0] < -3 &&
                                 ptr_node_ch->ptr_box[box_idxNbr_im1_j0_km1] < -3 &&
                                 ptr_node_ch->ptr_box[box_idxNbr_im1_jm1_km1] < -3)
                        {
                            is_bder_grid_point = true;
                        }
                        //** Backward connection   **/
                        else if (ptr_node_ch->ptr_box[box_idxNbr_i0_j0_k0] < -3 &&
                                ptr_node_ch->ptr_box[box_idxNbr_im1_j0_k0] < -3 &&
                                ptr_node_ch->ptr_box[box_idxNbr_i0_j0_km1] < -3 &&
                                ptr_node_ch->ptr_box[box_idxNbr_im1_j0_km1] < -3)
                        {
                            is_bder_grid_point = true;
                        }
                        //** Forward connection   **/
                        else if (ptr_node_ch->ptr_box[box_idxNbr_i0_jm1_k0] < -3 &&
                                ptr_node_ch->ptr_box[box_idxNbr_im1_jm1_k0] < -3 &&
                                ptr_node_ch->ptr_box[box_idxNbr_i0_jm1_km1] < -3 &&
                                ptr_node_ch->ptr_box[box_idxNbr_im1_jm1_km1] < -3)
                        {
                            is_bder_grid_point = true;
                        }
                        //** Upward connection **/
                        else if (ptr_node_ch->ptr_box[box_idxNbr_i0_j0_k0] < -3 &&
                                 ptr_node_ch->ptr_box[box_idxNbr_im1_j0_k0] < -3 &&
                                 ptr_node_ch->ptr_box[box_idxNbr_i0_jm1_k0] < -3 &&
                                 ptr_node_ch->ptr_box[box_idxNbr_im1_jm1_k0] < -3)
                        {
                            is_bder_grid_point = true;
                        }
                        //** Down connection **/
                        else if (ptr_node_ch->ptr_box[box_idxNbr_i0_j0_km1] < -3 &&
                                 ptr_node_ch->ptr_box[box_idxNbr_im1_j0_km1] < -3 &&
                                 ptr_node_ch->ptr_box[box_idxNbr_i0_jm1_km1] < -3 &&
                                 ptr_node_ch->ptr_box[box_idxNbr_im1_jm1_km1] < -3)
                        {
                            is_bder_grid_point = true;
                        }

                        //** >> Adding the grid point **/

                        //** >> Border grid point**/
                        if(is_bder_grid_point == true)
                        {
                            //** >> Space checking of border grid points array**/
                            if (space_check(&(ptr_node_ch->grid_bder_cap), ptr_node_ch->grid_bder_size + 1, "p1i1", &(ptr_node_ch->ptr_grid_bder)) == _FAILURE_)
                            {
                                printf("Error, in space_check function\n");
                                return _FAILURE_;
                            }

                            //** >> Adding the grid point to the border array **/
                            ptr_node_ch->ptr_grid_bder[ptr_node_ch->grid_bder_size] = box_grid_idx;
                            ptr_node_ch->grid_bder_size += 1; // Increasing the number of interior grid points in the array
                        }
                        //** Interior grid point **/
                        else
                        {
                            //** >> Space checking of border grid points array**/
                            if (space_check(&(ptr_node_ch->grid_intr_cap), ptr_node_ch->grid_intr_size + 1, "p1i1", &(ptr_node_ch->ptr_grid_intr)) == _FAILURE_)
                            {
                                printf("Error, in space_check function\n");
                                return _FAILURE_;
                            }

                            //** >> Adding the grid point to the interior array **/
                            ptr_node_ch->ptr_grid_intr[ptr_node_ch->grid_intr_size] = box_grid_idx;
                            ptr_node_ch->grid_intr_size += 1; // Increasing the number of interior grid points in the array
                        }
                    }
                }
            }
        }

        //** >> Refinement Criterion **/
        cap = ptr_node_ch->box_real_dim_x * ptr_node_ch->box_real_dim_y * ptr_node_ch->box_real_dim_z;
        ptr_node_ch->ptr_box_mass = (vtype *)calloc(cap, sizeof(vtype));
        //** >> Local mass **/
        for (int j = 0; j < ptr_node_ch->ptcl_size; j++)
        {
            ptr_node_ch->local_mass += GL_ptcl_mass[ptr_node_ch->ptr_ptcl[j]];
        }

        //* Potential, Acceleration and density of the grid **/
        cap = (ptr_node_ch->box_real_dim_x + 1) * (ptr_node_ch->box_real_dim_y + 1) * (ptr_node_ch->box_real_dim_z + 1);
        ptr_node_ch->ptr_pot = (vtype *)calloc(cap , sizeof(vtype));
        ptr_node_ch->ptr_ax = (vtype *)calloc(cap , sizeof(vtype));
        ptr_node_ch->ptr_ay = (vtype *)calloc(cap, sizeof(vtype));
        ptr_node_ch->ptr_az = (vtype *)calloc(cap, sizeof(vtype));
        ptr_node_ch->ptr_d = (vtype *)calloc(cap, sizeof(vtype));

        //** >> Tree structure **/
        ptr_node_pt->pptr_chn[i] = ptr_node_ch; // Parent node pointing to child node i
        ptr_node_ch->ptr_pt = ptr_node_pt;      // Child node i pointig to parent node

    } // End filling child nodes

    int sum;

    for (int i = 0; i < ptr_node_pt->zones_size; i++)
    {
        sum = 0;
        for (int j = 0; j < ptr_node_pt->ptr_zone_size[i]; j++)
        {
            box_idx_x = ptr_node_pt->ptr_cell_idx_x[ptr_node_pt->pptr_zones[i][j]] - ptr_node_pt->box_ts_x;
            box_idx_y = ptr_node_pt->ptr_cell_idx_y[ptr_node_pt->pptr_zones[i][j]] - ptr_node_pt->box_ts_y;
            box_idx_z = ptr_node_pt->ptr_cell_idx_z[ptr_node_pt->pptr_zones[i][j]] - ptr_node_pt->box_ts_z;
            box_idx = box_idx_x + box_idx_y * ptr_node_pt->box_real_dim_x + box_idx_z * ptr_node_pt->box_real_dim_x * ptr_node_pt->box_real_dim_y;
            sum += ptr_cell_ptcl_size[box_idx];
        }
    }

    return _SUCCESS_;
} // end function fill_child_nodes

static int fill_tentacles(const struct node *ptr_node_pt)
{
    int lv;
    int cap;
    int size;

    struct node **pptr_node_aux; // Auxiliary pointer in realloc

    lv = ptr_node_pt->lv - lmin + 1; // Children level

    cap = GL_tentacles_cap[lv];
    size = GL_tentacles_size[lv];
    // printf("lv = %d, ID pt = %d, size tenta = %d, cap tenta = %d\n", lv + lmin, ptr_node_pt->ID, size, cap);
    if (size + ptr_node_pt->chn_size > cap)
    {
        //** >> Increase the capacity **/
        cap = 2 * (size + ptr_node_pt->chn_size);

        if (GL_tentacles_old[lv] == NULL)
        {
            GL_tentacles_old[lv] = (struct node **)malloc(cap * sizeof(struct node *));
            GL_tentacles_new[lv] = (struct node **)malloc(cap * sizeof(struct node *));
        }
        else
        {
            //** >> Old tentacle **/
            pptr_node_aux = NULL;
            pptr_node_aux = (struct node **)realloc(GL_tentacles_old[lv], cap * sizeof(struct node *));
            if (pptr_node_aux == NULL)
            {
                printf("Error in reallocate the GL_tentacles_old[%d] array\n", lv);
                return _FAILURE_;
            }
            GL_tentacles_old[lv] = pptr_node_aux;

            //** >> New tentacle **/
            pptr_node_aux = NULL;
            pptr_node_aux = (struct node **)realloc(GL_tentacles_new[lv], cap * sizeof(struct node *));
            if (pptr_node_aux == NULL)
            {
                printf("Error in reallocate the GL_tentacles_new[%d] array\n", lv);
                return _FAILURE_;
            }
            GL_tentacles_new[lv] = pptr_node_aux;
        }

        //** >> Initializing new allocations **/
        for (int i = size; i < cap; i++)
        {
            GL_tentacles_old[lv][i] = NULL;
            GL_tentacles_new[lv][i] = NULL;
        }

        //** >> Memory computing **/
        TOTAL_MEMORY_TENTACLES += 2 * (2 * (size + ptr_node_pt->chn_size) - GL_tentacles_cap[lv]) * sizeof(struct node *);

        //** >> New value for the capacity of the tentacles at level lv**/
        GL_tentacles_cap[lv] = cap;
    }

    //** >> Filling the new tentacles for the next level of refinement **/
    for (int i = 0; i < ptr_node_pt->chn_size; i++)
    {
        //** >> Putting elements in the new tentacles **/
        GL_tentacles_old[lv][size + i] = ptr_node_pt->pptr_chn[i];
        GL_tentacles_new[lv][size + i] = ptr_node_pt->pptr_chn[i];
    }
    //** Increasing the number of structs in the level lv **/
    GL_tentacles_size[lv] += ptr_node_pt->chn_size;
    //** >> Increasing the deepth of levels of the tentacles **/
    if (GL_tentacles_level_max < lv)
    {
        GL_tentacles_level_max++;
    }

    return _SUCCESS_;
} // end function fill_tentacles

int tree_construction()
{

    //** >> Working in the refinement zones **/
    if (lmin < lmax)
    {
        struct node *ptr_node; // Parent node
        int **pptr_cell_ptcl;    // Particles indexes in each cell of the main node
        int *ptr_cell_ptcl_cap;  // Maximum number (capacity) of particles in each cell of the box
        int *ptr_cell_ptcl_size; // Number of particles in each cell of the box
        int space_cell_ptcl;     // Total space of the pptr_cell_ptcl array
        int required_space_cell_ptcl; // New total space required to the array pptr_cell_ptcl

        int no_pts; // Number of parents in the cycle

        int aux_cap; // Old capacity of the space cell particles array. Used in space checking

        int lv;

        //** >> Initializing cell particles pointers **/
        pptr_cell_ptcl = (int **)malloc(box_side_lmin_pow3 * sizeof(int *)); // The default value is equal to the number of cell in the coarsest level
        ptr_cell_ptcl_cap = (int *)calloc(box_side_lmin_pow3, sizeof(int));  // Initializing the capacity at 0 to say there is no allocation in that cell
        ptr_cell_ptcl_size = (int *)calloc(box_side_lmin_pow3, sizeof(int)); // By default there are 0 particles in the cells
        space_cell_ptcl = box_side_lmin_pow3;
        //** >> Initializing pointers of each cell **/
        for (int i = 0; i < box_side_lmin_pow3; i++)
        {
            pptr_cell_ptcl[i] = NULL;
        }

        // Comienza el ciclo for:
        //** >> For cycle over refinement levels **/
        //for (int lv = 0; lv < lmax - lmin + 1; lv++)
        lv = 0;
        while (lv < lmax - lmin + 1 && lv <= GL_tentacles_level_max)
        {
            no_pts = GL_tentacles_size[lv];
            //** >> For cycle over parent nodes **/
            for (int i = 0; i < no_pts; i++)
            {
                ptr_node = GL_tentacles_old[lv][i];

                required_space_cell_ptcl = ptr_node->box_real_dim_x * ptr_node->box_real_dim_y * ptr_node->box_real_dim_z;

                //** >> Space checking of ptcl_cell arrays**/

                if (space_cell_ptcl < required_space_cell_ptcl)
                {
                    aux_cap = space_cell_ptcl;
                    if (space_check(&(space_cell_ptcl), required_space_cell_ptcl, "p3i2i1i1", &(pptr_cell_ptcl), &(ptr_cell_ptcl_cap), &(ptr_cell_ptcl_size)) == _FAILURE_)
                    {
                        printf("Error, in space_check function\n");
                        return _FAILURE_;
                    }

                    //** >> Initiazling the new allocated values **/
                    for (int j = 0; j < space_cell_ptcl - aux_cap; j++)
                    {
                        pptr_cell_ptcl[space_cell_ptcl - j - 1] = NULL;
                        ptr_cell_ptcl_cap[space_cell_ptcl - j - 1] = 0;
                    }
                }

                for (int j = 0; j < space_cell_ptcl; j++)
                {
                    ptr_cell_ptcl_size[j] = 0;
                }

                //** >> Filling the cell particles arrays **/
                if (particle_distribution_in_cell(pptr_cell_ptcl, ptr_cell_ptcl_cap, ptr_cell_ptcl_size, ptr_node) == _FAILURE_)
                {
                    printf("Error at function particle_distribution_in_cell() \n");
                    return _FAILURE_;
                }

                //** >> Filling the mass box **/
                filling_mass_box(pptr_cell_ptcl, ptr_cell_ptcl_size, ptr_node);

                //** >> Computing the zones of refinement **/
                if (lv < lmax - lmin)
                {
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
                    if(ptr_node->zones_size > 0)
                    {
                        if (fill_child_nodes(pptr_cell_ptcl, ptr_cell_ptcl_size, ptr_node) == _FAILURE_)
                        {
                            printf("Error at function create_child_nodes()\n");
                            return _FAILURE_;
                        }
                    }

                    //** >> Filling Tentacles for the next cycle at level of refinement
                    if (ptr_node->chn_size > 0)
                    {
                        if (fill_tentacles(ptr_node) == _FAILURE_)
                        {
                            printf("Error at function fill_tentacles()\n");
                            return _FAILURE_;
                        }
                    }

                }// End computing new refinement zones and children
            } // End of cycle over parent nodes
            lv++;
        } // End of cycle over refinement levels


        //** >> Free memory **/
        for (int i = 0; i < space_cell_ptcl; i++)
        {
            if (pptr_cell_ptcl[i] != NULL)
            {
                free(pptr_cell_ptcl[i]);
            }
        }
        free(ptr_cell_ptcl_size);
        free(ptr_cell_ptcl_cap);
        free(pptr_cell_ptcl);
    } // Finalized tree structure

        return _SUCCESS_;
}