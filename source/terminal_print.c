/*
 * terminal_print.c
 *
 * Print variables in terminal
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

#include "terminal_print.h"

static void computing_memory()
{

    int no_pts;

    struct node *ptr_node;

    //** >> Tentacles **/
    TOTAL_MEMORY_TENTACLES += (lmax - lmin + 1) * sizeof(struct node **);

    for (int lv = 0; lv < lmax - lmin + 1; lv++)
    {
        TOTAL_MEMORY_TENTACLES += GL_tentacles_cap[lv] * sizeof(struct node *);
    }

    //** >> Global particles **/
    TOTAL_MEMORY_PARTICLES += GL_no_ptcl * (sizeof(bool) + 10 * sizeof(vtype)); // Global particles

    //** >> Node properties **/
    for (int lv = 0; lv < GL_tentacles_level_max + 1; lv++)
    {
        no_pts = GL_tentacles_size[lv];

        for (int i = 0; i < no_pts; i++)
        {
            ptr_node = GL_tentacles[lv][i];
            TOTAL_MEMORY_CELDAS += 3 * ptr_node->cell_cap * sizeof(int);
            
            if(lmin < lmax)
            {
                TOTAL_MEMORY_CELL_STRUCT += 2* ptr_node->box_cap * sizeof(struct cell_struct);
                if (ptr_node->ptr_cell_struct != NULL)
                {
                    for (int j = 0; j < ptr_node->box_cap; j++)
                    {
                        TOTAL_MEMORY_CELL_STRUCT += ptr_node->ptr_cell_struct[j].ptcl_cap * sizeof(int);
                    }
                }

                if (ptr_node->ptr_cell_struct_old != NULL)
                {
                    for (int j = 0; j < ptr_node->box_cap; j++)
                    {
                        TOTAL_MEMORY_CELL_STRUCT += ptr_node->ptr_cell_struct_old[j].ptcl_cap * sizeof(int);
                    }     
                }

            }
            TOTAL_MEMORY_CAJAS += 2 * ptr_node->box_cap * (sizeof(int) + sizeof(vtype)); //Boxes and mass boxes
            TOTAL_MEMORY_GRID_POINTS += 2 * (ptr_node->grid_bder_cap + ptr_node->grid_intr_cap) * sizeof(int); // Grid interior and border points
            TOTAL_MEMORY_GRID_PROPERTIES += 5 * ptr_node->grid_properties_cap * sizeof(vtype); // Grid properties, accelerations, potential and density
            TOTAL_MEMORY_AUX += ptr_node->zones_cap * sizeof(int *) + ptr_node->cell_ref_cap * sizeof(int);
            for (int j = 0; j < ptr_node->zones_cap;j++)
            {
                TOTAL_MEMORY_AUX += ptr_node->ptr_zone_cap[j] * sizeof(int);
            }
            TOTAL_MEMORY_AUX += ptr_node->aux_idx_cap * sizeof(int);
            
        }
        TOTAL_MEMORY_NODES += no_pts * sizeof(struct node);
    }

    //Stack of memory pool
    ptr_node = GL_pool_node_start;
    while(ptr_node != NULL)
    {
        TOTAL_MEMORY_STACK += 3 * ptr_node->cell_cap * sizeof(int);
        TOTAL_MEMORY_STACK += ptr_node->box_cap * sizeof(struct cell_struct);
        if(lmin < lmax)
        {
            TOTAL_MEMORY_STACK += 2 * ptr_node->box_cap * sizeof(struct cell_struct);
            if (ptr_node->ptr_cell_struct != NULL)
            {
                for (int j = 0; j < ptr_node->box_cap; j++)
                {
                    TOTAL_MEMORY_STACK += ptr_node->ptr_cell_struct[j].ptcl_cap * sizeof(int);
                }
            }

            if (ptr_node->ptr_cell_struct_old != NULL)
            {
                for (int j = 0; j < ptr_node->box_cap; j++)
                {
                    TOTAL_MEMORY_STACK += ptr_node->ptr_cell_struct_old[j].ptcl_cap * sizeof(int);
                }
            }
        }
        TOTAL_MEMORY_STACK += 2 * ptr_node->box_cap * (sizeof(int) + sizeof(vtype));                       // Boxes and mass boxes
        TOTAL_MEMORY_STACK += 2 * (ptr_node->grid_bder_cap + ptr_node->grid_intr_cap) * sizeof(int);       // Grid interior and border points
        TOTAL_MEMORY_STACK += 5 * ptr_node->grid_properties_cap * sizeof(vtype);                           // Grid properties, accelerations, potential and density
        TOTAL_MEMORY_STACK += ptr_node->zones_cap * sizeof(int *) + ptr_node->cell_ref_cap * sizeof(int);
        for (int j = 0; j < ptr_node->zones_cap; j++)
        {
            TOTAL_MEMORY_STACK += ptr_node->ptr_zone_cap[j] * sizeof(int);
        }
        TOTAL_MEMORY_STACK += ptr_node->aux_idx_cap * sizeof(int);
        TOTAL_MEMORY_STACK += sizeof(struct node);
        if (ptr_node == GL_pool_node_end)
        {
            ptr_node = NULL;
        }
        else
        {
            ptr_node = ptr_node->ptr_pt;
        }
    }
    
}

void terminal_print()
{

    computing_memory();

    double sum = TOTAL_MEMORY_NODES + TOTAL_MEMORY_CELDAS + TOTAL_MEMORY_PARTICLES + TOTAL_MEMORY_CELL_STRUCT + TOTAL_MEMORY_CAJAS + TOTAL_MEMORY_GRID_POINTS + TOTAL_MEMORY_GRID_PROPERTIES + TOTAL_MEMORY_AUX + TOTAL_MEMORY_TENTACLES + TOTAL_MEMORY_STACK;
    int max_memory_secction = 0;

    double Total_memory[10] = {TOTAL_MEMORY_NODES,
                              TOTAL_MEMORY_CELDAS,
                              TOTAL_MEMORY_PARTICLES,
                              TOTAL_MEMORY_CELL_STRUCT,
                              TOTAL_MEMORY_CAJAS,
                              TOTAL_MEMORY_GRID_POINTS,
                              TOTAL_MEMORY_GRID_PROPERTIES,
                              TOTAL_MEMORY_AUX,
                              TOTAL_MEMORY_TENTACLES,
                              TOTAL_MEMORY_STACK};

    for (int i = 0; i < 10; i++)
    {
        max_memory_secction = Total_memory[max_memory_secction] > Total_memory[i] ? max_memory_secction : i;
    }

    printf("\n\n%sMEMORY [MB]:%s\n\n", KRED, KNRM);

    char Memory_names[50][100] = {
        "Nodes",
        "Celdas",
        "Particulas",
        "Cell struct",
        "Cajas",
        "Grid Points",
        "Grid Properties",
        "Auxiliary",
        "Tentacles",
        "Memory Pool",
        };

    for (int i = 0; i < 10; i++)
    {
        if (i == max_memory_secction)
        {
            printf("%s%s = %f ~ %.1f %% %s\n", KCYN, Memory_names[i], Total_memory[i] / 1000000,Total_memory[i] * 100 /sum  , KNRM);
        }
        else
        {
            printf("%s = %f ~ %.1f %%\n", Memory_names[i], Total_memory[i] / 1000000,Total_memory[i] * 100 /sum);
        }
    }

    printf("\n%sTOTAL = %f MB%s\n\n", KMAG, sum / 1000000, KNRM);

    printf("\n\n%sTIEMPOS [s]:%s\n\n", KRED, KNRM);

    char Time_names[50][100] = {
        "Global variables",
        "Input",
        "Initialization",
        "Tree Construction",
        "Grid Density",
        "Total potential",
        "Grid Acceleration",
        "Particle acceleration",
        "Time-step computing",
        "Particle Updating A",
        "Tree Adaptation",
        "Reset",
        "Particle Updating B",
        "Observables",
        "Output Main Parameters",
        "Output Snapshots",
        ""};

    char Time_names_potential[50][100] = {
        "Head Potential",
        "Error Check in the Head",
        "Potential transfer in branches",
        "Filling Red and Black in branches",
        "Potential in branches",
        "Error Check in branches",
        ""};

    double TOTAL_TIME = 0;
    int MAX_TIME_SECCTION = 0;
    for (int i = 0; i < 20; i++)
    {
        TOTAL_TIME += GL_times[i];

        MAX_TIME_SECCTION = GL_times[MAX_TIME_SECCTION] > GL_times[i] ? MAX_TIME_SECCTION : i;
    }

    int MAX_TIME_SECCTION_POTENTIAL = 20;
    for (int i = 20; i < 26; i++)
    {
        MAX_TIME_SECCTION_POTENTIAL = GL_times[MAX_TIME_SECCTION_POTENTIAL] > GL_times[i] ? MAX_TIME_SECCTION_POTENTIAL : i;
    }

    //** >> Input time **/
    //** >> Initialization Time
    //** >> Tree Construction time **/
    //** >> Grid Density time **/
    for (int i = 0; i < 5; i++)
    {
        if (i == MAX_TIME_SECCTION)
        {
            printf("%s%s = %1.2e ~ %.1f %% %s\n", KCYN, Time_names[i], GL_times[i], GL_times[i] * 100 / TOTAL_TIME, KNRM);
        }
        else
        {
            printf("%s = %1.2e ~ %.1f %%\n", Time_names[i], GL_times[i], GL_times[i] * 100 / TOTAL_TIME);
        }
    }
    //** Potential times **/

    if (MAX_TIME_SECCTION == 5)
    {
        for (int i = 20; i < 26; i++)
        {
            if (MAX_TIME_SECCTION_POTENTIAL == i)
            {
                printf("%s%s = %1.2e ~ %.1f %% %s\n", KYEL, Time_names_potential[i-20], GL_times[i], GL_times[i] * 100 / TOTAL_TIME, KNRM);
            }
            else
            {
                printf("%s = %1.2e ~ %.1f %%\n", Time_names_potential[i-20], GL_times[i] , GL_times[i] * 100 / TOTAL_TIME);
            }
        }
    }
    else
    {
        for (int i = 20; i < 26; i++)
        {
            printf("%s = %1.2e ~ %.1f %%\n", Time_names_potential[i-20], GL_times[i], GL_times[i] * 100 / TOTAL_TIME);
        }
    }



    //** Total potential time **/
    //** >> Grid Acceleration time **/
    //** >> Particle acceleration time **/
    //** >> Time-step time **/
    //** >> Particle Updating time **/
    for (int i = 5; i < 14; i++)
    {
        if (i == MAX_TIME_SECCTION)
        {
            printf("%s%s = %1.2e ~ %.1f %% %s\n", KCYN, Time_names[i], GL_times[i], GL_times[i] * 100 / TOTAL_TIME ,KNRM);
        }
        else
        {
            printf("%s = %1.2e ~ %.1f %%\n", Time_names[i], GL_times[i], GL_times[i] * 100 / TOTAL_TIME );
        }
    }

    //** >> Output Time **/
    if (18 == MAX_TIME_SECCTION)
    {
        printf("%sOutput Main Parameters  = %1.2e ~ %.1f %% %s\n", KCYN, GL_times[18],GL_times[18] * 100 / TOTAL_TIME, KNRM);
    }
    else
    {
        printf("Output Main Parameters = %.2e ~ %.1f %%\n", GL_times[18],GL_times[18] * 100 / TOTAL_TIME);
    }

    //** >> Output Time **/
    if (19 == MAX_TIME_SECCTION)
    {
        printf("%sOutput Snapshots = %.2e ~ %.1f %% %s\n", KCYN, GL_times[19],GL_times[19] * 100 / TOTAL_TIME, KNRM);
    }
    else
    {
        printf("Output Snapshots = %.2e ~ %.1f %%\n", GL_times[19],GL_times[19] * 100 / TOTAL_TIME);
    }

    //** >> Total time **/
    printf("\n%sTOTAL = %.3e s%s\n", KMAG, TOTAL_TIME, KNRM);

    char Time_names_tree_adaptation[50][100] = {
        "updating_cell_struct",
        "initialization_node_boxes",
        "initialization_ref_aux",
        "fill_cell_ref",
        "fill_zones_ref",
        "create_links",
        "create_links_2",
        "remov_cells_nolonger_require_refinement",
        "adapt_child_box_and_cells",
        "create_new_child_nodes",
        "moving_old_child_to_new_child",
        "moving_new_zones_to_new_child",
        "update_border_child_boxes",
        "reorganization_child_node",
        "reorganization_grandchild_node",
        "moved_unused_child_node_to_memory_pool",
        "updating_ref_zones_grandchildren",
        "update_child_grid_points",
        "exchange_box_aux_to_box",
        "tentacles_updating",
        "update_chn_size",
        "updating_tentacles_max_lv",
        "Creation_and_free_links_pointers",
        ""};

    printf("\n\n");
    //** >> TREE ADAPTATION TIME **/
    printf("\n\nTREE ADAPTATION TIME [s], Percentage over tree adaptation time\n\n");
    for (int i = 30; i < 53; i++)
    {
        printf("%d: %s = %1.2e ~ %.1f %%\n",i-30, Time_names_tree_adaptation[i - 30], GL_times[i],GL_times[i] * 100 / GL_times[10]);
    }
    printf("\n\n");
}
