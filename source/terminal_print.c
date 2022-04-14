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

void terminal_print()
{
    double sum = TOTAL_MEMORY_NODES + TOTAL_MEMORY_CELDAS + TOTAL_MEMORY_PARTICULAS + TOTAL_MEMORY_CAJAS + TOTAL_MEMORY_AUX + TOTAL_MEMORY_TENTACLES + TOTAL_MEMORY_OTROS;
    int max_memory_secction = 0;

    double Total_memory[7] = {TOTAL_MEMORY_NODES,
                              TOTAL_MEMORY_CELDAS,
                              TOTAL_MEMORY_PARTICULAS,
                              TOTAL_MEMORY_CAJAS,
                              TOTAL_MEMORY_AUX,
                              TOTAL_MEMORY_TENTACLES,
                              TOTAL_MEMORY_OTROS};

    for (int i = 0; i < 7; i++)
    {
        max_memory_secction = Total_memory[max_memory_secction] > Total_memory[i] ? max_memory_secction : i;
    }

    printf("\n\n%sMEMORY INITIAL TREE [MB]:%s\n\n", KRED, KNRM);

    char Memory_names[50][100] = {
        "Nodes",
        "Celdas",
        "Particulas",
        "Cajas",
        "Auxiliary",
        "Tentacles",
        "OTROS",
        ""};

    for (int i = 0; i < 7; i++)
    {
        if (i == max_memory_secction)
        {
            printf("%s%s = %f%s\n", KCYN, Memory_names[i], Total_memory[i] / 1000000, KNRM);
        }
        else
        {
            printf("%s = %1.3e\n", Memory_names[i], Total_memory[i]);
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
            printf("%s%s = %1.3e%s\n", KCYN, Time_names[i], GL_times[i], KNRM);
        }
        else
        {
            printf("%s = %1.3e\n", Time_names[i], GL_times[i]);
        }
    }
    //** Potential times **/

    if (MAX_TIME_SECCTION == 5)
    {
        for (int i = 20; i < 26; i++)
        {
            if (MAX_TIME_SECCTION_POTENTIAL == i)
            {
                printf("%s%s = %1.3e%s\n", KYEL, Time_names_potential[i-20], GL_times[i], KNRM);
            }
            else
            {
                printf("%s = %1.3e\n", Time_names_potential[i-20], GL_times[i]);
            }
        }
    }
    else
    {
        for (int i = 20; i < 26; i++)
        {
            printf("%s = %1.3e\n", Time_names_potential[i-20], GL_times[i]);
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
            printf("%s%s = %1.3e%s\n", KCYN, Time_names[i], GL_times[i], KNRM);
        }
        else
        {
            printf("%s = %1.3e\n", Time_names[i], GL_times[i]);
        }
    }

    //** >> Output Time **/
    if (18 == MAX_TIME_SECCTION)
    {
        printf("%sOutput Main Parameters  = %1.3e%s\n", KCYN, GL_times[18], KNRM);
    }
    else
    {
        printf("Output Main Parameters = %.3e\n", GL_times[18]);
    }

    //** >> Output Time **/
    if (19 == MAX_TIME_SECCTION)
    {
        printf("%sOutput Snapshots = %1.3e%s\n", KCYN, GL_times[19], KNRM);
    }
    else
    {
        printf("Output Snapshots = %.3e\n", GL_times[19]);
    }

    //** >> Total time **/
    printf("\n%sTOTAL = %.3e s%s\n", KMAG, TOTAL_TIME, KNRM);
}