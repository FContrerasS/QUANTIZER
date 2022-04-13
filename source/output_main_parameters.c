/*
 * output_main_parameters.c
 *
 * Export a .txt file with the main parameters used in the simulation
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

#include "output_main_parameters.h"

void output_main_parameters(vtype final_time, int Number_timesteps, int Number_outputs)
{
    char parameters_name[1100];
    sprintf(parameters_name, "%s/Parameters.dat", folder_name);

    FILE *file = fopen(parameters_name, "w");
    if (file == NULL)
    {
        printf("Error opening file!");
        exit(EXIT_FAILURE);
    }

    fprintf(file, "lmin %d\n", lmin);
    fprintf(file, "lmax %d\n", lmax);
    fprintf(file, "_MAX_NUMBER_OF_ITERATIONS_IN_POISSON_EQUATION_ %d\n", _MAX_NUMBER_OF_ITERATIONS_IN_POISSON_EQUATION_);
    fprintf(file, "_ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_ %1.5e\n", _ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_);
    fprintf(file, "check_poisson_error_method %d\n", check_poisson_error_method);
    fprintf(file, "multigrid_method %d\n", multigrid);
    fprintf(file, "solver %d\n", solver);
    fprintf(file, "_NiterPreS_ %d\n", _NiterPreS_);
    fprintf(file, "_NiterPostS_ %d\n", _NiterPostS_);
    fprintf(file, "_Niterfinddphic_ %d\n", _Niterfinddphic_);
    fprintf(file, "_Iter_branches_SOR_ %d\n", _Iter_branches_SOR_);
    fprintf(file, "_w_SOR_ %1.5e\n", _w_SOR_);
    fprintf(file, "_CFL_ %1.5e\n", _CFL_);
    fprintf(file, "_MAX_dt_ %1.5e\n", _MAX_dt_);

    fprintf(file, "Ngrid %d\n", no_lmin_cell + 1);
    fprintf(file, "Nparticles %d\n", GL_no_ptcl);
    fprintf(file, "Number_timesteps	%d\n", Number_timesteps);
    fprintf(file, "Number_outputs %d\n", Number_outputs);
    fprintf(file, "fr_output %d\n", fr_output);
    fprintf(file, "BoxSize %1.5e\n", (double)_User_BoxSize_);
    fprintf(file, "vtype %d\n", _VTYPE_);
    fprintf(file, "Maxdt %1.5e\n", (double)(final_time / _Mgyear_)); // In Mega years
    fprintf(file, "G %1.5e\n", (double)_G_);

    fclose(file);
}