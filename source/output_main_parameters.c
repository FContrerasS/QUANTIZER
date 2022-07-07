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

    FILE *file = NULL;
    file = fopen(parameters_name, "w");
    if (file == NULL)
    {
        printf("Error opening file!");
        exit(EXIT_FAILURE);
    }

    fprintf(file, "lmin %d\n", (int) lmin);
    fprintf(file, "lmax %d\n", (int) lmax);
    fprintf(file, "_MAX_NUMBER_OF_ITERATIONS_IN_POISSON_EQUATION_ %d\n", (int) _MAX_NUMBER_OF_ITERATIONS_IN_POISSON_EQUATION_);
    fprintf(file, "_ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_ %1.5e\n", (double) _ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_);
    fprintf(file, "_ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_2 %1.5e\n", (double)_ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_2);
    fprintf(file, "check_poisson_error_method %d\n", (int) check_poisson_error_method);
    fprintf(file, "multigrid_cycle %d\n", (int )multigrid_cycle);
    fprintf(file, "solverPreS %d\n", (int )solverPreS);
    fprintf(file, "solverfinddphic %d\n", (int) solverfinddphic);
    fprintf(file, "solverPostS %d\n",(int) solverPostS);
    fprintf(file, "_NiterPreS_ %d\n",(int) _NiterPreS_);
    fprintf(file, "_NiterPostS_ %d\n", (int) _NiterPostS_);
    fprintf(file, "_Niterfinddphic_ %d\n", _Niterfinddphic_);
    fprintf(file, "_Iter_branches_solver_ %d\n", (int) _Iter_branches_solver_);
    fprintf(file, "_w_SOR_ %1.5e\n", (double) _w_SOR_);
    fprintf(file, "_w_SOR_HEAD_ %1.5e\n", (double)_w_SOR_HEAD_);
    fprintf(file, "head_pot_method %d\n", (int)head_pot_method);
    fprintf(file, "_CFL_ %1.5e\n", (double)_CFL_);
    fprintf(file, "_MAX_dt_ %1.5e\n",(double) _MAX_dt_);
    fprintf(file, "ref_criterion_ptcl %d\n", (int) ref_criterion_ptcl);
    fprintf(file, "n_exp %d\n", (int) n_exp);
    fprintf(file, "force_stencil %d\n", (int)force_stencil);
    fprintf(file, "potential_energy_type %d\n", (int)potential_energy_type);
    fprintf(file, "GC_iter %d\n",(int) Garbage_Collector_iter);

    fprintf(file, "Ngrid %d\n",(int) (no_lmin_cell + 1));
    fprintf(file, "Nparticles %d\n", (int) GL_no_ptcl);
    fprintf(file, "Number_timesteps	%d\n", (int) Number_timesteps);
    fprintf(file, "Number_outputs %d\n", (int) Number_outputs);
    fprintf(file, "fr_output %d\n", (int) fr_output);
    fprintf(file, "BoxSize %1.5e\n", (double)_User_BoxSize_);
    fprintf(file, "vtype %d\n", _VTYPE_);
    fprintf(file, "Maxdt %1.5e\n", (double)(final_time / _Mgyear_)); // In Mega years
    // fprintf(file, "Length_unit kpc\n");
    // fprintf(file, "Mass_unit Solar_mass\n");
    fprintf(file, "Time_unit %f\n", (double)tt);
    fprintf(file, "G %1.5e\n", (double)_G_);



    fclose(file);
}