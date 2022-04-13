/*
 * input.c
 *
 * Input of the simulation. Mass, position and velocitiy
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

#include "input.h"

int input_plummer_model()
{

    printf("Input by plummer\n");
    // char input_name[100] = "../source/initial_conditions/"
    //                              "ultima_esperanza.txt"; 

    char input_name[100] = "../source/initial_conditions//"
                           "Plummer_Np(10000)_BoxSize(0.1)_Mass(1000000).txt";

    FILE *file = fopen(input_name, "r");
    if (file == NULL)
    {
        printf("Error opening plummer file!\n");
        exit(EXIT_FAILURE);
    }

#if _VTYPE_ == 1
    char aux[3] = "%f";
#elif _VTYPE_ == 2
    char aux[4] = "%lf";
#elif _VTYPE_ == 3
    char aux[4] = "%Lf";
#else
    char aux[4] = "%lf";
#endif

    for (int i = 0; i < GL_no_ptcl; i++)
    {
        GL_ptcl_mass[i] = 100; // Total mass over number of particles
        if (fscanf(file, aux, &GL_ptcl_x[i]) == 0)
        {
            printf("\nERROR: Fail to read the position x at plummer file\n");
        }
        if (fscanf(file, aux, &GL_ptcl_y[i]) == 0)
        {
            printf("\nERROR: Fail to read the position y at plummer file\n");
        }
        if (fscanf(file, aux, &GL_ptcl_z[i]) == 0)
        {
            printf("\nERROR: Fail to read the position z at plummer file\n");
        }
        if (fscanf(file, aux, &GL_ptcl_vx[i]) == 0)
        {
            printf("\nERROR: Fail to read the velocity x at plummer file\n");
        }
        if (fscanf(file, aux, &GL_ptcl_vy[i]) == 0)
        {
            printf("\nERROR: Fail to read the velocity y at plummer file\n");
        }
        if (fscanf(file, aux, &GL_ptcl_vz[i]) == 0)
        {
            printf("\nERROR: Fail to read the velocity z at plummer file\n");
        }
    }
    fclose(file);

    return _SUCCESS_;
}

int input_particle_initialization()
{
    printf("Input by random\n");

    //srand(88);
    srand(7878);

    vtype maxdistance = _User_BoxSize_ * 0.5;
    for (int i = 0; i < GL_no_ptcl; i++)
    {
        GL_ptcl_mass[i] = meanmass; // at solar mass;
        GL_ptcl_x[i] = ((vtype)rand() / RAND_MAX - 0.5) * maxdistance;
        GL_ptcl_y[i] = ((vtype)rand() / RAND_MAX - 0.5) * maxdistance;
        GL_ptcl_z[i] = ((vtype)rand() / RAND_MAX - 0.5) * maxdistance;

        GL_ptcl_vx[i] = 0.0;
        GL_ptcl_vy[i] = 0.0;
        GL_ptcl_vz[i] = 0.0;

        GL_ptcl_ax[i] = 0;
        GL_ptcl_ay[i] = 0;
        GL_ptcl_az[i] = 0;
    }
    
    return _SUCCESS_;
}

static int input_code_units()
{
    for (int i = 0; i < GL_no_ptcl; i++)
    {
        GL_ptcl_x[i] = GL_ptcl_x[i] / _User_BoxSize_ + 0.5;
        GL_ptcl_y[i] = GL_ptcl_y[i] / _User_BoxSize_ + 0.5;
        GL_ptcl_z[i] = GL_ptcl_z[i] / _User_BoxSize_ + 0.5;
        
        GL_ptcl_vx[i] = GL_ptcl_vx[i] * sqrt(_User_BoxSize_);
        GL_ptcl_vy[i] = GL_ptcl_vy[i] * sqrt(_User_BoxSize_);
        GL_ptcl_vz[i] = GL_ptcl_vz[i] * sqrt(_User_BoxSize_);

        GL_ptcl_ax[i] = GL_ptcl_ax[i] * _User_BoxSize_ * _User_BoxSize_;
        GL_ptcl_ay[i] = GL_ptcl_ay[i] * _User_BoxSize_ * _User_BoxSize_;
        GL_ptcl_az[i] = GL_ptcl_az[i] * _User_BoxSize_ * _User_BoxSize_;
    }
    return _SUCCESS_;
}

    int
    input()
{

    if (input_plummer_model() == _FAILURE_)
    {
        printf("\n\n Error running input_plummer_model function\n\n");
        return _FAILURE_;
    }

/*     if (input_particle_initialization() == _FAILURE_)
    {
        printf("\n\n Error running input_particle_initialization function\n\n");
        return _FAILURE_;
    } */

    if (input_code_units() == _FAILURE_)
    {
        printf("\n\n Error running input_codedimensions function\n\n");
        return _FAILURE_;
    }

    return _SUCCESS_;
}
