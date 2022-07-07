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

//** >> Local Functions
int input_galaxies_merger(void);
int input_plummer_model(void);
int input_particle_initialization(void);
static int input_code_units(void);

int input_galaxies_merger(void)
{

    printf("Input by galaxies merger\n");

    
    // char input_name[100] = "../source/initial_conditions/"
    //                        "galaxies_IC_low_ptcl.txt";

    // char input_name[100] = "../source/initial_conditions/"
    //                        "galaxies_IC_very_low_ptcl.txt";

    // char input_name[100] = "../source/initial_conditions/"
    //                        "galaxies_IC_high_ptcl.txt";

    // char input_name[100] = "../source/initial_conditions/"
    //                        "galaxy1_alone.txt";

    // char input_name[100] = "../source/initial_conditions/"
    //                        "galaxies_model2.txt";

    // char input_name[100] = "../source/initial_conditions/"
    //                        "galaxies_merger_v2_30k_particles.txt";

    // char input_name[100] = "../source/initial_conditions/"
    //                        "galaxies_merger_v2_300k_particles.txt";

    // char input_name[100] = "../source/initial_conditions/"
    //                        "galaxies_merger_v2_3000k_particles.txt";

    // char input_name[100] = "../source/initial_conditions/"
    //                        "galaxies_merger_v2_70k_particles_with_gas.txt";

    // char input_name[100] = "../source/initial_conditions/"
    //                        "galaxies_merger_v2_700k_particles_with_gas.txt";

    // char input_name[100] = "../source/initial_conditions/"
    //                        "spiral_galaxy.txt";

    // char input_name[100] = "../source/initial_conditions/"
    //                        "spiral_galaxy_stiff_universe.txt";



//** >> GALIC simulations
    // char input_name[100] = "../source/initial_conditions/Galic/"
    //                        "galaxy_alone_15k.txt";

    // char input_name[100] = "../source/initial_conditions/Galic/"
    //                        "galaxy_alone_30k.txt";

    // char input_name[100] = "../source/initial_conditions/Galic/"
    //                        "galaxies_merger_30k+15k.txt";

    // char input_name[100] = "../source/initial_conditions/Galic/"
    //                        "galaxies_merger_30k+15k_v2.txt";

    // char input_name[100] = "../source/initial_conditions/Galic/"
    //                        "galaxies_merger_30k+15k_v3.txt";

    // char input_name[100] = "../source/initial_conditions/Galic/"
    //                        "galaxies_merger_30k+15k_v4.txt";

    char input_name[100] = "../source/initial_conditions/Galic/"
                           "galaxies_merger_300k+150k.txt";

    FILE *file = NULL;
    file = fopen(input_name, "r");
    if (file == NULL)
    {
        printf("Error opening galaxies file!\n");
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

    //char aux[3] = "%lf";

    for (int i = 0; i < GL_no_ptcl; i++)
    {
        // GL_ptcl_mass[i] = 100; // Total mass over number of particles
        if (fscanf(file, aux, &GL_ptcl_x[i]) == 0)
        {
            printf("\nERROR: Fail to read the position x at galaxies file\n");
        }
        if (fscanf(file, aux, &GL_ptcl_y[i]) == 0)
        {
            printf("\nERROR: Fail to read the position y at galaxies file\n");
        }
        if (fscanf(file, aux, &GL_ptcl_z[i]) == 0)
        {
            printf("\nERROR: Fail to read the position z at galaxies file\n");
        }
        if (fscanf(file, aux, &GL_ptcl_vx[i]) == 0)
        {
            printf("\nERROR: Fail to read the velocity x at galaxies file\n");
        }
        if (fscanf(file, aux, &GL_ptcl_vy[i]) == 0)
        {
            printf("\nERROR: Fail to read the velocity y at galaxies file\n");
        }
        if (fscanf(file, aux, &GL_ptcl_vz[i]) == 0)
        {
            printf("\nERROR: Fail to read the velocity z at galaxies file\n");
        }
        if (fscanf(file, aux, &GL_ptcl_mass[i]) == 0)
        {
            printf("\nERROR: Fail to read the mass at galaxies file\n");
        }
    }
    fclose(file);


    return _SUCCESS_;
}

int input_plummer_model(void)
{

    printf("Input by plummer\n");
    // char input_name[100] = "../source/initial_conditions/"
    //                        "plummer_center_moved.txt";

    // char input_name[100] = "../source/initial_conditions//"
    //                        "Plummer_Np(10000)_BoxSize(0.1)_Mass(1000000).txt";

    char input_name[100] = "../source/initial_conditions/"
                           "Plummer_v2.txt";

    FILE *file = NULL;
    file = fopen(input_name, "r");
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

int input_particle_initialization(void)
{
    printf("Input by random\n");

    //srand(88);
    srand(7878);

    vtype maxdistance = _User_BoxSize_ * 0.5;
    for (int i = 0; i < GL_no_ptcl; i++)
    {
        GL_ptcl_mass[i] = meanmass; // at solar mass;
        GL_ptcl_x[i] = ((vtype)rand() / RAND_MAX - 0.5L) * maxdistance;
        GL_ptcl_y[i] = ((vtype)rand() / RAND_MAX - 0.5L) * maxdistance;
        GL_ptcl_z[i] = ((vtype)rand() / RAND_MAX - 0.5L) * maxdistance;

        GL_ptcl_vx[i] = 0.0;
        GL_ptcl_vy[i] = 0.0;
        GL_ptcl_vz[i] = 0.0;

        GL_ptcl_ax[i] = 0;
        GL_ptcl_ay[i] = 0;
        GL_ptcl_az[i] = 0;
    }
    
    return _SUCCESS_;
}

static int input_code_units(void)
{
    for (int i = 0; i < GL_no_ptcl; i++)
    {
        GL_ptcl_x[i] = GL_ptcl_x[i] / _User_BoxSize_ + 0.5L;
        GL_ptcl_y[i] = GL_ptcl_y[i] / _User_BoxSize_ + 0.5L;
        GL_ptcl_z[i] = GL_ptcl_z[i] / _User_BoxSize_ + 0.5L;
        
        GL_ptcl_vx[i] = GL_ptcl_vx[i] * sqrt(_User_BoxSize_) ;
        GL_ptcl_vy[i] = GL_ptcl_vy[i] * sqrt(_User_BoxSize_) ;
        GL_ptcl_vz[i] = GL_ptcl_vz[i] * sqrt(_User_BoxSize_) ;

        GL_ptcl_ax[i] = GL_ptcl_ax[i] * _User_BoxSize_ * _User_BoxSize_;
        GL_ptcl_ay[i] = GL_ptcl_ay[i] * _User_BoxSize_ * _User_BoxSize_;
        GL_ptcl_az[i] = GL_ptcl_az[i] * _User_BoxSize_ * _User_BoxSize_;
    }
    return _SUCCESS_;
}

int input(void)
{
    // if (input_plummer_model() == _FAILURE_)
    // {
    //     printf("\n\n Error running input_plummer_model function\n\n");
    //     return _FAILURE_;
    // }

    if (input_galaxies_merger() == _FAILURE_)
    {
        printf("\n\n Error running input_galaxies_merger function\n\n");
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
