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

//* >> Local Functions
int input_galaxies_merger(void);
int input_plummer_model(void);
int input_particle_initialization(void);
void computing_cm_and_initial_mass(void);
static int input_code_units(void);

int input_galaxies_merger(void)
{

  printf("Input by galaxies merger\n");

  // char input_name[100] = "../input/Galaxy_merger_M1/"
  //                        "Galaxy_alone_G1_translated.csv";



  // char input_name[100] = "../input/Galaxy_merger_M1/"
  //                        "Galaxy_alone_G2_translated.csv";



  char input_name[100] = "../input/Galaxy_merger_M1/"
                         "Galaxy_merger_model_M1.csv";



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

  for (int i = 0; i < GL_no_ptcl_initial; i++)
  {
    if (fscanf(file, aux, &GL_ptcl_x[i]) == 0)
    {
      printf("\nERROR: Fail to read the position x at galaxies file\n");
      return _FAILURE_;
    }
    if (fscanf(file, aux, &GL_ptcl_y[i]) == 0)
    {
      printf("\nERROR: Fail to read the position y at galaxies file\n");
      return _FAILURE_;
    }
    if (fscanf(file, aux, &GL_ptcl_z[i]) == 0)
    {
      printf("\nERROR: Fail to read the position z at galaxies file\n");
      return _FAILURE_;
    }
    if (fscanf(file, aux, &GL_ptcl_vx[i]) == 0)
    {
      printf("\nERROR: Fail to read the velocity x at galaxies file\n");
      return _FAILURE_;
    }
    if (fscanf(file, aux, &GL_ptcl_vy[i]) == 0)
    {
      printf("\nERROR: Fail to read the velocity y at galaxies file\n");
      return _FAILURE_;
    }
    if (fscanf(file, aux, &GL_ptcl_vz[i]) == 0)
    {
      printf("\nERROR: Fail to read the velocity z at galaxies file\n");
      return _FAILURE_;
    }
    if (fscanf(file, aux, &GL_ptcl_mass[i]) == 0)
    {
      printf("\nERROR: Fail to read the mass at galaxies file\n");
      return _FAILURE_;
    }
  }


  return _SUCCESS_;
}

int input_plummer_model(void)
{

  printf("Input by plummer\n");

  // Plummer model using a scale parameter a = 1 kpc
  // char input_name[100] = "../input/"
  //                        "Plumer_Model/Plummer_Model_np_(1000)_a_(1000.0 pc)_M_(1.0e+08 Msun).csv";

    // char input_name[100] = "../input/"
    //                      "Plumer_Model/Plummer_Model_np_(10000)_a_(1000.0 pc)_M_(1.0e+08 Msun).csv";

    char input_name[100] = "../input/"
                         "Plumer_Model/Plummer_Model_np_(100000)_a_(1000.0 pc)_M_(1.0e+08 Msun).csv";                         



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

  for (int i = 0; i < GL_no_ptcl_initial; i++)
  {
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
    if (fscanf(file, aux, &GL_ptcl_mass[i]) == 0)
    {
      printf("\nERROR: Fail to read the mass at plummer file\n");
    }
  }
  fclose(file);

  return _SUCCESS_;
}

int input_particle_initialization(void)
{
  vtype sum = 0.0;
  printf("Input by random\n");

  // srand(88);
  srand(7878);


  //CIRCULAR ORBIT
  //Center particle
  GL_ptcl_mass[0] = 1.0; // at solar mass;
  sum += GL_ptcl_mass[0]; 
  GL_ptcl_x[0] = 0.0;// + _User_BoxSize_ / (1 << 7);
  GL_ptcl_y[0] = 0.0; // + _User_BoxSize_ / (1 << 7);
  GL_ptcl_z[0] = 0.0; // + _User_BoxSize_ / (1 << 7);
  GL_ptcl_vx[0] = 0.0; // 1.0227299814943555e-07 [kpc/year] = 100 km/s
  GL_ptcl_vy[0] = 0.0; //
  GL_ptcl_vz[0] = 0.0; //

  //Test particle
  
  GL_ptcl_mass[1] = 0.0; // 3.00261437e-06; // 3.0026143790849676e-06 = solar mass at units of Msun
  sum += GL_ptcl_mass[1]; 
  GL_ptcl_x[1] = 		4.8476302883992e-9 ; // 1 au = 4.8476302883992e-9L  kpc 
  GL_ptcl_y[1] = 0.0;
  GL_ptcl_z[1] = 0.0;
  vtype G_physical = 6.67430e-11 / (_kpc_to_m_ * _kpc_to_m_ * _kpc_to_m_) * _Msolar_to_kg_ * _year_to_s_ * _year_to_s_; // G [kpc^3 Msun^-1 year^-2]
  printf("physical G = %1.12e\n",(double)G_physical);
  //vtype distance = mysqrt(GL_ptcl_x[1]*GL_ptcl_x[1] + GL_ptcl_y[1]*GL_ptcl_y[1] + GL_ptcl_z[1]*GL_ptcl_z[1]);
  vtype distance = GL_ptcl_x[1] ;
  GL_ptcl_vx[1] = 0.0;
  GL_ptcl_vy[1] = mysqrt(G_physical * GL_ptcl_mass[0]/distance); // velocity at kpc/year = 6.283066641487499
  GL_ptcl_vz[1] = 0.0;
  printf("GL_ptcl_vy[1] = %1.12e\n",(double) GL_ptcl_vy[1]);

  vtype period = 2.0 * _PI_ * mysqrt(GL_ptcl_x[1] * GL_ptcl_x[1] * GL_ptcl_x[1] /(G_physical * GL_ptcl_mass[0] ) );
  printf("Period of the Test particle = %1.12e years\n", (double)period);
  printf("expected number of cycles in the tota time = %f\n",(double) ((Maxdt/(_conversion_time_))/period));


  if(sum == 0)
  {
    printf("Error at function input(), the total mass is cero\n");
    return _FAILURE_;
  }

  return _SUCCESS_;
}

void computing_cm_and_initial_mass(void)
{

  //Initial mass
  GL_total_mass_initial = 0.0;
  for (int i = 0; i < GL_no_ptcl_initial; i++)
  {
    GL_total_mass_initial += GL_ptcl_mass[i];
  }

  //Center of Mass
  GL_cm_x = 0.0;
  GL_cm_y = 0.0;
  GL_cm_z = 0.0;
  for (int i = 0; i < GL_no_ptcl_final; i++)
  {
    GL_cm_x += GL_ptcl_mass[i] * GL_ptcl_x[i];
    GL_cm_y += GL_ptcl_mass[i] * GL_ptcl_y[i];
    GL_cm_z += GL_ptcl_mass[i] * GL_ptcl_z[i];
  }
  //* Normalizing per the total mass *//
  GL_cm_x = GL_cm_x / GL_total_mass_initial;
  GL_cm_y = GL_cm_y / GL_total_mass_initial;
  GL_cm_z = GL_cm_z / GL_total_mass_initial;

  //Center of velocity
  GL_cm_vx = 0.0;
  GL_cm_vy = 0.0;
  GL_cm_vz = 0.0;
  for (int i = 0; i < GL_no_ptcl_final; i++)
  {
    GL_cm_vx += GL_ptcl_mass[i] * GL_ptcl_vx[i];
    GL_cm_vy += GL_ptcl_mass[i] * GL_ptcl_vy[i];
    GL_cm_vz += GL_ptcl_mass[i] * GL_ptcl_vz[i];
  }
  //* Normalizing per the total mass *//
  GL_cm_vx = GL_cm_vx / GL_total_mass_initial;
  GL_cm_vy = GL_cm_vy / GL_total_mass_initial;
  GL_cm_vz = GL_cm_vz / GL_total_mass_initial;
}


static int input_code_units(void)
{

  computing_cm_and_initial_mass();
  // printf("INITIAL VELOCITIES \n\n\n");
  for (int i = 0; i < GL_no_ptcl_initial; i++)
  {

    GL_ptcl_x[i] = (GL_ptcl_x[i] - GL_cm_x) * _conversion_dist_ + 0.5;
    GL_ptcl_y[i] = (GL_ptcl_y[i] - GL_cm_y) * _conversion_dist_ + 0.5;
    GL_ptcl_z[i] = (GL_ptcl_z[i] - GL_cm_z) * _conversion_dist_ + 0.5;

    GL_ptcl_vx[i] = (GL_ptcl_vx[i] - GL_cm_vx) * _conversion_velocity_;
    GL_ptcl_vy[i] = (GL_ptcl_vy[i] - GL_cm_vy) * _conversion_velocity_;
    GL_ptcl_vz[i] = (GL_ptcl_vz[i] - GL_cm_vz) * _conversion_velocity_;

    if ((0.0 >= GL_ptcl_x[i]) || (1.0 <= GL_ptcl_x[i]) || (0.0 >= GL_ptcl_y[i]) || (1.0 <= GL_ptcl_y[i]) || (0.0 >= GL_ptcl_z[i]) || (1.0 <= GL_ptcl_z[i]) )
    {
      printf("\nERROR: The input particle ID = %d is outside of the simulation box \n",i);
      printf("\nParticle position at code coordinates: x = %f, y = %f, z = %f, Size of the box = %f^3\n",(double) GL_ptcl_x[i],(double) GL_ptcl_y[i],(double) GL_ptcl_z[i], (double)_User_BoxSize_);
      return _FAILURE_;
    }
  }
  return _SUCCESS_;
}

int input(void)
{
  // if (input_plummer_model() == _FAILURE_)
  // {
  //   printf("\n\n Error running input_plummer_model function\n\n");
  //   return _FAILURE_;
  // }

  // if (input_galaxies_merger() == _FAILURE_)
  // {
  //     printf("\n\n Error running input_galaxies_merger function\n\n");
  //     return _FAILURE_;
  // }

 if (input_particle_initialization() == _FAILURE_)
  {
    printf("\n\n Error running input_particle_initialization function\n\n");
    return _FAILURE_;
  } 

  if (input_code_units() == _FAILURE_)
  {
    printf("\n\n Error running input_codedimensions function\n\n");
    return _FAILURE_;
  }

  return _SUCCESS_;
}


