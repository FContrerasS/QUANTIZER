/*
 * output_snapshots.c
 *
 * Export snapshots of the simulation
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

#include "output_snapshots.h"

//* >> Local Functions
static int output_folder_names(void);
static int create_output_folders(void);

char file_data_name[1000];
char folder_name[1000];

// Here we should read a file with the name information
const char folder_main_output_name_read[] = " ";
const char folder_Simulation_output_name_read[] = " ";
const char file_snashot_output_name_read[] = " ";

char folder_main_output_name[100];
char folder_Simulation_output_name[100];
char file_snashot_output_name[100];

// Creating the simulation directory
char address_simulation_directory[300];

static int output_folder_names(void)
{
  //----------------------------------------------------------------------------//

  // Puting the name of the Main folder
  if (strcmp(folder_main_output_name_read, " ") == 0)
  {
    // Default value
    strcpy(folder_main_output_name, "../output");
  }
  else
  {
    // User value
    sprintf(folder_main_output_name, "../%s", folder_main_output_name_read);
  }

  //----------------------------------------------------------------------------//

  // Puting the name of the Main folder
  if (strcmp(folder_Simulation_output_name_read, " ") == 0)
  {
    // Default value (total time, Ng, Np, )
    char aux_char[100];
    sprintf(aux_char, "Simulation_Maxdt[My](%1.1f)_lv(%d+%d)_Np(%d)", (double)(Maxdt / _Mgyear_), lmin, lmax - lmin, GL_no_ptcl_initial);
    strcpy(folder_Simulation_output_name, aux_char);
  }
  else
  {
    // User value
    strcpy(folder_Simulation_output_name, folder_Simulation_output_name_read);
  }

  //----------------------------------------------------------------------------//

  // Puting the name of the data snapshots
  if (strcmp(file_snashot_output_name_read, " ") == 0)
  {
    // Default value
    strcpy(file_snashot_output_name, "Data");
  }
  else
  {
    // User value
    strcpy(file_snashot_output_name, file_snashot_output_name_read);
  }
  return _SUCCESS_;
}

static int create_output_folders(void)
{
  char aux_char[500];
  char aux_char2[600];

  // Creating Main Directory
  mkdir(&folder_main_output_name[0], 0777);

  // adding Output main + Simulation folders
  sprintf(address_simulation_directory, "%s/%s", folder_main_output_name, folder_Simulation_output_name);

  // Creating the folder if check == 0, if not add a Release version
  int check = mkdir(&address_simulation_directory[0], 0777);
  if (check == -1)
  {

    strcat(address_simulation_directory, "_Release(");
    int k = 0;
    while (check == -1)
    {
      k++;
      sprintf(aux_char, "%s%d)", address_simulation_directory, k);
      check = mkdir(&aux_char[0], 0777);
    }
    strcpy(address_simulation_directory, aux_char);
  }

  // Address of the Snapshots data folder
  sprintf(aux_char2, "%s/Data", address_simulation_directory);

  // Creating the data directory inside the simulation for storing snapshots
  mkdir(&aux_char2[0], 0777);

  // Defining the Final address of the Folder Simulation
  strcpy(folder_name, address_simulation_directory);

  // Address of the snapshots files without the number of the snapshot. See
  // Export_Snapshot_Data.c
  sprintf(file_data_name, "%s/%s_", aux_char2, file_snashot_output_name);

  folder_created = true;

  return _SUCCESS_;
}

int output_snapshots(const vtype *energies, vtype actualtime, int snapshot)
{

  int size = GL_no_ptcl_final * sizeof(vtype);
  char snapshot_name[1500];

  vtype time_Megayear = actualtime / _Mgyear_; // Actual time in mega years

  if (folder_created == false)
  { // The folder are created only one time
    if (output_folder_names() == _FAILURE_)
    {
      printf("\n\n Error running output_dir_name function\n\n");
      return _FAILURE_;
    }

    if (create_output_folders() == _FAILURE_)
    {
      printf("\n\n Error running output_dir_name function\n\n");
      return _FAILURE_;
    }
  }

  vtype _one_over_User_sqrt_BoxSize_ = 1.0L / sqrt(_User_BoxSize_);

  vtype *GL_ptcl_x_conversion;
  vtype *GL_ptcl_y_conversion;
  vtype *GL_ptcl_z_conversion;
  vtype *GL_ptcl_vx_conversion;
  vtype *GL_ptcl_vy_conversion;
  vtype *GL_ptcl_vz_conversion;

  GL_ptcl_x_conversion = (vtype *)malloc(size);
  GL_ptcl_y_conversion = (vtype *)malloc(size);
  GL_ptcl_z_conversion = (vtype *)malloc(size);
  GL_ptcl_vx_conversion = (vtype *)malloc(size);
  GL_ptcl_vy_conversion = (vtype *)malloc(size);
  GL_ptcl_vz_conversion = (vtype *)malloc(size);

  // Returning to user units
  for (int i = 0; i < GL_no_ptcl_final; i++)
  {
    GL_ptcl_x_conversion[i] = (GL_ptcl_x[i] - 0.5) * _User_BoxSize_;
    GL_ptcl_y_conversion[i] = (GL_ptcl_y[i] - 0.5) * _User_BoxSize_;
    GL_ptcl_z_conversion[i] = (GL_ptcl_z[i] - 0.5) * _User_BoxSize_;

    GL_ptcl_vx_conversion[i] = GL_ptcl_vx[i] * _one_over_User_sqrt_BoxSize_;
    GL_ptcl_vy_conversion[i] = GL_ptcl_vy[i] * _one_over_User_sqrt_BoxSize_;
    GL_ptcl_vz_conversion[i] = GL_ptcl_vz[i] * _one_over_User_sqrt_BoxSize_;
  }

  sprintf(snapshot_name, "%s%d.bin", file_data_name, snapshot);

  FILE *file = NULL;
  file = fopen(snapshot_name, "w");

  if (file == NULL)
  {
    printf("Error opening snapshot file!\n");
    exit(EXIT_FAILURE);
  }

  fwrite(GL_ptcl_mass, size, 1, file);
  fwrite(GL_ptcl_x_conversion, size, 1, file);
  fwrite(GL_ptcl_y_conversion, size, 1, file);
  fwrite(GL_ptcl_z_conversion, size, 1, file);
  fwrite(GL_ptcl_vx_conversion, size, 1, file);
  fwrite(GL_ptcl_vy_conversion, size, 1, file);
  fwrite(GL_ptcl_vz_conversion, size, 1, file);
  fwrite(energies, 3 * sizeof(vtype), 1, file);
  fwrite(&time_Megayear, sizeof(vtype), 1, file);
  fclose(file);

  // Free memory
  free(GL_ptcl_x_conversion);
  free(GL_ptcl_y_conversion);
  free(GL_ptcl_z_conversion);
  free(GL_ptcl_vx_conversion);
  free(GL_ptcl_vy_conversion);
  free(GL_ptcl_vz_conversion);

  return _SUCCESS_;
}
