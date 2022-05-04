/*
 * poisson-error.c
 *
 * Compute the error of the solution of the Poisson equation potential
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

#include "poisson_error.h"

static bool interior_grid_point(const struct node *ptr_node, int i)
{
	bool check;

	check = true;

	int box_idx_x;
	int box_idx_y;
	int box_idx_z;
	int box_idx;

	box_idx_x = ptr_node->ptr_cell_idx_x[i] - ptr_node->box_ts_x;
	box_idx_y = ptr_node->ptr_cell_idx_y[i] - ptr_node->box_ts_y;
	box_idx_z = ptr_node->ptr_cell_idx_z[i] - ptr_node->box_ts_z;
	box_idx = box_idx_x + box_idx_y * ptr_node->box_real_dim_x + box_idx_z * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;

	//** Checking the nearest 3 neighbors of face
	if (ptr_node->ptr_box[box_idx - 1] < -3)
	{
		check = false;
	}
	else if (ptr_node->ptr_box[box_idx - ptr_node->box_real_dim_x] < -3)
	{
		check = false;
	}
	else if (ptr_node->ptr_box[box_idx - ptr_node->box_real_dim_x * ptr_node->box_real_dim_y] < -3)
	{
		check = false;
	}

	return check;
}

static bool poisson_error_mehod_0(const struct node *ptr_node)
{
/* 	The condition used is option 1 (default one), and it
   uses the root mean square. In general, the error in the grid of the Poisson
   equation solution is computed as follows:

   error = Sqrt[ 1/Ng^3 * Sum[{( density_i - Laplacian(pot_i) )/rhomean }^2 ] ]

   The criterion is accepted when the error is less than the threshold. See the
   _ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_ parameter in
   Precision_Parameters.h. Note that error is an adimensional variable. */

	vtype error; // Error of the potential solution
	vtype diff;	 // Diference between the Laplacian of the potential solution and the density

	vtype H;
	vtype one_over_H_pow_2;

	bool check; // Check if the grid point is an interior grid point

	int size; // Number of cells in the node
	int cntr; // Counter of the number of interior grid points

	int box_grid_idx_x;
	int box_grid_idx_y;
	int box_grid_idx_z;
	int box_grid_idx;

	vtype rhomean_times_4piG;

	H = 1.0L / (1 << ptr_node->lv);
	one_over_H_pow_2 = 1.0L / (H * H);
	size = ptr_node->cell_size;

	//printf("local_mass = %1.3e\n", ptr_node->local_mass);
	rhomean_times_4piG = 4 * _G_ * _PI_ * ptr_node->local_mass / (size * H * H * H);

	error = 0;
	cntr = 0;
	
	// Computing the root mean square normalized to the mean density rhomean
	for (int i = 0; i < size; i++)
	{
		//** >> Checking if the grid point is a interior grid point **/
		check = interior_grid_point(ptr_node,i);
		if (check == true)
		{
			box_grid_idx_x = ptr_node->ptr_cell_idx_x[i] - ptr_node->box_ts_x;
			box_grid_idx_y = ptr_node->ptr_cell_idx_y[i] - ptr_node->box_ts_y;
			box_grid_idx_z = ptr_node->ptr_cell_idx_z[i] - ptr_node->box_ts_z;
			box_grid_idx = box_grid_idx_x + box_grid_idx_y * (ptr_node->box_real_dim_x + 1) + box_grid_idx_z * (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1);

			diff = ptr_node->ptr_d[box_grid_idx] + one_over_H_pow_2 * 
			(
			6.0 * ptr_node->ptr_pot[box_grid_idx] 
			- ptr_node->ptr_pot[box_grid_idx + 1] 
			- ptr_node->ptr_pot[box_grid_idx - 1] 
			- ptr_node->ptr_pot[box_grid_idx + (ptr_node->box_real_dim_x + 1)] 
			- ptr_node->ptr_pot[box_grid_idx - (ptr_node->box_real_dim_x + 1)]  
			- ptr_node->ptr_pot[box_grid_idx + (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1) ] 
			- ptr_node->ptr_pot[box_grid_idx - (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1) ] 
			);
			error += diff * diff;
			cntr++;
			
		}
	}
	error = error / cntr;
	error = sqrt(error);
	error = error / rhomean_times_4piG;

	//** >> If the precision condition satisfied **/
	if (error < _ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_)
	{
		return true;
	}
	else
	{
		return false;
	}
}

static bool poisson_error_mehod_1(const struct node *ptr_node)
{

/*    This module is responsible for checking the condition for the solution of the
   Poisson equation. The condition used is option 2, and it is the absolute
   error between the density and the laplacian for each grid point. In general,
   the error in each grid point of the Poisson equation solution is computed as
   follows:

   error_i =  |density_i - Laplacian(pot_i)| /rhomean

   The criterion is accepted only if all errors in the grid points are less than
   the threshold. See the _ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_ parameter in
   Precision_Parameters.h. Note that error is an adimensional variable. */

	vtype error; // Error of the potential solution

	vtype H;
	vtype one_over_H_pow_2;

	bool check; // Check if the grid point is an interior grid point

	int size; // Number of cells in the node
	int cntr; // Counter of the number of interior grid points

	int box_grid_idx_x;
	int box_grid_idx_y;
	int box_grid_idx_z;
	int box_grid_idx;

	vtype rhomean_times_4piG;

	H = 1.0L / (1 << ptr_node->lv);
	one_over_H_pow_2 = 1.0L / (H * H);
	size = ptr_node->cell_size;

	rhomean_times_4piG = 4 * _G_ * _PI_ * ptr_node->local_mass / (size * H * H * H);

	error = 0;
	cntr = 0;

	// Computing the root mean square normalized to the mean density rhomean
	for (int i = 0; i < size; i++)
	{
		//** >> Checking if the grid point is a interior grid point **/
		check = interior_grid_point(ptr_node, i);
		if (check == true)
		{
			box_grid_idx_x = ptr_node->ptr_cell_idx_x[i] - ptr_node->box_ts_x;
			box_grid_idx_y = ptr_node->ptr_cell_idx_y[i] - ptr_node->box_ts_y;
			box_grid_idx_z = ptr_node->ptr_cell_idx_z[i] - ptr_node->box_ts_z;
			box_grid_idx = box_grid_idx_x + box_grid_idx_y * (ptr_node->box_real_dim_x + 1) + box_grid_idx_z * (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1);

			error = ptr_node->ptr_d[box_grid_idx] + one_over_H_pow_2 *
			(
			6.0 * ptr_node->ptr_pot[box_grid_idx] 
			- ptr_node->ptr_pot[box_grid_idx + 1] 
			- ptr_node->ptr_pot[box_grid_idx - 1] 
			- ptr_node->ptr_pot[box_grid_idx + (ptr_node->box_real_dim_x + 1)] 
			- ptr_node->ptr_pot[box_grid_idx - (ptr_node->box_real_dim_x + 1)] 
			- ptr_node->ptr_pot[box_grid_idx + (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1)] 
			- ptr_node->ptr_pot[box_grid_idx - (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1)]
			);

			error = myabs(error) / rhomean_times_4piG;
			cntr++;
			if (error >= _ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_)
			{
				return false;
			}
		}
	}

	return true;
}

static bool poisson_error_mehod_2(const struct node *ptr_node)
{
	/* 	   This module is responsible for checking the condition for the solution of the
   Poisson equation. The condition used is option 3, and it uses the root mean
   square. In general, the error in the grid of the Poisson equation solution
   is computed as follows:

   error =

   Sqrt[ 1/Ng^3 * Sum[{( density_i - Laplacian(pot_i) )/rhomean }^2 ] ] * 1/N^2

   The criterion is accepted when the error is less than the threshold. See the
   _ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_ parameter in
   Precision_Parameters.h. Note that error is an adimensional variable. */

	vtype error; // Error of the potential solution
	vtype diff;	 // Diference between the Laplacian of the potential solution and the density

	vtype H;
	vtype one_over_H_pow_2;

	bool check; // Check if the grid point is an interior grid point

	int size; // Number of cells in the node
	int cntr; // Counter of the number of interior grid points

	int box_grid_idx_x;
	int box_grid_idx_y;
	int box_grid_idx_z;
	int box_grid_idx;

	vtype rhomean_times_4piG;

	H = 1.0L / (1 << ptr_node->lv);
	one_over_H_pow_2 = 1.0L / (H * H);
	size = ptr_node->cell_size;

	rhomean_times_4piG = 4 * _G_ * _PI_ * ptr_node->local_mass / (size * H * H * H);

	error = 0;
	cntr = 0;

	// Computing the root mean square normalized to the mean density rhomean
	for (int i = 0; i < size; i++)
	{
		//** >> Checking if the grid point is a interior grid point **/
		check = interior_grid_point(ptr_node, i);
		if (check == true)
		{
			box_grid_idx_x = ptr_node->ptr_cell_idx_x[i] - ptr_node->box_ts_x;
			box_grid_idx_y = ptr_node->ptr_cell_idx_y[i] - ptr_node->box_ts_y;
			box_grid_idx_z = ptr_node->ptr_cell_idx_z[i] - ptr_node->box_ts_z;
			box_grid_idx = box_grid_idx_x + box_grid_idx_y * (ptr_node->box_real_dim_x + 1) + box_grid_idx_z * (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1);

			diff = ptr_node->ptr_d[box_grid_idx] + one_over_H_pow_2 *
													   (6.0 * ptr_node->ptr_pot[box_grid_idx] - ptr_node->ptr_pot[box_grid_idx + 1] - ptr_node->ptr_pot[box_grid_idx - 1] - ptr_node->ptr_pot[box_grid_idx + (ptr_node->box_real_dim_x + 1)] - ptr_node->ptr_pot[box_grid_idx - (ptr_node->box_real_dim_x + 1)] - ptr_node->ptr_pot[box_grid_idx + (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1)] - ptr_node->ptr_pot[box_grid_idx - (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1)]);
			error += diff * diff;
			cntr++;
		}
	}
	error = error / cntr;
	error = sqrt(error);
	error = error / rhomean_times_4piG;
	// Following line is the only difference with poisson_error_mehod_0
	error = error / (cbrt(cntr) * cbrt(cntr));

	//** >> If the precision condition satisfied **/
	if (error < _ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_)
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool poisson_error(struct node *ptr_node)
{

	bool check = false;

	if (check_poisson_error_method == 0)
	{
		check = poisson_error_mehod_0(ptr_node);
	}
	else if (check_poisson_error_method == 1)
	{
		check = poisson_error_mehod_1(ptr_node);
	}
	else if (check_poisson_error_method == 2)
	{
		check = poisson_error_mehod_2(ptr_node);
	}
	else
	{
		printf("Error: check_poisson_error_method is different to 0, 1 or 2\n");
		exit(EXIT_FAILURE);
	}

	return check;
}