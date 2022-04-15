/*
 * potential_heada_node.c
 *
 * Compute the potential in the head node or coarsest level of refinement
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

#include "potential_head_node.h"

static vtype *multigrid_restriction(const vtype *phi, int cgridsize)
{

	int gridsize = 2 * cgridsize - 1;
	int gridsize_pow_2 = gridsize * gridsize;
	int cgridsize_pow_2 = cgridsize * cgridsize;
	int Index_A, Index_B;

	vtype *cphi;
	cphi = (vtype *)calloc(cgridsize * cgridsize * cgridsize, sizeof(vtype));

	for (int lz = 1; lz < cgridsize - 1; lz++)
	{
		for (int ly = 1; ly < cgridsize - 1; ly++)
		{
			for (int lx = 1; lx < cgridsize - 1; lx++)
			{
				Index_A = lz * cgridsize_pow_2 + ly * cgridsize + lx;
				Index_B = 2 * lz * gridsize_pow_2 + 2 * ly * gridsize + 2 * lx;
				cphi[Index_A] = phi[Index_B];
			}
		}
	}
	return cphi;
}

static vtype *multigrid_projection(const vtype *dphic, int cgridsize)
{

	vtype *dphi;
	int Index_A, Index_B;

	int gridsize = 2 * cgridsize - 1;
	int gridsize_pow_2 = gridsize * gridsize;
	int cgridsize_pow_2 = cgridsize * cgridsize;

	dphi = (vtype *)calloc(gridsize * gridsize * gridsize, sizeof(vtype));

	for (int lz = 1; lz < cgridsize - 1; lz++)
	{
		for (int ly = 1; ly < cgridsize - 1; ly++)
		{
			for (int lx = 1; lx < cgridsize - 1; lx++)
			{
				Index_A = 2 * lz * gridsize_pow_2 + 2 * ly * gridsize + 2 * lx;
				Index_B = lz * cgridsize_pow_2 + ly * cgridsize + lx;
				// There are 4 cycles. The size are 1,6,12 and 8.
				// They are organized according the distance to the coarse grid point.

				// Type 1
				dphi[Index_A] = dphic[Index_B];
				// Type 2
				dphi[Index_A + 1] += 0.5 * dphic[Index_B];
				dphi[Index_A - 1] += 0.5 * dphic[Index_B];
				dphi[Index_A + gridsize] += 0.5 * dphic[Index_B];
				dphi[Index_A - gridsize] += 0.5 * dphic[Index_B];
				dphi[Index_A + gridsize_pow_2] += 0.5 * dphic[Index_B];
				dphi[Index_A - gridsize_pow_2] += 0.5 * dphic[Index_B];
				// Type 3
				dphi[Index_A + gridsize + 1] += 0.25 * dphic[Index_B];
				dphi[Index_A + gridsize - 1] += 0.25 * dphic[Index_B];
				dphi[Index_A - gridsize + 1] += 0.25 * dphic[Index_B];
				dphi[Index_A - gridsize - 1] += 0.25 * dphic[Index_B];
				dphi[Index_A + gridsize_pow_2 + 1] += 0.25 * dphic[Index_B];
				dphi[Index_A + gridsize_pow_2 - 1] += 0.25 * dphic[Index_B];
				dphi[Index_A - gridsize_pow_2 + 1] += 0.25 * dphic[Index_B];
				dphi[Index_A - gridsize_pow_2 - 1] += 0.25 * dphic[Index_B];
				dphi[Index_A + gridsize_pow_2 + gridsize] += 0.25 * dphic[Index_B];
				dphi[Index_A + gridsize_pow_2 - gridsize] += 0.25 * dphic[Index_B];
				dphi[Index_A - gridsize_pow_2 + gridsize] += 0.25 * dphic[Index_B];
				dphi[Index_A - gridsize_pow_2 - gridsize] += 0.25 * dphic[Index_B];
				// Type 4
				dphi[Index_A + gridsize_pow_2 + gridsize + 1] += 0.125 * dphic[Index_B];
				dphi[Index_A + gridsize_pow_2 + gridsize - 1] += 0.125 * dphic[Index_B];
				dphi[Index_A + gridsize_pow_2 - gridsize + 1] += 0.125 * dphic[Index_B];
				dphi[Index_A + gridsize_pow_2 - gridsize - 1] += 0.125 * dphic[Index_B];
				dphi[Index_A - gridsize_pow_2 + gridsize + 1] += 0.125 * dphic[Index_B];
				dphi[Index_A - gridsize_pow_2 + gridsize - 1] += 0.125 * dphic[Index_B];
				dphi[Index_A - gridsize_pow_2 - gridsize + 1] += 0.125 * dphic[Index_B];
				dphi[Index_A - gridsize_pow_2 - gridsize - 1] += 0.125 * dphic[Index_B];
			}
		}
	}
	return dphi;
}

static void gauss_saidel(vtype *phi, const vtype *rho, int gridsize, int iter_max)
{

	vtype H_pow2 = (1.0L / (gridsize - 1)) * (1.0L / (gridsize - 1));

	int grid_idx;

	int gridsize_pow_2 = gridsize * gridsize;
	int C1 = (gridsize - 1);
	int C2 = gridsize * (gridsize - 1);
	int C3 = gridsize_pow_2 * (gridsize - 1);

	for (int iter = 0; iter < iter_max; iter++)
	{
		for (int k = gridsize_pow_2; k < C3; k += gridsize_pow_2)
		{
			for (int j = gridsize; j < C2; j += gridsize)
			{
				for (int i = 1; i < C1; i++)
				{
					grid_idx = k + j + i;
					phi[grid_idx] = _Onesixth_ * (-H_pow2 * rho[grid_idx] + phi[grid_idx + 1] + phi[grid_idx - 1] + phi[grid_idx + gridsize] + phi[grid_idx - gridsize] + phi[grid_idx + gridsize_pow_2] + phi[grid_idx - gridsize_pow_2]);
				}
			}
		}
	}
}

// Computation Jacobi function
static void compute_jacobi(vtype *phi_old, vtype* phi_new, const vtype *rho,
						  int gridsize, int iter_max)
{

	vtype *pointer_auxiliary;
	vtype C = -(1.0L / (gridsize - 1)) * (1.0L / (gridsize - 1));
	int gridsize_pow_2 = gridsize * gridsize;
	int C1 = (gridsize - 1);
	int C2 = gridsize * (gridsize - 1);
	int C3 = gridsize_pow_2 * (gridsize - 1);
	int Index;

	for (int k = 0; k < iter_max; k++)
	{
		for (int lz = gridsize_pow_2; lz < C3; lz += gridsize_pow_2)
		{
			for (int ly = gridsize; ly < C2; ly += gridsize)
			{
				for (int lx = 1; lx < C1; lx++)
				{
					Index = lz + ly + lx;
					phi_new[Index] = _Onesixth_ * (C * rho[Index] + phi_old[Index + 1] + phi_old[Index - 1] + phi_old[Index + gridsize] + phi_old[Index - gridsize] + phi_old[Index + gridsize_pow_2] + phi_old[Index - gridsize_pow_2]);
				}
			}
		}
		pointer_auxiliary = phi_old;
		phi_old = phi_new;
		phi_new = pointer_auxiliary;
	}
}

static void jacobi(vtype *phi, const vtype *rho, int gridsize, int iter_max)
{

	vtype *phi_copy;
	int gridsize_pow_3 = gridsize * gridsize * gridsize;

	// Copy of phi
	phi_copy = (vtype *)malloc(gridsize_pow_3 * sizeof(vtype));
	memcpy(phi_copy, phi, gridsize_pow_3 * sizeof(vtype));

	// It depends on what buff update first
	if (iter_max % 2 == 0)
	{
		compute_jacobi(phi, phi_copy, rho, gridsize, iter_max);
		
	}
	else
	{
		compute_jacobi(phi_copy, phi, rho, gridsize, iter_max + 1);
	}
	free(phi_copy);
}

static void multigrid_solver(vtype *phi, const vtype *rho, int gridsize, int type_multigrid)
{
	int grid_idx;
	vtype *dphi;
	vtype *dphic;
	vtype *rest;
	vtype *restcoarse;

	int gridsize_pow_2 = gridsize * gridsize;
	int gridsize_pow_3 = gridsize_pow_2 * gridsize;
	int gridsizecoarse = (int)(gridsize + 1) / 2;
	int gridsizecoarse_pow_3 = gridsizecoarse * gridsizecoarse * gridsizecoarse;
	vtype one_over_H2 = (gridsize - 1) * (gridsize - 1);
	int C1 = (gridsize - 1);
	int C2 = gridsize * (gridsize - 1);
	int C3 = gridsize_pow_2 * (gridsize - 1);

	dphic = (vtype *)calloc(gridsizecoarse_pow_3, sizeof(vtype));
	rest = (vtype *)calloc(gridsize_pow_3, sizeof(vtype));

	//** >> Pre - Smoothing **/
	if(solver == 0)
	{
		gauss_saidel(phi, rho, gridsize, _NiterPreS_);
	}
	else if (solver == 1)
	{
		jacobi(phi, rho, gridsize, _NiterPreS_);
	}

	//** >> Finding rest **/
	// rest_i = rho_i - Lap[phi_i]
	for (int k = gridsize_pow_2; k < C3; k += gridsize_pow_2)
	{
		for (int j = gridsize; j < C2; j += gridsize)
		{
			for (int i = 1; i < C1; i++)
			{
				grid_idx = k + j + i;
				rest[grid_idx] = rho[grid_idx] + one_over_H2 * (

															 6 * phi[grid_idx] - phi[grid_idx + 1] - phi[grid_idx - 1] - phi[grid_idx + gridsize] - phi[grid_idx - gridsize] - phi[grid_idx + gridsize_pow_2] - phi[grid_idx - gridsize_pow_2]);
			}
		}
	}

	//** >> Restriction: Interpoling fine to coarse **/
	restcoarse = multigrid_restriction(rest, gridsizecoarse);

	//** >> Finding dphic in coarse level **/
	if (gridsizecoarse < 4)
	{
		//** >> SOLVER **/
		if (solver == 0)
		{
			gauss_saidel(phi, rho, gridsize, _Niterfinddphic_);
		}
		else if (solver == 1)
		{
			jacobi(phi, rho, gridsize, _Niterfinddphic_);
		}
	}
	else
	{
		// Here is the onj difference between the multigrid methods
		if(type_multigrid == 0)
		{
			multigrid_solver(dphic, restcoarse, gridsizecoarse, 0);
		}
		else if (type_multigrid == 1)
		{
			multigrid_solver(dphic, restcoarse, gridsizecoarse, 1);
			multigrid_solver(dphic, restcoarse, gridsizecoarse, 0);
		}
		else if (type_multigrid == 2)
		{
			multigrid_solver(dphic, restcoarse, gridsizecoarse, 2);
			multigrid_solver(dphic, restcoarse, gridsizecoarse, 2);
		}
		
	}

	//** >> Projection: Interpolation coarse to fine **/
	dphi = multigrid_projection(dphic, gridsizecoarse);

	//** >> Adding dphi in solution phi = phi + dphi **/
	for (int k = gridsize_pow_2; k < C3; k += gridsize_pow_2)
	{
		for (int j = gridsize; j < C2; j += gridsize)
		{
			for (int i = 1; i < C1; i++)
			{
				grid_idx = k + j + i;
				phi[grid_idx] += dphi[grid_idx];
			}
		}
	}

	//** >> Post - Smoothing **/
	if (solver == 0)
	{
		gauss_saidel(phi, rho, gridsize, _NiterPostS_);
	}
	else if (solver == 1)
	{
		jacobi(phi, rho, gridsize, _NiterPostS_);
	}

	free(dphic);
	free(rest);
	free(restcoarse);
	free(dphi);
}

static int compute_potential_head_node(struct node *ptr_node_pt)
{
	int lv;

	lv = ptr_node_pt->lv;

	if(lv == lmin)
	{
		int grid_side;	// Number of grid points per side
		int grid_side_pow2;
		int grid_limit; // Superior border in the internal grid box

		vtype *aux_pot;
		vtype *aux_d;
		int size;
		int size_pow2;
		int size_pow3;

		int aux_idx;
		int box_grid_idx;

		size = no_lmin_cell + 1;
		size_pow2 = size * size;
		size_pow3 = size_pow2 * size;
		aux_pot = (vtype *)malloc(size_pow3 * sizeof(vtype));
		aux_d = (vtype *)malloc(size_pow3 * sizeof(vtype));

		grid_side = box_side_lmin + 1;
		grid_side_pow2 = grid_side * grid_side;
		grid_limit = grid_side - bder_os_sim - 1;

		//** >> Pasing information from real parameters to auxiliary parameters **/
		for (int k = bder_os_sim; k < grid_limit + 1; k++)
		{
			for (int j = bder_os_sim; j < grid_limit + 1; j++)
			{
				for (int i = bder_os_sim; i < grid_limit + 1; i++)
				{
					box_grid_idx = i + j * grid_side + k * grid_side_pow2;
					aux_idx = (i - bder_os_sim) + (j - bder_os_sim) * size + (k - bder_os_sim) * size_pow2;
					aux_pot[aux_idx] = ptr_node_pt->ptr_pot[box_grid_idx];
					aux_d[aux_idx] = ptr_node_pt->ptr_d[box_grid_idx];
				}
			}
		}

		int type_multigrid = multigrid;
		//** >> Poisson Solver **/
		// Potential-Density-gridsize-type of multigrid
		multigrid_solver(aux_pot, aux_d, size, type_multigrid);

		//** >> Pasing information from auxiliary parameters to real parameters **/
		for (int k = bder_os_sim; k < grid_limit + 1; k++)
		{
			for (int j = bder_os_sim; j < grid_limit + 1; j++)
			{
				for (int i = bder_os_sim; i < grid_limit + 1; i++)
				{
					box_grid_idx = i + j * grid_side + k * grid_side_pow2;
					aux_idx = (i - bder_os_sim) + (j - bder_os_sim) * size + (k - bder_os_sim) * size_pow2;
					ptr_node_pt->ptr_pot[box_grid_idx] = aux_pot[aux_idx];
					ptr_node_pt->ptr_d[box_grid_idx] = aux_d[aux_idx];
				}
			}
		
		}
		free(aux_pot);
		free(aux_d);
	}

	return _SUCCESS_;
}

int potential_head_node()
{

	struct node *ptr_node_pt;

	ptr_node_pt = GL_ptr_tree;

	//** >> Coarsest potential **/
	if(compute_potential_head_node(ptr_node_pt) == _FAILURE_)
	{
		printf("Error at function potential\n");
		return _FAILURE_;
	}

	return _SUCCESS_;
}