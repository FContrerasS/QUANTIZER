/*
 * multigrid.c
 *
 * Solve the Discrete Poisson equation using the multigrid method
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

#include "multigrid.h"

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

void multigrid_solver(vtype *phi, const vtype *rho, int gridsize, int type_multigrid)
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
    if (solverPreS == 0)
    {
        gauss_saidel(phi, rho, gridsize, _NiterPreS_);
    }
    else if (solverPreS == 1)
    {
        jacobi(phi, rho, gridsize, _NiterPreS_);
    }
    else if (solverPreS == 2)
    {
        conjugate_gradient_multigrid(phi, rho, gridsize, _NiterPreS_);
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
                rest[grid_idx] = rho[grid_idx] + one_over_H2 * (6 * phi[grid_idx] - phi[grid_idx + 1] - phi[grid_idx - 1] - phi[grid_idx + gridsize] - phi[grid_idx - gridsize] - phi[grid_idx + gridsize_pow_2] - phi[grid_idx - gridsize_pow_2]);
            }
        }
    }

    //** >> Restriction: Interpoling fine to coarse **/
    restcoarse = multigrid_restriction(rest, gridsizecoarse);

    //** >> Finding dphic in coarse level **/
    if (gridsizecoarse < 4)
    {
        //** >> SOLVER **/
        if (solverfinddphic == 0)
        {
            gauss_saidel(phi, rho, gridsize, _Niterfinddphic_);
        }
        else if (solverfinddphic == 1)
        {
            jacobi(phi, rho, gridsize, _Niterfinddphic_);
        }
        else if (solverfinddphic == 2)
        {
            conjugate_gradient_multigrid(phi, rho, gridsize, _Niterfinddphic_);
        }
    }
    else
    {
        // Here is the onj difference between the multigrid methods
        if (type_multigrid == 0)
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
    if (solverPostS == 0)
    {
        gauss_saidel(phi, rho, gridsize, _NiterPostS_);
    }
    else if (solverPostS == 1)
    {
        jacobi(phi, rho, gridsize, _NiterPostS_);
    }
    else if (solverPostS == 2)
    {
        conjugate_gradient_multigrid(phi, rho, gridsize, _NiterPostS_);
    }

    free(dphic);
    free(rest);
    free(restcoarse);
    free(dphi);
}
