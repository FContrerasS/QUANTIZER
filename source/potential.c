/*
 * potential.c
 *
 * Potential computation in each node
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

#include "potential.h"

int potential(void)
{

    //**>> POTENTIAL HEAD NODE **/
    if (potential_head_node() == _FAILURE_)
    {
        printf("\n\n Error running potential_head_node() function \n\n ");
        return _FAILURE_;
    }

    //** >> POTENTIAL BRANCH NODES **/
    if (lmin < lmax)
    {
        //**>> POTENTIAL HEAD NODE **/
        if (potential_branches() == _FAILURE_)
        {
            printf("\n\n Error running potential_branches() function \n\n ");
            return _FAILURE_;
        }
    }

    return _SUCCESS_;
}