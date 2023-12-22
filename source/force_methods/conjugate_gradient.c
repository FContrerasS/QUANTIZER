/*
 * conjugate_gradient.c
 *
 * Solve the Discrete Poisson equation using the conjugate gradient method
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

#include "conjugate_gradient.h"

int conjugate_gradient_oldversion(struct node *ptr_node)
{

  int number_iter = ptr_node->grid_intr_size;

  //printf("\nccomputing conjugate gradient\n");


  //int box_grid_idx_vec[8] = {150381,150382,150448,150449,154870,154871,154937,154938};
  // for(int hh = 0; hh < 8; hh++)
  // {
  //   printf("box_grid_idx = %d, pot = %.12f\n",box_grid_idx_vec[hh],(double)ptr_node->ptr_pot[box_grid_idx_vec[hh]]);
  // }



  //int number_of_threads = 8;

  //omp_set_num_threads(NUMBER_OF_THEADS);

  //clock_t aux_clock;
  //double aux_clock_omp;
  struct timespec start, finish;

  int box_grid_idx;

  vtype H = 1.0 / (1 << ptr_node->lv);
  vtype one_over_H_pow_2 = 1.0 / (H * H);
  int size = (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1) * (ptr_node->box_real_dim_z + 1);

  vtype *r, *p, *Ap;
  r = (vtype *)calloc(size, sizeof(vtype));
  p = (vtype *)calloc(size , sizeof(vtype));
  Ap = (vtype *)malloc(size * sizeof(vtype));

  vtype alpha, rsold, rsnew;

  int grid_box_real_dim_X = ptr_node->box_real_dim_x + 1;
  int grid_box_real_dim_X_times_Y = (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1);

  rsold = 0.0;
  

  
  //clock_gettime( CLOCK_REALTIME, &GL_start);
  #pragma omp parallel for private(box_grid_idx) reduction(+:rsold) 
  for (int i = 0; i < number_iter; i++)
  {
    //printf("my id = %d\n", omp_get_thread_num());
    box_grid_idx = ptr_node->ptr_intr_box_grid_idx[i];
    // Computing r
    // r[box_grid_idx] = ptr_node->ptr_d[box_grid_idx] + one_over_H_pow_2 *
    //                                                       (6.0 * ptr_node->ptr_pot[box_grid_idx] - ptr_node->ptr_pot[box_grid_idx + 1] - ptr_node->ptr_pot[box_grid_idx - 1] - ptr_node->ptr_pot[box_grid_idx + grid_box_real_dim_X] - ptr_node->ptr_pot[box_grid_idx - grid_box_real_dim_X] - ptr_node->ptr_pot[box_grid_idx + grid_box_real_dim_X_times_Y] - ptr_node->ptr_pot[box_grid_idx - grid_box_real_dim_X_times_Y]);
    r[box_grid_idx] = ptr_node->ptr_d[box_grid_idx] + one_over_H_pow_2 *
                                                      (6.0 * ptr_node->ptr_pot[box_grid_idx] - ptr_node->ptr_pot[box_grid_idx + 1] - ptr_node->ptr_pot[box_grid_idx - 1] - ptr_node->ptr_pot[box_grid_idx + grid_box_real_dim_X] - ptr_node->ptr_pot[box_grid_idx - grid_box_real_dim_X] - ptr_node->ptr_pot[box_grid_idx + grid_box_real_dim_X_times_Y] - ptr_node->ptr_pot[box_grid_idx - grid_box_real_dim_X_times_Y]);
  
    // Copy r to p
    p[box_grid_idx] = r[box_grid_idx];
    // Computing rsold
    rsold += r[box_grid_idx] * r[box_grid_idx];
  }
  //clock_gettime( CLOCK_REALTIME, &GL_finish);
  //GL_times[50] += ( GL_finish.tv_sec - GL_start.tv_sec ) + ( GL_finish.tv_nsec - GL_start.tv_nsec )/ 1000000000.;

  // // Copy r to p
  // #pragma omp parallel for
  // for(int i = 0; i < size; i++)
  // {
  //   p[i] = r[i];
  // }
  //memcpy(p, r, size * sizeof(vtype));
  

  // Computing rsold
  // rsold = 0;
  // #pragma omp parallel for private(box_grid_idx) reduction(+:rsold)
  // for (int i = 0; i < number_iter; i++)
  // {
  //   box_grid_idx = ptr_node->ptr_intr_box_grid_idx[i];
  //   rsold += r[box_grid_idx] * r[box_grid_idx];
  // }

  //printf("rsold = %f\n",rsold);

  // //* >> CHEKING ERROR SOLUTION CONDITION *//
  // aux_clock = clock();
  // if (poisson_error(ptr_node, r, rsold, 1) == true)
  // {
  //     GL_times[25] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;
  //     free(r);
  //     free(p);
  //     free(Ap);
  //     return _SUCCESS_;
  // }
  // GL_times[25] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;


  // printf("\n before the cycle begins\n");
  // for(int hh = 0; hh < 8; hh++)
  // {
  //   printf("box_grid_idx = %d, pot = %.12f\n",box_grid_idx_vec[hh],(double)ptr_node->ptr_pot[box_grid_idx_vec[hh]]);
  // }

  int cycle_size = size < _MAX_NUMBER_OF_ITERATIONS_IN_POISSON_EQUATION_ ? size : _MAX_NUMBER_OF_ITERATIONS_IN_POISSON_EQUATION_;
  // For cycle over mutually conjugate vectors p
  for (int i = 0; i < cycle_size; i++)
  {


    // Computing Ap and alpha
    alpha = 0.0;
    #pragma omp parallel for private(box_grid_idx) reduction(+:alpha) 
    for (int j = 0; j < number_iter; j++)
    {
      box_grid_idx = ptr_node->ptr_intr_box_grid_idx[j];
      Ap[box_grid_idx] = -one_over_H_pow_2 * (6.0 * p[box_grid_idx] - p[box_grid_idx + 1] - p[box_grid_idx - 1] - p[box_grid_idx + grid_box_real_dim_X] - p[box_grid_idx - grid_box_real_dim_X] - p[box_grid_idx + grid_box_real_dim_X_times_Y] - p[box_grid_idx - grid_box_real_dim_X_times_Y]);
      alpha += p[box_grid_idx] * Ap[box_grid_idx];
    }




    //printf("p[boxgrid] = %.12f\n",p[ptr_node->ptr_intr_box_grid_idx[5]] );

    // // Computing Ap
    // #pragma omp parallel for private(box_grid_idx)
    // for (int j = 0; j < number_iter; j++)
    // {
    //   box_grid_idx = ptr_node->ptr_intr_box_grid_idx[j];
    //   Ap[box_grid_idx] = -one_over_H_pow_2 * (6.0 * p[box_grid_idx] - p[box_grid_idx + 1] - p[box_grid_idx - 1] - p[box_grid_idx + grid_box_real_dim_X] - p[box_grid_idx - grid_box_real_dim_X] - p[box_grid_idx + grid_box_real_dim_X_times_Y] - p[box_grid_idx - grid_box_real_dim_X_times_Y]);
    //   // printf("Ap = %1.3Le\n", Ap[box_grid_idx]);
    // }

    // // Computing alpha
    // alpha = 0;
    // #pragma omp parallel for private(box_grid_idx) reduction(+:alpha)
    // for (int j = 0; j < number_iter; j++)
    // {
    //   box_grid_idx = ptr_node->ptr_intr_box_grid_idx[j];
    //   alpha += p[box_grid_idx] * Ap[box_grid_idx];
    // }
    
    //clock_gettime( CLOCK_REALTIME, &GL_start);
    alpha = rsold / alpha;
    
    rsnew = 0.0;

    // if(ptr_node->lv == 6)
    // {
    //   printf("\n cycle size = %d \n",i);
    //   // for(int hh = 0; hh < 8; hh++)
    //   // {
    //   //   printf("box_grid_idx = %d, pot = %.12f\n",box_grid_idx_vec[hh],(double)ptr_node->ptr_pot[box_grid_idx_vec[hh]]);
    //   // }
    // }
    

    #pragma omp parallel for private(box_grid_idx) reduction(+:rsnew)
    //#pragma omp parallel for private(box_grid_idx)
    for (int j = 0; j < number_iter; j++)
    {
      box_grid_idx = ptr_node->ptr_intr_box_grid_idx[j];
      // New solution
      ptr_node->ptr_pot[box_grid_idx] += alpha * p[box_grid_idx];
      // New rest
      r[box_grid_idx] -= alpha * Ap[box_grid_idx];
      // // Computing rsnew
      rsnew += r[box_grid_idx] * r[box_grid_idx];
    }

    printf("rsnew = %1.16e\n",(double)rsnew);


    //printf("hola = %d, rsnew = %f, alpha = %.12f, beta = %.12f\n",i,rsnew,alpha,rsnew/rsold);

    //clock_gettime( CLOCK_REALTIME, &GL_finish);
    //GL_times[50] += ( GL_finish.tv_sec - GL_start.tv_sec ) + ( GL_finish.tv_nsec - GL_start.tv_nsec )/ 1000000000.;






    // // New solution
    // #pragma omp parallel for private(box_grid_idx)
    // for (int j = 0; j < number_iter; j++)
    // {
    //   box_grid_idx = ptr_node->ptr_intr_box_grid_idx[j];
    //   ptr_node->ptr_pot[box_grid_idx] += alpha * p[box_grid_idx];
    // }

    // // New rest
    // #pragma omp parallel for private(box_grid_idx)
    // for (int j = 0; j < number_iter; j++)
    // {
    //   box_grid_idx = ptr_node->ptr_intr_box_grid_idx[j];
    //   r[box_grid_idx] -= alpha * Ap[box_grid_idx];
    // }

    // // Computing rsnew
    // rsnew = 0;
    // #pragma omp parallel for private(box_grid_idx) reduction(+:rsnew)
    // for (int j = 0; j < number_iter; j++)
    // {
    //   box_grid_idx = ptr_node->ptr_intr_box_grid_idx[j];
    //   rsnew += r[box_grid_idx] * r[box_grid_idx];
    // }











    //* >> CHEKING ERROR SOLUTION CONDITION *//
    if ((i+1) % iter_between_check_potential_solution == 0)
    {
      //printf("\ni = %d\n", i);
      //aux_clock = clock();
      //aux_clock_omp = omp_get_wtime();
      clock_gettime( CLOCK_REALTIME, &start);
      if (poisson_error(ptr_node, r, rsnew, 1) == true)
      {
          
        // if(ptr_node->lv == 6)
        // {
        //   printf("\n end of the cycle\n");
        //   // for(int hh = 0; hh < 8; hh++)
        //   // {
        //   //   printf("box_grid_idx = %d, pot = %.12f\n",box_grid_idx_vec[hh],(double)ptr_node->ptr_pot[box_grid_idx_vec[hh]]);
        //   // }
        // }
        
        //printf("\ni = %d\n", i);
        //GL_times[25] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;
        //GL_times[25] += omp_get_wtime() - aux_clock_omp;
        clock_gettime( CLOCK_REALTIME, &finish);
        GL_times[25] += ( finish.tv_sec - start.tv_sec ) + ( finish.tv_nsec - start.tv_nsec )/ 1000000000.;
        //clock_gettime( CLOCK_REALTIME, &GL_start);
        //printf("rsnew = %1.12e; rsold = %1.12e, alpha = %1.12e\n",rsnew,rsold,alpha);
        free(r);
        free(p);
        free(Ap);
        //printf("iters = %d\n",i);
        //clock_gettime( CLOCK_REALTIME, &GL_finish);
        //GL_times[50] += ( GL_finish.tv_sec - GL_start.tv_sec ) + ( GL_finish.tv_nsec - GL_start.tv_nsec )/ 1000000000.;
        return _SUCCESS_;
      }
      //GL_times[25] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;
      //GL_times[25] += omp_get_wtime() - aux_clock_omp;
      clock_gettime( CLOCK_REALTIME, &finish);
      GL_times[25] += ( finish.tv_sec - start.tv_sec ) + ( finish.tv_nsec - start.tv_nsec )/ 1000000000.;
    }
    // printf("\ni = %d, rsnew = %1.3Le, alpha = %1.3Le\n", i, rsnew,alpha);
    //  Updating p
    //clock_gettime( CLOCK_REALTIME, &GL_start);
        // printf("p[boxgrid] = %.12f, Ap = %.12f, pot = %.12f, r = %.12f\n",p[ptr_node->ptr_intr_box_grid_idx[5]] , Ap[ptr_node->ptr_intr_box_grid_idx[5]],
        //                                                     ptr_node->ptr_pot[ptr_node->ptr_intr_box_grid_idx[5]],r[ptr_node->ptr_intr_box_grid_idx[5]]);
    #pragma omp parallel for private(box_grid_idx)
    for (int j = 0; j < number_iter; j++)
    {
      box_grid_idx = ptr_node->ptr_intr_box_grid_idx[j];
      p[box_grid_idx] = r[box_grid_idx] + rsnew / rsold * p[box_grid_idx];
    }
    // clock_gettime( CLOCK_REALTIME, &GL_finish);
    // GL_times[50] += ( GL_finish.tv_sec - GL_start.tv_sec ) + ( GL_finish.tv_nsec - GL_start.tv_nsec )/ 1000000000.;
    
    // Updating rsold
    rsold = rsnew;
    
  }


  // printf("\n end of the cycle begins\n");
  // for(int hh = 0; hh < 8; hh++)
  // {
  //   printf("box_grid_idx = %d, pot = %.12f\n",box_grid_idx_vec[hh],(double)ptr_node->ptr_pot[box_grid_idx_vec[hh]]);
  // }

  
  free(r);
  free(p);
  free(Ap);
  
  return _FAILURE_;
}



int conjugate_gradient(struct node *ptr_node)
{



  int number_iter = ptr_node->grid_intr_size;

  struct timespec start, finish;

  int box_grid_idx;

  vtype H = 1.0 / (1 << ptr_node->lv);
  vtype one_over_H_pow_2 = 1.0 / (H * H);
  int size = (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1) * (ptr_node->box_real_dim_z + 1);

  vtype *r, *p, *Ap;
  r = (vtype *)calloc(size, sizeof(vtype));
  p = (vtype *)calloc(size , sizeof(vtype));
  Ap = (vtype *)malloc(size * sizeof(vtype));


  vtype alpha, rsold, rsnew;

  int grid_box_real_dim_X = ptr_node->box_real_dim_x + 1;
  int grid_box_real_dim_X_times_Y = (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1);

  rsold = 0.0;
  


  
  vtype alpha_aux;

  // int threads_to_use;
  // int auxxxx = number_iter / 250 > 1 ?  number_iter / 100 : 1; 
  // threads_to_use = NUMBER_OF_THEADS <  auxxxx ? NUMBER_OF_THEADS: auxxxx; 



  //#pragma omp parallel for private(box_grid_idx) reduction(+:rsold) num_threads(threads_to_use)
  #pragma omp parallel for private(box_grid_idx) reduction(+:rsold)
  for (int i = 0; i < number_iter; i++)
  {
    box_grid_idx = ptr_node->ptr_intr_box_grid_idx[i];
    // Computing r
    r[box_grid_idx] = ptr_node->ptr_d[box_grid_idx] + one_over_H_pow_2 *
                                                      (6.0 * ptr_node->ptr_pot[box_grid_idx] - ptr_node->ptr_pot[box_grid_idx + 1] - ptr_node->ptr_pot[box_grid_idx - 1] - ptr_node->ptr_pot[box_grid_idx + grid_box_real_dim_X] - ptr_node->ptr_pot[box_grid_idx - grid_box_real_dim_X] - ptr_node->ptr_pot[box_grid_idx + grid_box_real_dim_X_times_Y] - ptr_node->ptr_pot[box_grid_idx - grid_box_real_dim_X_times_Y]);
  
    // Copy r to p
    p[box_grid_idx] = r[box_grid_idx];
    // Computing rsold
    rsold += r[box_grid_idx] * r[box_grid_idx];
  }

  
  int cycle_size = size < _MAX_NUMBER_OF_ITERATIONS_IN_POISSON_EQUATION_ ? size : _MAX_NUMBER_OF_ITERATIONS_IN_POISSON_EQUATION_;
  // For cycle over mutually conjugate vectors p
  for (int i = 0; i < cycle_size; i++)
  {

    // Computing Ap and alpha
    alpha_aux = 0.0;
    rsnew = 0.0;




    // clock_gettime( CLOCK_REALTIME, &GL_start);
    //#pragma omp parallel shared(alpha) private(box_grid_idx)  num_threads(threads_to_use)
    #pragma omp parallel shared(alpha) private(box_grid_idx)
    {
      #pragma omp for reduction(+:alpha_aux) 
      for (int j = 0; j < number_iter; j++)
      {
        box_grid_idx = ptr_node->ptr_intr_box_grid_idx[j];
        Ap[box_grid_idx] = -one_over_H_pow_2 * (6.0 * p[box_grid_idx] - p[box_grid_idx + 1] - p[box_grid_idx - 1] - p[box_grid_idx + grid_box_real_dim_X] - p[box_grid_idx - grid_box_real_dim_X] - p[box_grid_idx + grid_box_real_dim_X_times_Y] - p[box_grid_idx - grid_box_real_dim_X_times_Y]);
        alpha_aux += p[box_grid_idx] * Ap[box_grid_idx];
      }

      if(omp_get_thread_num() == 0)
      {
        alpha = rsold / alpha_aux;
      }
      #pragma omp barrier

      #pragma omp for reduction(+:rsnew) 
      for (int j = 0; j < number_iter; j++)
      {
        box_grid_idx = ptr_node->ptr_intr_box_grid_idx[j];
        // New solution
        ptr_node->ptr_pot[box_grid_idx] += alpha * p[box_grid_idx];
        // New rest
        r[box_grid_idx] -= alpha * Ap[box_grid_idx];
        // // Computing rsnew
        rsnew += r[box_grid_idx] * r[box_grid_idx];
      }
    }
    // clock_gettime( CLOCK_REALTIME, &GL_finish);
    // GL_times[50] += ( GL_finish.tv_sec - GL_start.tv_sec ) + ( GL_finish.tv_nsec - GL_start.tv_nsec )/ 1000000000.;




    //* >> CHEKING ERROR SOLUTION CONDITION *//
    if ((i+1) % iter_between_check_potential_solution == 0)
    {

      clock_gettime( CLOCK_REALTIME, &start);
      if (poisson_error(ptr_node, r, rsnew, 1) == true)
      {
        clock_gettime( CLOCK_REALTIME, &finish);
        GL_times[25] += ( finish.tv_sec - start.tv_sec ) + ( finish.tv_nsec - start.tv_nsec )/ 1000000000.;

        free(r);
        free(p);
        free(Ap);
        return _SUCCESS_;
      }
      clock_gettime( CLOCK_REALTIME, &finish);
      GL_times[25] += ( finish.tv_sec - start.tv_sec ) + ( finish.tv_nsec - start.tv_nsec )/ 1000000000.;
    }

    
    //#pragma omp parallel for private(box_grid_idx) num_threads(threads_to_use)
    #pragma omp parallel for private(box_grid_idx)
    for (int j = 0; j < number_iter; j++)
    {
      box_grid_idx = ptr_node->ptr_intr_box_grid_idx[j];
      p[box_grid_idx] = r[box_grid_idx] + rsnew / rsold * p[box_grid_idx];
    }
    

    rsold = rsnew;
    
  }
  
  free(r);
  free(p);
  free(Ap);

  return _FAILURE_;
}


// int conjugate_gradient_S(struct node *ptr_node)
// {

//   int number_iter = ptr_node->grid_intr_size;



//   //int number_of_threads = 8;

//   //omp_set_num_threads(NUMBER_OF_THEADS);

//   //clock_t aux_clock;
//   //double aux_clock_omp;
//   struct timespec start, finish;

//   int box_grid_idx;

//   vtype H = 1.0L / (1 << ptr_node->lv);
//   vtype one_over_H_pow_2 = 1.0L / (H * H);
//   int size = (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1) * (ptr_node->box_real_dim_z + 1);

//   vtype *r, *p, *Ap, *Ar;
//   r = (vtype *)calloc(size, sizeof(vtype));
//   p = (vtype *)calloc(size , sizeof(vtype));
//   Ap = (vtype *)malloc(size * sizeof(vtype));
//   Ar = (vtype *)malloc(size * sizeof(vtype));

//   vtype alpha, beta, rsold, rsnew,Ar_times_r;

//   int grid_box_real_dim_X = ptr_node->box_real_dim_x + 1;
//   int grid_box_real_dim_X_times_Y = (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1);

//   rsold = 0;
//   #pragma omp parallel for private(box_grid_idx) reduction(+:rsold)
//   for (int i = 0; i < number_iter; i++)
//   {
//     //printf("my id = %d\n", omp_get_thread_num());
//     box_grid_idx = ptr_node->ptr_intr_box_grid_idx[i];
//     // Computing r
//     r[box_grid_idx] = ptr_node->ptr_d[box_grid_idx] + one_over_H_pow_2 *
//                                                           (6.0 * ptr_node->ptr_pot[box_grid_idx] - ptr_node->ptr_pot[box_grid_idx + 1] - ptr_node->ptr_pot[box_grid_idx - 1] - ptr_node->ptr_pot[box_grid_idx + grid_box_real_dim_X] - ptr_node->ptr_pot[box_grid_idx - grid_box_real_dim_X] - ptr_node->ptr_pot[box_grid_idx + grid_box_real_dim_X_times_Y] - ptr_node->ptr_pot[box_grid_idx - grid_box_real_dim_X_times_Y]);
  
//     // Copy r to p
//     //p[box_grid_idx] = r[box_grid_idx];
//     // Computing rsold
//     rsold += r[box_grid_idx] * r[box_grid_idx];
//   }




//   alpha = 0;
//   #pragma omp parallel for private(box_grid_idx) reduction(+:alpha)
//   for (int i = 0; i < number_iter; i++)
//   {
//     //printf("my id = %d\n", omp_get_thread_num());
//     box_grid_idx = ptr_node->ptr_intr_box_grid_idx[i];
//     // Computing Ar0
//     Ar[box_grid_idx] = - one_over_H_pow_2 * (6.0 * r[box_grid_idx] - r[box_grid_idx + 1] - r[box_grid_idx - 1] - r[box_grid_idx + grid_box_real_dim_X] - r[box_grid_idx - grid_box_real_dim_X] - r[box_grid_idx + grid_box_real_dim_X_times_Y] - r[box_grid_idx - grid_box_real_dim_X_times_Y]);
//     // Computing alpha0 = <r0,r0> / <r0,Ar0>
//     alpha += r[box_grid_idx] * Ar[box_grid_idx];
//   }

//   // Computing alpha0
//   alpha = rsold / alpha;

//   //printf("alpha0 =  %.12f\n",alpha);

//   beta = 0;
//   int cycle_size = size < _MAX_NUMBER_OF_ITERATIONS_IN_POISSON_EQUATION_ ? size : _MAX_NUMBER_OF_ITERATIONS_IN_POISSON_EQUATION_;
//   // For cycle over mutually conjugate vectors p
//   for (int i = 0; i < cycle_size; i++)
//   {
    
//     rsnew = 0;
//     #pragma omp parallel for private(box_grid_idx) reduction(+:rsnew) 
//     for (int j = 0; j < number_iter; j++)
//     {
//       box_grid_idx = ptr_node->ptr_intr_box_grid_idx[j];
//       //Computing p
//       p[box_grid_idx] = r[box_grid_idx] + beta * p[box_grid_idx];
//       //Computing Ap
//       Ap[box_grid_idx] = Ar[box_grid_idx] + beta * Ap[box_grid_idx];
//       //Computing new solution
//       //printf("p[box_grid_idx] , r[box_grid_idx] = %f ; %f \n",p[box_grid_idx],r[box_grid_idx]);
//       ptr_node->ptr_pot[box_grid_idx] += alpha * p[box_grid_idx];
//       //Computing r
//       r[box_grid_idx] -=  alpha * Ap[box_grid_idx];
//       //Computing r^2
//       rsnew += r[box_grid_idx] * r[box_grid_idx];
//     }

    

//     //* >> CHEKING ERROR SOLUTION CONDITION *//
//     if ((i+1) % iter_between_check_potential_solution == 0)
//     {
//       //printf("\ni = %d\n", i);
//       //aux_clock = clock();
//       //aux_clock_omp = omp_get_wtime();
//       clock_gettime( CLOCK_REALTIME, &start);
//       if (poisson_error(ptr_node, r, rsnew, 1) == true)
//       {
//         //printf("\ni = %d\n", i);
//         //GL_times[25] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;
//         //GL_times[25] += omp_get_wtime() - aux_clock_omp;
//         clock_gettime( CLOCK_REALTIME, &finish);
//         GL_times[25] += ( finish.tv_sec - start.tv_sec ) + ( finish.tv_nsec - start.tv_nsec )/ 1000000000.;
//         //clock_gettime( CLOCK_REALTIME, &GL_start);
//         free(r);
//         free(p);
//         free(Ap);
//         free(Ar);
//         //clock_gettime( CLOCK_REALTIME, &GL_finish);
//         //GL_times[50] += ( GL_finish.tv_sec - GL_start.tv_sec ) + ( GL_finish.tv_nsec - GL_start.tv_nsec )/ 1000000000.;
//         return _SUCCESS_;
//       }
//       //GL_times[25] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;
//       //GL_times[25] += omp_get_wtime() - aux_clock_omp;
//       clock_gettime( CLOCK_REALTIME, &finish);
//       GL_times[25] += ( finish.tv_sec - start.tv_sec ) + ( finish.tv_nsec - start.tv_nsec )/ 1000000000.;
//     }

//     Ar_times_r = 0;
//     #pragma omp parallel for private(box_grid_idx) reduction(+:Ar_times_r) 
//     for (int j = 0; j < number_iter; j++)
//     {
//       box_grid_idx = ptr_node->ptr_intr_box_grid_idx[j];
//       //Computing Ar
//       Ar[box_grid_idx] = - one_over_H_pow_2 * (6.0 * r[box_grid_idx] - r[box_grid_idx + 1] - r[box_grid_idx - 1] - r[box_grid_idx + grid_box_real_dim_X] - r[box_grid_idx - grid_box_real_dim_X] - r[box_grid_idx + grid_box_real_dim_X_times_Y] - r[box_grid_idx - grid_box_real_dim_X_times_Y]);
//       //Computing beta = <Ar_{k+1},r_{k+1}>/<r_k,r_k>
//       Ar_times_r += Ar[box_grid_idx] * r[box_grid_idx];
//     }

//     //Computing beta = <Ar_{k+1},r_{k+1}>/<r_k,r_k>
//     beta = rsnew / rsold;
//     //computing alpha = <r_{k+1},r_{k+1}>/(<Ar_{k+1},r_{k+1}> - (beta/alpha)<r_k,r_k>)
//     alpha = rsnew / (Ar_times_r - beta/alpha * rsnew);

//     // Updating rsold
//     rsold = rsnew;
    
//     //printf("hola = %d, rsnew = %f, alpha = %.12f, beta = %.12f\n",i,rsnew,alpha,beta);
//     // printf("p[boxgrid] = %.12f, Ap = %.12f, pot = %.12f, r = %.12f\n",p[ptr_node->ptr_intr_box_grid_idx[5]] , Ap[ptr_node->ptr_intr_box_grid_idx[5]],
//     //                                                         ptr_node->ptr_pot[ptr_node->ptr_intr_box_grid_idx[5]],r[ptr_node->ptr_intr_box_grid_idx[5]]);

//   }

  
//   free(r);
//   free(p);
//   free(Ap);
//   free(Ar);
  
//   return _FAILURE_;
// }
