/*
 * space_check.c
 *
 * Allocate or reallocate variable/s if it is necessary
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

#include "space_check.h"

int space_check(int *ptr_cap, int size, float increase, const char *format, ...)
{
  if (*ptr_cap < size)
  {
    int cap = *ptr_cap;

    if (*format != 'p')
    {
      printf("Error in the realocation of the pointer\n");
      printf("The format string must starts with the letter p\n");
      return _FAILURE_;
    }

    ++format; // Next character in the string array format

    if (*format < '0' || *format > '9')
    {
      printf("Error in the realocation of the pointer\n");
      printf("The format string must have a digit in its second letter\n");
      return _FAILURE_;
    }
    int cntr = *format - '0'; // Number of variables inserted in the space checking function

    *ptr_cap = cntr > 0 ? ((int)(increase * size)) : *ptr_cap; // Duplicated capacity

    ++format; // Type of variable: i = int, v = vtype

    va_list args;
    va_start(args, format);

    char format_aux;

    for (int i = 0; i < cntr; i++)
    {
      format_aux = *format;
      ++format; // Type of pointer: 1 = *, 2 = **
      if (format_aux == 'i' && *format == '1')
      {
        int **pptr_int = NULL;
        pptr_int = va_arg(args, int **);
        if (*pptr_int == NULL)
        {
          *pptr_int = (int *)malloc((*ptr_cap) * sizeof(int));
        }
        else
        {
          int *ptr_aux = NULL;
          ptr_aux = (int *)realloc(*pptr_int, (*ptr_cap) * sizeof(int));
          if (ptr_aux == NULL)
          {
            printf("Error in the realocation of ptr_int\n");
            return _FAILURE_;
          }
          else
          {
            *pptr_int = ptr_aux;
          }
          ptr_aux = NULL;
        }
        for (int j = cap; j < *ptr_cap; j++)
        {
          (*pptr_int)[j] = 0;
        }
      }
      else if (format_aux == 'i' && *format == '2')
      {
        int ***ppptr_int = NULL;
        ppptr_int = va_arg(args, int ***);
        if (*ppptr_int == NULL)
        {
          *ppptr_int = (int **)malloc((*ptr_cap) * sizeof(int *));
        }
        else
        {
          int **pptr_aux = NULL;
          pptr_aux = (int **)realloc(*ppptr_int, (*ptr_cap) * sizeof(int *));
          if (pptr_aux == NULL)
          {
            printf("Error in the realocation of ppptr_int\n");
            return _FAILURE_;
          }
          else
          {
            *ppptr_int = pptr_aux;
          }
          pptr_aux = NULL;
        }
        for (int j = cap; j < *ptr_cap; j++)
        {
          (*ppptr_int)[j] = NULL;
        }
      }
      else if (format_aux == 'v' && *format == '1')
      {
        vtype **pptr_vtype = NULL;
        pptr_vtype = va_arg(args, vtype **);
        if (*pptr_vtype == NULL)
        {
          *pptr_vtype = (vtype *)malloc((*ptr_cap) * sizeof(vtype));
        }
        else
        {
          vtype *ptr_aux = NULL;
          ptr_aux = (vtype *)realloc(*pptr_vtype, (*ptr_cap) * sizeof(vtype));
          if (ptr_aux == NULL)
          {
            printf("Error in the realocation of ptr_vtype\n");
            return _FAILURE_;
          }
          else
          {
            *pptr_vtype = ptr_aux;
          }
          ptr_aux = NULL;
        }
        for (int j = cap; j < *ptr_cap; j++)
        {
          (*pptr_vtype)[j] = 0.0;
        }
      }
      else if (format_aux == 'v' && *format == '2')
      {
        vtype ***ppptr_vtype = NULL;
        ppptr_vtype = va_arg(args, vtype ***);
        if (*ppptr_vtype == NULL)
        {
          *ppptr_vtype = (vtype **)malloc((*ptr_cap) * sizeof(vtype *));
        }
        else
        {
          vtype **pptr_aux = NULL;
          pptr_aux = (vtype **)realloc(*ppptr_vtype, (*ptr_cap) * sizeof(vtype *));
          if (pptr_aux == NULL)
          {
            printf("Error in the realocation of pptr_vtype\n");
            return _FAILURE_;
          }
          else
          {
            *ppptr_vtype = pptr_aux;
          }
          pptr_aux = NULL;
        }
        for (int j = cap; j < *ptr_cap; j++)
        {
          (*ppptr_vtype)[j] = NULL;
        }
      }
      else if (format_aux == 'n' && *format == '1')
      {
        struct node **pptr_node = NULL;
        pptr_node = va_arg(args, struct node **);
        if (*pptr_node == NULL)
        {
          *pptr_node = (struct node *)malloc((*ptr_cap) * sizeof(struct node));
        }
        else
        {
          struct node *ptr_aux = NULL;
          ptr_aux = (struct node *)realloc(*pptr_node, (*ptr_cap) * sizeof(struct node));
          if (ptr_aux == NULL)
          {
            printf("Error in the realocation of ptr_node\n");
            return _FAILURE_;
          }
          else
          {
            *pptr_node = ptr_aux;
          }
          ptr_aux = NULL;
        }
      }
      else if (format_aux == 'n' && *format == '2')
      {
        struct node ***ppptr_node = NULL;
        ppptr_node = va_arg(args, struct node ***);
        if (*ppptr_node == NULL)
        {
          *ppptr_node = (struct node **)malloc((*ptr_cap) * sizeof(struct node *));
        }
        else
        {
          struct node **pptr_aux = NULL;
          pptr_aux = (struct node **)realloc((*ppptr_node), (*ptr_cap) * sizeof(struct node *));
          if (pptr_aux == NULL)
          {
            printf("Error in the realocation of pptr_node\n");
            return _FAILURE_;
          }
          else
          {
            *ppptr_node = pptr_aux;
          }
          pptr_aux = NULL;
        }
        for (int j = cap; j < *ptr_cap; j++)
        {
          (*ppptr_node)[j] = NULL;
        }
      }
      else if (format_aux == 'c' && *format == '1')
      {
        struct cell_struct **pptr_cell = NULL;
        pptr_cell = va_arg(args, struct cell_struct **);
        if (*pptr_cell == NULL)
        {
          *pptr_cell = (struct cell_struct *)malloc((*ptr_cap) * sizeof(struct cell_struct));
        }
        else
        {
          struct cell_struct *ptr_aux = NULL;
          ptr_aux = (struct cell_struct *)realloc(*pptr_cell, (*ptr_cap) * sizeof(struct cell_struct));
          if (ptr_aux == NULL)
          {
            printf("Error in the realocation of ptr_cell\n");
            return _FAILURE_;
          }
          else
          {
            *pptr_cell = ptr_aux;
          }
          ptr_aux = NULL;
        }
        for (int j = cap; j < *ptr_cap; j++)
        {
          initialize_cell_struct(&((*pptr_cell)[j]));
        }
      }
      else if (format_aux == 'b' && *format == '1')
      {
        bool **pptr_bool = NULL;
        pptr_bool = va_arg(args, bool **);
        if (*pptr_bool == NULL)
        {
          *pptr_bool = (bool *)malloc((*ptr_cap) * sizeof(bool));
        }
        else
        {
          bool *ptr_aux = NULL;
          ptr_aux = (bool *)realloc(*pptr_bool, (*ptr_cap) * sizeof(bool));
          if (ptr_aux == NULL)
          {
            printf("Error in the realocation of ptr_bool\n");
            return _FAILURE_;
          }
          else
          {
            *pptr_bool = ptr_aux;
          }
          ptr_aux = NULL;
        }
        for (int j = cap; j < *ptr_cap; j++)
        {
          (*pptr_bool)[j] = false;
        }
      }
      else
      {
        printf("Error in the realocation of the pointer\n");
        printf("The format string must have a convination of i or v, plus a number at position = %d \n", cntr);
        return _FAILURE_;
      }
      ++format; // Moving to the next type of variable
    }

    va_end(args);
  }

  return _SUCCESS_;
}
