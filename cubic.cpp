/* File: cubic_extrema.c -- Contains FindCubicExtrema() and supporting routines
   that implement the cubic spline extrema algorithm. Given a set of x,y data
   points, determine in the cubic spline sense the relative extrema of the
   function describing the data. -- by Mike J. Courtney
*/
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include "cubic.h"
/*  Primary routine that implements the cubic spline extrema algorithm. Calls
   ComputeSecDerivs() to compute the second derivatives, computes quadratic
   coefficients, calls FindQuadRoots() to solve quadratic roots, determines if
   roots are valid abscissa, and calls ComputeY() to compute ordinates.
*/
BOOL FindCubicExtrema (
  unsigned int num_pnts,    /* input - address of number of x & y points */
  float *x,                 /* input - array of x values */
  float *y,                 /* input - array of y values */
  struct point *extr)       /* output - singly linked list of extrema */
{
  float a, b, c;        /* coefficients of quadratic equation */
  float x1, x2;         /* roots of quadratic equation to be computed */
  unsigned int i;       /* array index */
  float *sec_deriv;     /* array of second derivatives of each data interval */
  BOOL root_stat;       /* computation status returned by FindRoots() */
  BOOL valid_flag;      /* TRUE if at least one valid root found */
  BOOL first_flag;      /* TRUE if current root is the first one */

  /* allocate array for second derivatives */
  //sec_deriv = malloc(num_pnts * sizeof(float));
  /* compute the second derivatives */
  ComputeSecDerivs(num_pnts, x, y, sec_deriv);
  /* initialize extrema flags */
  valid_flag = FALSE;
  first_flag = TRUE;
  /* loop through all the input points and find the extrema */
  for (i = 0; i < num_pnts - 1; i++)
    {
    /* compute the quadratic coefficients */
    a = 3.0 * (sec_deriv[i+1] - sec_deriv[i]);
    b = 6.0 * (x[i+1] * sec_deriv[i] - x[i] * sec_deriv[i+1]);
    c = 6.0 * (y[i+1] - y[i]) - (2.0 * x[i+1] * x[i+1] - x[i] * x[i] + 2 *
        x[i] * x[i+1]) * sec_deriv[i];
    c -= (x[i+1] * x[i+1] - 2.0 * x[i] * x[i+1] - 2.0 * x[i] * x[i]) *
        sec_deriv[i+1];
    /* determine the roots of the cubic extrema quadratic equation */
    root_stat = FindQuadRoots(a, b, c, &x1, &x2);
    if (root_stat != FAILURE)
      {
      /* if root x1 was calculated 1fully */
      if (root_stat != FAILURE1)
        {
        /* Determine if root is within the interval */
        if ((x1 > x[i]) && (x1 < x[i+1]))
          {
          /* first root (extremum) */
          if (first_flag == TRUE)
            first_flag = FALSE;
          /* beyond first valid root so allocate next extremum structure */
          else
            {
            extr->next = (struct point *)malloc(sizeof(struct point));
            extr = extr->next;
            }
          extr->next = NULL;
          extr->x = x1;
          /* compute the corresponding value of y at the extreme x value */
          ComputeY(i, x, extr->x, y, sec_deriv, &extr->y);
          valid_flag = TRUE;
          }
        }
     /* if root x2 was calculated 1fully */
     if (root_stat != FAILURE2)
        {
        /* Determine if root is within the current interval */
        if ((x2 > x[i]) && (x2 < x[i+1]))
          {
          /* first root (extremum) */
          if (first_flag == TRUE)
            first_flag = FALSE;
          /* beyond first valid root so allocate next extremum structure */
          else
            {
            extr->next = (struct point *)malloc(sizeof(struct point));
            extr = extr->next;
            }
          extr->next = NULL;
          extr->x = x2;
          /* compute the corresponding value of y at the extreme x value */
          ComputeY(i, x, extr->x, y, sec_deriv, &extr->y);
          valid_flag = TRUE;
          }
        }
      } /* end of if(root_stat ! = FAILURE) */
    } /* end of for(i) */
  free(sec_deriv);
  if (valid_flag == TRUE)
    return 1;
  else
    {
    /*  Set next to NULL just in case it was not set in the loop - this is
        so that free loop will operate properly upon return  */
    extr->next = NULL;
    return FAILURE;
    }
  }
/* Use input x,y data to form tridiagonal matrix and compute second
   derivatives of function in the cubic spline sense. */
BOOL ComputeSecDerivs (
  unsigned int num_pnts,    /* input - number of x & y points */
  float *x,                 /* input - array of x values */
  float *y,                 /* input - array of y values */
  float *sec_deriv)         /* output - array of 2nd derivatives of intervals */
{
  unsigned int   i; /* index */
  float ftemp;      /* temporary float */
  float main_diag[num_pnts -2]; /* ptr to matrix main diagonal array */
  float *diag;      /* ptr to matrix diagonal array */
  float *right;     /* ptr to array of right sides of matrix equations */
  //main_diag = malloc((num_pnts - 2) * sizeof(float));
  //diag = malloc((num_pnts - 2) * sizeof(float));
  //right = malloc((num_pnts - 2) * sizeof(float));

  /* compute the matrix main and off-diagonal values */
  /* even though the calling program is suppose to have guaranteed that the
     x values are increasing, assert that neither of the diagonal
     differences are zero to avoid a divide by zero condition */
  for (i = 1; i < num_pnts - 3; i++)
    {
    main_diag[i-1] = 2.0 * (x[i+1] - x[i-1]);
    std::cout << x[i] << "\n";
    //assert(main_diag[i-1] > 0);
    }
  for (i = 0; i < num_pnts - 1; i++)
    {
    diag[i-1] = x[i+1] - x[i];
    assert(diag[i-1] > 0);
    }
  /* compute right hand side of equation */
  for (i = 1; i < num_pnts - 1; i++)
    right[i-1] = 6.0 * ((y[i+1]-y[i])/diag[i-1]-(y[i]-y[i-1])/diag[i-2]);
  /* forward eliminate tridiagonal */
  sec_deriv[0] = 0.0;
  sec_deriv[num_pnts - 1] = 0.0;
  for (i = 1; i < num_pnts - 2; i++)
    {
    ftemp = diag[i] / main_diag[i];
    right[i] -= (right[i-1] * ftemp);
    main_diag[i] -= (diag[i-1] * ftemp);
    }
  /* backward substitution to solve for second derivative at each knot */
  for (i = num_pnts - 2; i > 0; i--)
    sec_deriv[i] = (right[i-1] - diag[i-1] * sec_deriv[i+1]) / main_diag[i-1];
  free(main_diag);
  free(diag);
  free(right);
  return 1;
}
/*  Solve for roots x1 and x2 of a quadratic equation of the form
    a * (x * x) + b * x + c = 0 using the following formula x1 = d / a  and
    x2 = c / d, where d = -0.5 * [b + sgn(b) * sqrt(b*b - 4ac)].
   This algorithm is particularly good at yielding accurate results when
   a and/or c are small values.
*/
BOOL FindQuadRoots (
  float a,      /* input - coefficient a of quadratic equation */
  float b,      /* input - coefficient b of quadratic equation */
  float c,      /* input - coefficient c of quadratic equation */
  float *x1,    /* output - first root computed */
  float *x2)    /* output - second root computed */
{
  float d;      /* root algorithm variable */
  BOOL  root_stat;  /* status of root computations */

  d = b * b - 4 * a *c;
  if (d < 0)
    return FAILURE;
  else
    {
    d = (float)sqrt((double)d);
    /* make the result of sqrt the sign of b */
    if (b < 0 )
      d = -d;
    d = -0.5 * (b + d);
    /* solve for the roots of the quadratic */
    /* if both root computations will yield divide by zero ... forget it! */
    if ( (a == 0) && (d == 0) )
      return FAILURE;

    root_stat = 1;
    /* compute first root if denominator a is not zero */
    if (a == 0)
      root_stat = FAILURE1;
    else
      *x1 = d / a;
    /* compute second root if denominator d is not zero */
    if (d == 0)
      root_stat = FAILURE2;
    else
      *x2 = c / d;
    return root_stat;
    }
  }
/*  Given an abscissa (x) location, computes the corresponding cubic spline
   ordinate (y) value.
*/
void ComputeY (
  unsigned int i,   /* input - array index */
  float *x,         /* input - array of x values */
  float x_value,    /* input - x value at which to solve for y */
  float *y,         /* input - array of y values */
  float *sec_deriv, /* input - array of second derivatives of each data interval */
  float *y_value)   /* output - address of y extreme value at x */
{
  float A, B, C, D; /* cubic spline coefficients */
  float ftemp;      /* temporary float */
  /* compute the standard cubic spline coefficients */
  A = (x[i + 1] - x_value) / (x[i + 1] - x[i]);
  B = 1 - A;
  ftemp = (float) pow((double)(x[i + 1] - x[i]), 2.0) / 6.0;
  C = (A * A * A - A) * ftemp;
  D = (B * B * B - B) * ftemp;
  /* compute the ordinate value at the abscissa location */
  *y_value = A * y[i] + B * y[i + 1] + C * sec_deriv[i] + D *
            sec_deriv[i + 1];
  return;
}
