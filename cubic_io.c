/* File: cubic_io.c -- Input/output program for the FindCubicExtrema() routine.
   Reads an ascii file of x,y data, calls the routine, and lists the
   returned extrema. -- by Mike J. Courtney
*/
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "cubic.h"
main()
{
  /* the following are required for main() */
  unsigned int i;           /* array indices */
  unsigned int num_extr;    /* number of extrema found */
  BYTE len;                 /* input data parsing variable */
  char data_file[80];       /* ascii data file name */
  char str[40], line[80];   /* input data parsing variables */
  FILE *fp;                 /* ascii data file pointer */
  struct point *first_extr; /* pointer to the first extreme structure */
  /* the following are required for FindCubicExtrema() */
  unsigned int num_pnts;    /* number of x,y pairs */
  float *x_in, *y_in;       /* arrays of x and y input data values */
  struct point *extr;       /* pointer to current extreme structure */
  /* open input data file and determine number of input points */
  printf("\n Enter name of input file: ");
  gets(data_file);
  if((fp = fopen(data_file, "r")) != NULL)
    {
    num_pnts = 0;
    while(!feof(fp))
      {
      fgets(line, 80, fp);
      num_pnts++;
      }
    num_pnts -= 1;
    printf("\n The number of input points was %d.",num_pnts);

    /* allocate the input data arrays */
    x_in   = malloc(num_pnts * sizeof(float));
    y_in   = malloc(num_pnts * sizeof(float));

    /* read in the each data line and parse into x and y arrays */
    rewind(fp);
    for (i = 0; i < num_pnts; i++)
      {
      /* get a line */
      fgets(line, 80, fp);
      len = strcspn(line,",");
      /* get the x value */
      strncpy(str, line, len);
      /* NULL terminate the x value */
      str[len] = '\0';
      x_in[i] = atof(str);
      /* get the y value */
      strcpy(str, line+len+1);
      y_in[i] = atof(str);
      }
    fclose(fp);
    /* allocate first structure of linked list of output extrema */
    extr = (struct point *)malloc(sizeof(struct point));
    /* save the address of the first structure in the list */
    first_extr = extr;
    /* call the routine that computes the extrema */
    if (FindCubicExtrema(num_pnts, x_in, y_in, extr) == FAILURE)
      printf("\n\7 No extrema found !\n");
    else
      {
      /* print the linked list of extrema */
      printf("\n\n Relative extrema computed:");
      extr = first_extr;
      num_extr = 0;
      while (extr)
        {
        printf("\n  %d.  y = %f at x = %f", num_extr+1, extr->y, extr->x);
        extr = extr->next;
        num_extr++;
        }
      }
    free(x_in);
    free(y_in);
    /* free the linked list of extreme structures */
    do
      {
      /* point to first structure */
      extr = first_extr;
      /* save address of next extreme structure so that when we free the
         current one, we still have a pointer to the next one */
      first_extr = extr->next;
      free(extr);
      }
    while (first_extr != NULL);
    }
  else
    printf("\n Couldn't open %s !", data_file);
}
