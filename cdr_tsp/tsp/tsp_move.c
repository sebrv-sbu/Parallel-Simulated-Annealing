
/*************************************************
* tsp_move.c  Created 4-01 by: Lorraine Greenwald
* Handles tour, allocates and frees memory 
* for tour arrays, calculates cost for entire tour
***************************************************/
/* Edited greatly by Seb */
#ifdef ICC
#include <mathimf.h>
#else
#include <math.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include "tsp_move.h"  /* need prototypes */ 
#include "edge_wt.h" /* needed for distance functions */

int *curr_tour = NULL;  /* pointer to tour of city ids in current order visited*/ 
int *curr_position = NULL;  /*maintains the current tour position for each city*/
                    /*index is city id-1 element is curr_tour index */ 
double curr_cost = 0.0;  /*cost of the current tour*/ 

/***********************************************
* allocates dynamic memory for the tour arrays.*
* Tours are 1 by dimension data structures.    *
* called by input_lib in tsp_lib.c              *
************************************************/

int tour_allocate(void)
{  /* begin tour_allocate */
   curr_tour = (int*) malloc (sizeof (int) * ncities);
   if (curr_tour == NULL)
      return 1;

   curr_position = (int*) malloc (sizeof (int) * ncities);
   if (curr_position == NULL)
      return 1;

   return 0;
}  /* end tour_allocate */

/**************************************************
* deallocates dymanic memory for the tour arrays  *
* called by main                 *
****************************************************/

void tour_deallocate(void)
{  /* begin tour_deallocate */

   free (curr_position);

   free (curr_tour);
}  /* end tour_deallocate */


/***************************************
* in _cost initializes tour costs to 0 *
****************************************/
/* why does this exist? */

/*********************************************
* Tour cost calculates the cost of the       *
* entire tour for all city pairs.  Any tour  *
* is passed in and its cost returned.        *
**********************************************/

double tour_cost(int *tour_pointer)
{  /* begin tour_cost*/
  /* Edited Aug 2, 2024 by Seb RV. Now it uses the GetDistance
   * function on the list of node coordinates. */
  int i; /*loop counter */
 double cost = 0;  /*cost of tour*/
for ( i=0; i<ncities-1; i++)
  {/*begin for */
	  cost += GetDistance(node_coords[tour_pointer[i]], node_coords[tour_pointer[i+1]]) ;
  } /*end for */
cost += GetDistance(node_coords[tour_pointer[0]],node_coords[tour_pointer[ncities-1]]) ;
 return cost;
}  /* end tour_cost*/


/********************************************
* Start Tour                                *
* come up with a starting tour by mindlessly*
* assigning cities is numeric order         *
********************************************/
double start_tour( void )
{  /* begin start_tour*/
int i; /*loop counter */
/* double  *cost_pointer;   dummy pointer for array parameter */
//init_cost(); /* why is this here? What problem is this solving? */
/* curr_tour contains city ids */
/* curr_position contains index into tour*/
/* city id is one more than index to start out */

/* changing this because it doesn't make sense to have   *
 * curr_tour[i]=i+1 at any point other than in printing? *
 * Seb RV Aug 2 2024 */
for ( i=0; i<ncities; i++)
  {/*begin for */
    curr_tour [i] = i;       /* curr_tour contains city_id */
    curr_position [i] = i;   /* curr_position contains an index (not city_id) */

 } /*end for */

/* let's get the cost */

curr_cost = tour_cost(curr_tour);



return curr_cost;

}  /* end start_tour*/


/*********************************************
* Print the current tour and its cost.       *
*********************************************/

void print_curr_tour(void)
{  /* begin print current tour */
int i;
printf("\n current cost is: %f\n", curr_cost);  
printf ("current tour is:\n");  

for ( i=0; i<ncities ; i++)
  {
          printf("%d\t", curr_tour[i]+1);  
  }  /*end for print */
printf ("\n");  
}  /* end print current tour */





