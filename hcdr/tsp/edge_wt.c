/* Seb RV Jul 31 2024
 * Starting edge weight from scratch. */

#include "edge_wt.h"
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <sys/param.h>


double GetDistance(coord A, coord B){
  /* Pythagoras */
  return sqrt(pow(A.coord_x-B.coord_x,2)+pow(A.coord_y-B.coord_y,2));
}


int GetJthNearestNeighbour(unsigned short home, unsigned short n){
  coord nhome = node_coords[home];
  for(int i=0; i<ncities; i++){
    distances[i].city     = i;
    distances[i].distance = GetDistance(node_coords[i], nhome);
    /* if we want to get into this, might be worth trying different
     * method here of just updating distance --- will look into
     * this but I'm not sure it will actually be faster... */
  }
  FloydRivestSelect( 0, ncities-1, n);
  return distances[n].city;
}

#define FRCONSTANT1 600
#define FRCONSTANT2 0.5
#define FRCONSTANT3 0.5
void FloydRivestSelect(unsigned short left, unsigned short right, unsigned short k){
  unsigned short i;
  unsigned short n;
  unsigned short j;
  unsigned short pivot_index;
  double pivot;
  double sd;
  double z;
  double s;
  while (right>left){
  /*  if(right-left>600){
       n = right-left + 1;
      i = k-left +1;
      z=log(n);
      s=FRCONSTANT2*exp(2*z/3);
      if (i>n/2)
        sd = FRCONSTANT3*sqrt(z*s*(n-s)/n);
      else
        sd = -FRCONSTANT3*sqrt(z*s*(n-s)/n);
      left  = MAX(left, (int)(ceil(k-i*s/n+sd)));
      right = MIN(right, (int)floor(k+(n-i)*s/n+sd));
      FloydRivestSelect(left, right, k);
      return;
    }*/
    /*Will fix later*/

      pivot_index=(left+right)/2;
      pivot = distances[pivot_index].distance;
      i = left;
      j = right;
      swap_inplace(left, pivot_index);
      if (distances[right].distance > pivot)
        swap_inplace(right, left);
      while (i <= j){
        while(distances[i].distance<pivot)
          i++;
 	while(distances[j].distance>pivot)
          j--;
 	if (i <= j){
          swap_inplace(i,j);
 	  i++;
 	  j--;
 	}
      }     
      if (k <= j)       
        right=j;     
      else if (k >= i) /* Necessary in case n=k */
	left=i;     
      else  return;
  }
}

void swap_inplace(unsigned short left, unsigned short right){
  dist_and_city _  = distances[left];
  distances[left]  = distances[right];
  distances[right] = _;
}
