/* Seb RV Jul 31 2024
 * Starting edge weight from scratch. */

#include "edge_wt.h"
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <sa.h>
#ifndef MPI_INCLUDED
#include <mpi.h>
#include "MPI.h"
#endif

double GetDistance(coord A, coord B){
  /* Pythagoras */
  return sqrt(pow(A.coord_x-B.coord_x,2)+pow(A.coord_y-B.coord_y,2));
}

int GetJthNearestNeighbour(int home, int j){
  
  return sorted_lists[(long)home*(long)(ncities-1)+(long)j-1];
}


void MergeSort(int home, int *final_result){
  int city;
  coord nhome=node_coords[home];
  int cross_city=0;
  for (city=0; city<ncities; city++){
    if(city != home){
    distances[city-cross_city].distance = GetDistance(node_coords[city], nhome);
    distances[city-cross_city].city     = city;
    }
    else
      cross_city++;
  }
  int i;
  int k;
  dist_and_city *_;

  for (i=1; i<ncities-1; i*=2){
    for (k=0; k<ncities-1; k+=i*2){
       if (ncities-1<k+i*2)
         MergeLim(k, i);
       else
         MergeInc(k, i);                                                                                                               
      }
    _ = temp;
    temp = distances;
    distances = _;
  }
    for(i=0;i<ncities-1;i++) 
     final_result[i]=distances[i].city;
}  

void MergeLim(int left, int increment){
  int ind_left=left;
  int ind_right=left+increment;
  int ind_temp=left;
  int i=0;
  /* if outside <= ind_right then all we would do is copy the array.
   * Thus, we include this step */
  if(ncities-1>ind_right){
    int size=ncities-1-left;
  /* This is all written in an extremely odd way to reduce branching as much as
   * possible */
      while ((ind_left<left+increment) && (ind_right<ncities-1)){ /* i<outside-left is guaranteed by 
                                                                   the two other conditions. */
       if(distances[ind_left].distance < distances[ind_right].distance){
         temp[ind_temp] = distances[ind_left];
         ind_left++;
         }
       else{
         temp[ind_temp] = distances[ind_right];
         ind_right++;
       }
       ind_temp++;
       i++;
     }
     if (ind_left==left+increment){/* all the left has been copied */
       while(i<size){
         temp[ind_temp] = distances[ind_right];
         ind_right++;
         ind_temp++;
         i++;
       }
     }
   else{ /* all the right has been copied */
     while (i<size){
       temp[ind_temp] = distances[ind_left];
       ind_left++;
       ind_temp++;
       i++;
       }
   }
   return;
  }
  else{ /* this means that outside <= ind_right and we just copy the array */
   while (ind_left<ncities-1){
       temp[ind_temp] = distances[ind_left];
       ind_left++;
       ind_temp++;
       i++;
       }
   }
  return;
}

void MergeInc(int left, int increment){
  int size=increment*2;
  int ind_left=left;
  int ind_right=left+increment;
  int ind_temp=left;
  int i=0;
  int outside=left+2*increment;
  while((ind_left-left<increment) && (ind_right<outside) && (i<size)){
    if(distances[ind_left].distance < distances[ind_right].distance){
      temp[ind_temp] = distances[ind_left];
      ind_left++;
    }
    else{
      temp[ind_temp] = distances[ind_right];
      ind_right++;
    }
    ind_temp++;
    i++;
  }
if (i==size)
  return;
if(ind_left-left==increment){/*this means all of the left has been copied */
  while(i<size){
    temp[ind_temp] = distances[ind_right];
    ind_right++;
    ind_temp++;
    i++;
    }
  return;
}
else{ /* all the right has been copied */
  while(i<size){
    temp[ind_temp] = distances[ind_left];
    ind_left++;
    ind_temp++;
    i++;
    }
  }
  return;
}

void FreeSortedLists(void)
{
  free(sorted_lists);
}
