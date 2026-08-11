#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <global.h>
#include "error.h"
#include <stdbool.h>
#include "move.h"
#include "random.h"
#include "initialize.h"

/*** STATIC VARIABLES ******************************************************/



/*********************************************
* Tour cost calculates the cost of the       *
* entire tour for all city pairs.  Any tour  *
* is passed in and its cost returned.        *
**********************************************/


/***************************************************************************
* ReadTSP  reads in the data from the input file from TSPLIB either:       *
*              distances stored in lower traingular form   OR              *
*              converts to it lower triangular form                        *  
*              (coordinate data is not used).                              * 
****************************************************************************
* Formerly:input_lib routine Created 10-99 by: Lorraine Greenwald          *
* LG 03-03 use two files params and instance. Needed for BIG problems      *
* LG 08-02 all writing to output file done by WriteResults and FinalMove   *
* LG: 08-01 do not write problem instance everytime to output file.        *
*           input file template has the distances, eliminates redundancy   *
*  LG: 07-05-00 set criterion to 0 in TSP input file to "hardwire" frozen  *
* for this application                                                     *
* Assumptions:the in file are already open                                 *
****************************************************************************/
void ReadIsing (FILE *infile) { /* begin ReadTSP instance file */
  int _n_points;
  int _n_edges;
  int site_index, neighbour_index;
  double weight;
  int num_used;

  fscanf(infile, "%d %d", &_n_points, &_n_edges);
  /* For now we assume regular toroidal */
  InitConfig(_n_points, _n_edges);

  bool *used = calloc(sizeof(*used), n_points*degree);
  if (used == NULL){ printf("Malloc error"); exit(1); }

  char line_data[180];
  line_data[0] = 1;
  while (fgets(line_data, sizeof(line_data), infile) != NULL){
    if ( sscanf(line_data, "%d %d %lf", 
      &site_index, &neighbour_index, &weight) != 3
      ) continue;
    site_index --;
    neighbour_index -- ;

    for (num_used = 0; used[site_index*degree + num_used]; ++num_used);

    VERTEX_EDGE(site_index, num_used) = 
      (edge_weight){.neighbour = neighbour_index, .weight = weight};

    used[site_index*degree+num_used]=true;
    
    for (num_used = 0; used[neighbour_index*degree + num_used]; ++num_used);

    VERTEX_EDGE(neighbour_index, num_used) =
      (edge_weight){.neighbour = site_index, .weight = weight};
    used[neighbour_index*degree+num_used]=true;
  }
  /**************************************************************
   * Structure:                                                 *
   * [spin][edge][edge][edge]...[spin][edge][edge][edge]...     *
   * giving 1+degree elements per vertex. move.h contains       *
   * relevant macros.                                           *
   **************************************************************/
  free(used);
}

/*** FindSection: This function finds a given section of the input file & **
 *                returns a pointer positioned to the first record of that *
 *                section. Section titles should be passed without the pre-*
 *                ceding '$'. If it can't find the right section, the func-*
 *                tion returns NULL.                                       *
 ***************************************************************************/

FILE *FindSection(FILE *fp, char *input_section)
{
  int      c;                      /* input happens character by character */
  int      nsought;                       /* holds length of section title */
  char     *base;                              /* string for section title */
  long     looksite;                            /* file position indicator */
  
  rewind(fp);                /* start looking at the beginning of the file */
  
  nsought = strlen(input_section);
  base = (char *)calloc(MAX_RECORD, sizeof(char));
  
/*** while loop goes through the file character by character... ************/
   while ( (c=getc(fp)) != EOF) {               /* ...until if finds a '$' */
    
    if ( c == '$') {                  /* found a sectioning control string */
      looksite = ftell(fp);                               /* where are we? */
      base = fgets(base,MAX_RECORD,fp);  /* get full sect title (without $)*/
      if ( !(strncmp(base, input_section, nsought)) ) {
        fseek(fp, looksite, 0);        /* found the sought string: reposi- */
        fscanf(fp, "%*s\n");              /* tion the pointer to '$', then */
	free(base);
        return(fp);                    /* advance the pointer to after the */
      }                                     /* section title and return it */
      else {                          /* didn't find it: skip this control */
        fseek(fp, looksite, 0);                    /* record, keep looking */
        fscanf(fp, "%*s");                 /* NOTE: "%*s" advances pointer */
      }                                              /* without assignment */
    }
  }
   
  free(base);
 
  return(NULL);                         /* couldn't find the right section */
}

int *ReadStartingPoints(FILE *start_file){
  /* Currently Non-Functional */
  return 0;
} 

double RandStartingPoints(){
  for (int i = 0 ; i < n_points; ++i){
    VERTEX_SPIN(i) = RandomInt(2) ? 1:-1;
  }
  curr_cost = CalcCost();
  curr_cost += SumWeights();
  return curr_cost;
}

double UnifStartingPoints(){
  for (int i = 0 ; i < n_points; ++i){
    VERTEX_SPIN(i) = 1; 
  }
  curr_cost = CalcCost();
  curr_cost += SumWeights();
  return curr_cost;
}

void InitIsing(FILE *infile) {
   long pid;            /* processor id for to seed drand48 (not really used)*/

   /* initialize random number generator drand48 with pid as seed */
   pid=getpid();  /* get proc id  */
   srand48(pid);  /* seed drand48 */

   /* read in the tsp data from file */

   ReadIsing(infile);
}

double SumWeights(){
  double sum=0;
  for (int i = 0 ; i < n_points; ++i){
    for (int j = 0 ; j < degree; ++j){
      sum += fabs(VERTEX_EDGE(i,j).weight);
    }
  }
  return sum/2;
}
