#include "edge_wt.h"
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <global.h>
#include "error.h"
#include <stdbool.h>
#ifndef MPI_INCLUDED
#include <mpi.h>
#include "MPI.h"
#endif

/*** STATIC VARIABLES ******************************************************/



/*********************************************
* Tour cost calculates the cost of the       *
* entire tour for all city pairs.  Any tour  *
* is passed in and its cost returned.        *
**********************************************/

double tour_cost(unsigned short*tour_pointer)
{  /* begin tour_cost*/
int i; /*loop counter */
 double cost = 0;  /*cost of tour*/
for ( i=0; i<ncities-1; i++)
  {/*begin for */
  cost += GetDistance(node_coords[tour_pointer[i]], node_coords[tour_pointer[i+1]]) ;
  } /*end for */
cost += GetDistance(node_coords[tour_pointer[0]],node_coords[tour_pointer[ncities-1]]) ;
 return cost;
}  /* end tour_cost*/

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

/* ReadTSP - July 31 2024. Will be completely changing ReadTSP to use      *
 * coordinate format. - Seb RV.                                            */
void ReadTSP (FILE *infile) { /* begin ReadTSP instance file */
  char DimensionString [80];
  int loop2= 0, loop1 = 0;
  bool ReadingNodeCoordData = false;
  bool ReadingCoordinateData = false;
  /*   double floatdata = 0.0; memory hog */
  float x,y;
  int i = 0;
  int _index;
  char stringdata[180];

/* heading names for TSP LIB files */
#define NAME_HEADER 						"NAME : "
#define TYPE_HEADER						"TYPE : "
#define COMMENT_HEADER					"COMMENT : "
#define DIMENSION_HEADER				"DIMENSION : "
#define EDGE_WEIGHT_TYPE_HEADER 		"EDGE_WEIGHT_TYPE : "
#define EDGE_WEIGHT_FORMAT_HEADER   "EDGE_WEIGHT_FORMAT : "
#define DISPLAY_DATA_TYPE_HEADER    "DISPLAY_DATA_TYPE : "
#define EDGE_WEIGHT_SECTION_HEADER  "EDGE_WEIGHT_SECTION"
#define NODE_COORD_SECTION_HEADER   "NODE_COORD_SECTION"

  if( !infile )
    error("ReadTSP: could not locate input file\n");

   stringdata[0]=0;
   while ((strstr(stringdata, "EOF")== NULL) && (strstr(stringdata, "$$$")==NULL))
   {
      if ((ReadingNodeCoordData == false) && (ReadingCoordinateData == false))
      { /* reading text parts of file into string data such as */
         /* the Name, Comments, Dimension */
      fgets (stringdata,180, infile);
      } /* reading text parts of file */

/* Don't use coordinate data--just read and print to output file as strings*/
/* can do same as "REWD" stuff to use it in future. */

      else if (ReadingNodeCoordData == true)
      {  /* reading in edge weight data as numeric values */
         /* for edge weights in lower triangular form */
         for (i = 0; i < ncities; ++i)
         { 
           fgets(stringdata,180,infile);
           sscanf(stringdata, "%d %f %f", &_index, &x, &y);
           node_coords[i].coord_x=x;
           node_coords[i].coord_y=y;
            
         }   /*  endfor edge weights in lower triangular form */
         ReadingNodeCoordData = false;
      }  /* reading in edge weight data as numeric values */

      if (strstr (stringdata, DIMENSION_HEADER) != NULL)
      {/* get the dimension */
        loop1 = 0;
        /* extract the dimension value from the string*/

         for (loop2 = strlen (DIMENSION_HEADER); loop2< strlen(stringdata); loop2++)
         { /* pick off the problem dimension from the input string */
          DimensionString[loop1] = stringdata[loop2];
                loop1++;
         }/* done picking off the problem dimension from the input string */

         DimensionString[loop1] = '\0';  /* null character */

          /* convert the problem dimension to integer */
         ncities = atoi (DimensionString);

         node_coords=calloc(ncities, sizeof(coord));
         if(node_coords==NULL){
          error("tsp_sa: node_coords could not be allocated");
         } 
      }
      if (strstr (stringdata, NODE_COORD_SECTION_HEADER) != NULL)
      {  /* we're up to the edge weight section */
         /* we need to empty the string inorder to keep */
         /* us from doing this again */
         /* if not we will do the comparison and it will */
         /* be true forever!! */
        stringdata[0] = '\0';
        ReadingNodeCoordData= true;
      }  /* done with edges */

      } /* end while (infile != EOF $$$ ) */
}  /* end of ReadTSP */

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
