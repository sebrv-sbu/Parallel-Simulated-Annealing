#include <stdio.h>
/*********************************************
* Tour cost calculates the cost of the       *
* entire tour for all city pairs.  Any tour  *
* is passed in and its cost returned.        *
**********************************************/
double TourCost(unsigned short* tour_pointer);

/***************************************************************************
* ReadTSP  reads in the data from the input file from TSPLIB either:       *
*              distances stored in lower traingular form   OR              *
*              converts to it lower triangular form                        *  
*              (coordinate data is not used).                              * 
****************************************************************************/
void ReadIsing (FILE *infile);

/*** FindSection: This function finds a given section of the input file & **
 *                returns a pointer positioned to the first record of that *
 *                section; section titles should be passed without the pre-*
 *                ceding '$'; if it can't find the right section, the func-*
 *                tion returns NULL                                        *
 ***************************************************************************/

FILE *FindSection(FILE *fp, char *input_section);


/* Reads starting points from files - Dec 7 2025*/
/* SRV */
int *ReadStartingPoints(FILE *start_file);

void InitIsing(FILE *infile);

double UnifStartingPoints();

double RandStartingPoints();

double SumWeights();
