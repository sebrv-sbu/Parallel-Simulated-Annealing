/*******************************************************************************
* I am calc_ave_error_bar.c and I alone become the calc_ave_error_bar
* executable.  
*
*1)    Usage: calc_ave_error_bar filename.data 
*      Input: probname.lambda.data 
*      Output: one line file named lambda.ave
*
* input file should be problemname.lambda.data and have all the results 
* from every lambda' runs in it.  
* output file gets overwritten with every call to this program.
*
**Version 3.0 data FORMAT****LG:9-02********************************
* 8 col format:  no longer  just 80 col                            *
* moved Filename after Energy and Iterations                       *
* Version Seed Lambda Gain Distribution Energy Iterations Filename *
********************************************************************
*
*    Count does not matter: may be 100 or 1000 runs--the number is counted here  
* 
*2)    At completion,  a 1 line output file is generated lambda.ave 
*      has the following format:
*
*  mean cost | mean iter | 
*              min cost | max cost | min iter | max iter
*
* The last 4 numbers will be the 4 corners of the error bars.
*
* I also used erand48 for recreate-ablility
*
* 3)    These lambda.ave files are read by the AVERAGE_results_LINUX3.0 script
*       and  each one is over-written with every call to this program.  
*       This program 
*       is called once for every problemname.lambda.data file.    
*       lambda.ave is deleted at completion of the script.
*
********************************************************************************/
#ifdef ICC
#include <mathimf.h>
#else
#include <math.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define LOWBITS 0x330E    /* for drand (32 bits) to erand (48 bits) compatibility */
#define BYTESIZE 8

long idum;

int main (int argc, char ** argv)
{  /* begin main */
int randin(int min, int max);
double erand48(unsigned short int xsubi[3]);
FILE            *input,*output;
int             i, j,m,nn;
double          energy[2000], itnum[2000], temp_energy,temp_itnum;
double          version, lambda, gain; 
long            seed;
double          dist; 
char            filename[60]; 
double          min_mean_energy,max_mean_energy,min_mean_itnum,max_mean_itnum;
double          mean_energy,mean_itnum,sum_energy,sum_itnum;
char            fname2[80]="lambda.ave"; /*one line output file*/
char            s[180];
unsigned short xsubj [3]; /* xsubj gets changed with each call to erand48 */ 
long seedval;   
unsigned short left16, middle16;
int left;

	if (argc != 2)
	{
	   printf ("Usage:  calc_ave_error_bar <infile> \n");
	   return 1;
	}
        input = fopen(argv[1],"r");
        if (input == NULL)
              {
                printf ("TSP ERROR_BAR: Error in opening input file\n");
                exit (1);
               }
        output = fopen(fname2,"w");
        if (output == NULL)
              {
                printf ("TSP ERROR_BAR: Error in opening output file\n");
                exit (1);
              }

/* initialize erand48 seedval    */
  seedval  = 987557231 ;
  xsubj[0] = LOWBITS;
  middle16 = (unsigned short) seedval ;
  xsubj[1] = middle16;
  left = seedval >> (BYTESIZE * sizeof(unsigned short));
  left16 = (unsigned short) left;
  xsubj[2] = left16;
 
                 m = 0;
                 sum_energy = sum_itnum = 0.0;

/****************************************************************
* This part reads the input file into energy and itnum
****************************************************************/
               while(NULL!=fgets(s,120,input)) 
               {
                 m++;             
                 sscanf(s,"%lf %ld %lf %lf %lf %lf %lf %s \n", 
                        &version, &seed, &lambda, &gain, &dist, 
                        &energy[m],&itnum[m], filename);                 

/*********************************************************************
**************** comment out print statement, just in case************
**********************************************************************/

/*		 printf("%lf %ld %lf %lf %lf %lf %lf %s\n", 
 *                     version, seed, lambda, gain, dist,
 *		     energy[m],itnum[m], filename);                 
 **********************************************************************/


                 sum_energy += energy[m];
                 sum_itnum += itnum[m];
/*	        printf("%ld %lf %lf \n", m, sum_itnum, itnum[m]);*/
              }
/****************************************************************
* Now calculate mean energy and mean itnum
****************************************************************/
               mean_energy = sum_energy / m;
               mean_itnum = sum_itnum / m;

/***************************************************************
* Now initialize the min_mean and max_mean "yets" to the true means. 
****************************************************************/
               min_mean_energy = max_mean_energy = mean_energy;
               min_mean_itnum = max_mean_itnum = mean_itnum; 

/****************************************************************
* Bootstrap algorithm:
* Randomly create synthetic sample data sets of the same size as
* the original data set.  In this case 'm' is number of data points
* and we are generating j=100 synthetic sample data sets of size m.  
****************************************************************/
 
              for (j=1; j<=100; j++) 
		{	/* begin bootstrap for j */
                 sum_energy = sum_itnum = 0.0;

/*******************************************************************
* randomly pick with replacement from the original data set.  The
* replicate or synthetic data set will contain 1/e or approximately 37% 
* duplicate points.
*********************************************************************/
                 for (i=1; i<=m; i++) 
                  {    /* begin for i = m sample size */
           	      nn= (int) (m * erand48(xsubj)); 
                      if (nn == 0) nn=1;
         
                    sum_energy += energy[nn];
                    sum_itnum += itnum[nn];
                   }   /* end for i = m sample size */
          
           
/**************************************************************
* find the mean energy and itnum for this synthetic sample set
***************************************************************/
                 temp_energy = sum_energy / m;
                 temp_itnum = sum_itnum / m;

/**************************************************************
* find the extremal mean energy and itnum for all sample sets
***************************************************************/
                 if ( temp_energy < min_mean_energy ) min_mean_energy = temp_energy;
                 if ( temp_energy > max_mean_energy ) max_mean_energy = temp_energy;
                 if ( temp_itnum < min_mean_itnum ) min_mean_itnum = temp_itnum;
                 if ( temp_itnum > max_mean_itnum ) max_mean_itnum = temp_itnum;
                }	/* end bootstrap for j */


              fprintf(output,"%g %g %g %g %g %g\n",
		     mean_energy,mean_itnum,
              min_mean_energy,max_mean_energy,min_mean_itnum,max_mean_itnum);

              printf("%g %g %g %g  %g %g\n",
		      mean_energy,mean_itnum,
              min_mean_energy,max_mean_energy,min_mean_itnum,max_mean_itnum);


              fclose(output);
return 0;
} /* end main */











