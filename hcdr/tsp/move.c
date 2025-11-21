/******************************************************
* move.c  Created 10-99 by: Lorraine Greenwald        * 
* LG 03-03: use unsigned shorts for tour arrays       *
*           and floats rather than doubles for edges  * 
*           for  BIG tsp instances memory crunch      *
* LG 03-03: use two files for BIG tsp instance        *
* LG 08-02: Integrate lsa.c 9.3.1 into tsp code       *
* LG 08-01: add distributions for move generation     *
* LG 02-02: add pareto dist for move generation       *
* LG 03-02: add poisson dist for move generation      *
* LG 05-02: fix dist=7 to calc abs value 5-16-02      *
* LG 05-02: add generate_dev call for move generation *
*    visiting distribution control parameter          *
*    read from input file to DistP.q distributions.h  *
* Handles tours (current, new, global minimum)        *
* allocates and frees memory for tour arrays          *
* calculates cost for entire tour and altered tour    *
* generates, accepts, reject moves                    *
*******************************************************/

#ifdef ICC
#include <mathimf.h>
#else
#include <math.h>
#endif
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>

#include "error.h"
#include "move.h"   /* need acc_tab (AccStats) & ap (AParms)struct and prototypes */
                    /* also need pi for distributions  LG 03-02 */
#include "random.h"
#include "sa.h"                                        /* *ONLY* for flags */
#include "edge_wt.h" /* included for prototypes and neighbors*/
#include "distributions.h"   /* problem independent distributions */
#include "initialize.h"

#include <mpi.h>                     /* this is the official MPI interface */
#include "MPI.h"  /* our own structs and such only needed by parallel code */

static int tour_debug=0;  /*debug flag for tours and edges */




/*** STATIC VARIABLES ******************************************************/

static AParms    ap;                /* static copy of annealing parameters */

/* SRV 2025-11-16 made acc_tab no longer a pointer and made it into
 * the actual struct itself */
static AccStats  acc_tab;   /* struct to accumulate acceptance statistics */

static unsigned  int    nhits;    /* number of moves since the start of execution */
static unsigned int    nsweeps;     /* number of sweeps since start of execution */


static int       prolix;           /* flag for printing stat info (prolix) */
                                   /* value of zero means no trace */
static char      *prolixfile;                   /* filename of prolix file */

static MPI_Aint nbytes;
static MPI_Aint size_arr;

/*** TSP TOUR VARIABLES ******************************************************/

/* Edited by seb: Changed shorts to ints because some of the test problems are massive */
static unsigned short *curr_tour_and_pos = NULL;

static unsigned short *curr_tour = NULL;  /* pointer to tour of city ids 
                                  in current order visited*/ 
static unsigned short *curr_position = NULL;  /* maintains the current tour position 
                               for each city*/
                    /*index is city id-1 element is curr_tour index */ 
static unsigned short swap[2];  /* array of indices of swapped elements in tour*/ 
static int *min_tour = NULL;  /* pointer to the minimum annealing tour 
                                 found so far-not really SA methodology*/ 
static double curr_cost = 0.0;  /*cost of the current tour*/ 
static double new_cost = 0.0;  /* cost of the proposed new tour*/ 
static double min_cost = 0.0;  /* cost of the minimum tour found so far 
                                  for fun-not really SA methodology*/
static double tour_max = 0.0;  /* max index into tour array */
static double neighbor_max = 0.0; /* max number for furthest neighbor 
                                    this is the max move size */
static int neighbor_maxint = 0; /* integer value max number for 
                                   furthest neighbor (move size) */


/*** FUNCTIONS *************************************************************/


/*** INITIALIZING AND RESTORING FUNCTIONS **********************************/


/***************************************************************************
* InitTSP                                                                  *
*       no longer returns starting energy *p_chisq=InitTSP(infile)         *
*         done by calling StartTour from InitialMove in tsp_sa.c           *
*       initializes the random number generator                            *
*       reads in the data from the input file ReadTSP in tsp_sa.c          *
*       calls start_tour to get initial energy and starting tour           *
****************************************************************************/

/* SEB July 31 2024: Changing InitTSP to reflect coordinate representation of tsp. 
 * Coordinate representation will allow for less memory consumption which 
 * will help greatly for larger problems. */
void InitTSP (FILE *infile) { /* begin InitTSP */


   long pid;            /* processor id for to seed drand48 (not really used)*/

   /* initialize random number generator drand48 with pid as seed */
   pid=getpid();  /* get proc id */
   srand48(pid);  /* seed drand48 ( TSP uses erand48) */

   /* read in the tsp data from file */

   ReadTSP(infile);
 
}  /* end of InitTSP */



/***********************************************
* allocates dynamic memory for the tour arrays.*
* Tours are 1 by dimension data structures.    *
* called by input_lib in tsp_sa.c              *
************************************************/
 int tour_allocate(void)
{  /* begin tour_allocate */
    curr_tour_and_pos=(unsigned short*) malloc (sizeof (unsigned short) * 2* ncities);   
if (curr_tour_and_pos == NULL)
         return 1;
    curr_tour = curr_tour_and_pos;
    curr_position=curr_tour_and_pos+ncities;
   if (curr_position == NULL)
      return 1;
 if (tour_debug) 
{ /* begin debug */
      min_tour = ( unsigned short*) malloc (sizeof (unsigned short) * ncities);
   if (min_tour == NULL)
      return 1;
} /* end debug */
    return 0;
}  /* end tour_allocate */


/*******************************************
* come up with a starting tour and set min *
* mindlessly assign cities is numeric order*
* returns starting energy (cost)           *
********************************************/
double StartTour( void )
{  /* begin start_tour*/

int i; /*loop counter */
FILE *debugptr = NULL; /* pointer for tour debug file */

/*************************************************/
/* curr_tour contains city ids                   */
/* curr_position contains index into tour        */
/* city id is one more than index to start out   */
/*************************************************/
/* again, why is this function here? (Aug 2 2024 SRV) */
//init_cost();
if (tour_allocate()==1)
  error("Error in allocating memory for tour");

for ( i=0; i<ncities; i++)
  {/*begin for */
    curr_tour [i] = i;     /* curr_tour contains city_id */
                             /* index is one less than city_id */
    curr_position [i] = i;   /* curr_position contains an index (not city_id) */
  } /*end for */

 /* let's get the cost */
 curr_cost = tour_cost(curr_tour);

 if (tour_debug) 
   { /* begin debug */
     /* first tour is min tour found so far only for fun, so not for everytime*/
     min_cost = curr_cost ;
     for ( i=0; i<ncities; i++)
       {/*begin for */
           min_tour [i] = curr_tour [i];
       } /*end for */
    } /* end debug */

    /* set up max for random swap and neighbor generators */
    tour_max = (double) (ncities-1);      /* max index into tour array */
    neighbor_max = (double) (ncities-2);  /* for now use all neighbors */
    neighbor_maxint = ncities -2;  /* for now use all neighbors */

  if (tour_debug) 
    { /* begin debug */
    debugptr = fopen("tour_debug", "w");
    if ( !debugptr ) {
      perror("tsp_sa");
      exit(1);
    }
  /* Seb here. Removed some debug statements that no longer make sense
   * when using node coordinates. */

	/* let's print tour we have */
	fprintf(debugptr,"\n The Starting Tour and cost is:");
	print_curr_tour (debugptr);
        fclose(debugptr);
    } /* end debug */

 return curr_cost;
 }  /* end start_tour*/


/****************************************
* init_cost initializes tour costs to 0 *
*****************************************/
 void init_cost(void)
{  /* begin init_cost*/
   min_cost = 0;
   curr_cost = 0;
   new_cost = 0;
   swap[0] = swap[1] = 0;
}  /* end init_cost*/


/***************************************************************************
 *** InitMoves: initializes the following moves.c-specific stuff:      *****
 *              - static annealing parameter struct (ap)                   *
 *              - initializes random number generator in lsa.c             *
 *              - initializes acc_tab for acceptance statistics            *
 *              - initializes distances for nmin                           *
 *                                                                         *
 *              it then returns the initial temperature to the caller      *
 ***************************************************************************/

double InitMoves(FILE *fp, int Tau)
{  /* begin init moves */
   unsigned short *xsubj; 
   long seedval;
   unsigned short left16, middle16;
   int left;

      xsubj = (unsigned short *)calloc(3, sizeof(unsigned short));

/* read annealing paramters */

  ap           = ReadAParameters(fp);       /* ap: static annealing params */
  if(ap.interval % Tau || !ap.interval){
    error("Reading Annealing Parameters: interval not divisible by tau");
  }
  ap.max_count = 0; 

  if ( equil == 1 )   /* read equilibration params and put them into lsa.c */
      InitEquilibrate(fp);
 
/* initialze the random number generator, now erand48() */

#ifdef MPI
  seedval  = ap.seed + myid;       /* each processor gets a different seed */
#else
  seedval  = ap.seed;
#endif

  xsubj[0] = LOWBITS;

  middle16 = (unsigned short)seedval;
  xsubj[1] = middle16;

  left     = seedval >> (BYTESIZE * sizeof(unsigned short));
  left16   = (unsigned short)left;
  xsubj[2] = left16;

  InitERand(xsubj);        /* makes the xsubj array static to random.c */

/* acc_tab is for statistics like acceptance ratio etc. */

    distances = (dist_and_city *)calloc(ncities, sizeof(dist_and_city));
  acc_tab.theta_bar = (floor (ncities/2));
  /* THETA_INIT varies by tsp instance so cannot be a constant */
  /* so do not use: acc_tab->theta_bar = THETA_INIT; */
  acc_tab.hits      = 0;
  acc_tab.success   = 0;
 
  nhits   = 0;
  nsweeps = 0;

  /* Initialize some byte info that we will need later 
   * for state messages */
  size_arr=ncities*2*sizeof(unsigned short);
  nbytes = 2*sizeof(unsigned char) + sizeof(unsigned short) +
  size_arr+ 2*sizeof(unsigned int) + 2*sizeof(double);
/* Finally, return the start temperature. */
  return ap.start_tempr;
}  /* end init moves */

/*** RestoreMoves: restores move generator from state file *****************
 *           NOTE: InitMoves will be called before this function during    *
 *                 a restore                                               *
 ***************************************************************************/
 /* SRV 2025-11-16 - Not currently working */
void RestoreMoves(MoveState *moveptr){
//  int i;
//  ncities   = moveptr->ncities;
//  nhits     = moveptr->nhits;
//  nsweeps   = moveptr->nsweeps;
//  curr_cost = moveptr->curr_cost;
//
//  tour_max        = (double) (ncities-1);
//  neighbor_max    = (double) (ncities-2);  /* for now use all neighbors */
//  neighbor_maxint = ncities -2;  /* for now use all neighbors */

 // if ( tour_allocate() == 1 )
 //  error("Error in allocating memory for tour\n");

 // for(i=0; i<ncities; i++){
 //   curr_tour[i]=moveptr->curr_tour[i];
 //   curr_position[i]=moveptr->curr_position[i];
 // }

 // free(moveptr->curr_tour);
 // free(moveptr->curr_position);
 // free(acc_tab);

 // acc_tab = moveptr->acc_tab_ptr;

 // free(moveptr);
}
/* RestoreProlix: not functional.Dummy for now.
 * */
void RestoreProlix(void)
{
  return;
}

/*** FUNCTIONs FOR FINALIZING A RUN ***************************************/

/*** GetFinalInfo: collects stop energy and final count for output to the **
 *                 data file                                               *
 ***************************************************************************/

AParms GetFinalInfo(void)
{
  ap.stop_energy = curr_cost;
  ap.max_count   = nhits;

  return ap;
}

 /**************************************************
* deallocates dynamic memory for the tour arrays   *
* called by final move in tsp_s.c                  *
****************************************************/
 void tour_deallocate(void)
{  /* begin tour_deallocate */
if (tour_debug) 
{ /* begin debug */
   free (min_tour);
} /* end debug */
    free (curr_tour_and_pos);
}  /* end tour_deallocate */


/*** MOVE GENERATION - PART 1: FUNCS NEEDED IN LSA.C ***********************/

/* GenerateMove: wrapper for Move makes a move and returns difference of    *
 * energies before and after the move          *
 ***************************************************************************/

/**************************************************
* generate_move provides next move to try         *
* LG: 05-02 used generate_dev routine from        * 
* distributions.c to be consistent                *
* use 2-opt swap for new tour generation          *
* uniform random picks first index of tour to swap*
* use various distributions to pick second index  *
* check that both indices are unique and          *
* that they are not identical to previous pair    *
* return the difference in cost curr_cost-new_cost*
**************************************************/

double GenerateMove()
{  /* begin GenerateMove*/
  /* delta_e is the change in energy returned as function value */
 /* it is the difference in the cost of the current tour and the new tour*/
  double delta_e = 0.0;
  double xi;  /* uniform random variable for move control */
  double theta; /* move control variable to pick neighbors*/
  int i = 0;
  int j = 0;
  unsigned short city_id;  /*city id for swap */

/* make a move, get energy and return delta_e */



/* increase counters */
/* SRV Nov 19 2025 - switching nhits to be after nsweeps *
 * in order to force UpdateControl to only occur after   *
 * communication steps.                                  */
  acc_tab.hits++;

  nhits++;
  nsweeps=  (nhits) * lam_group_size;   /* calculate number of sweeps */
/* update statistics if interval passed & at least one sweep completed */

/* Seb RV August 1st 2024: Changed this to happen after nsweeps is 
 * calculated. It is what happens in the CDR code so, idk I don't really
 * think it makes a difference. */

 /* Removed all of this logic entirely. UpdateControl can only occur 
 during an update_stats step. Seb RV Nov 19 2025 */


/* now generate new micro state */
/* ThermoDynamics! (/s) Seb RV 2024 */
/*************************************************************************
* Let try different distributions to see impact on move generation        *
* distribution type is in DistP.distribution from distributions.h         *
* Need to pick i uniformly (always). Pick neighbor j using distributions: *
* once i and j are picked swap array is set up the same (see below)       *
***************************************************************************/

  /*pick first tour index randomly */
  i=(int) RandomInt(tour_max);  /* uniform dist i=[0, prob_dimension-1] */

  /* control the neighbor pick the lam way */
  theta = generate_dev(acc_tab.theta_bar, DistP.distribution, DistP.q); 

  if (DistP.distribution == 7) /* tsp needs a positive value for theta */
   {theta = fabs(theta);}

  if (theta > neighbor_max)  /* in this case Lam used uniform dist*/
    {xi=RandomReal();
     theta = xi * neighbor_max;} 

  /*****************************************************************/
  /* round robin method is another way tried once upon a time      */
  /*      if (theta > neighbor_max)                                */
  /*         { theta  = theta % neighbor_max; }  * modulo division */
  /*****************************************************************/

  /* pick a neighbor index;  0th neighbor is self so don't go there */
   j=  (int)(ceil(theta));

  swap[0]=i;   /* i is index into tour (just dandy) */
  /* j= index into the row of neighbor; neighbor.to_city has 
   *    value of the city_id -1 which is what we need */

  city_id = GetJthNearestNeighbour(curr_tour[i],j);  
/***********************************************************************
 *** city_id variable is actually one less because it is the address   *
 *** the address is what we need to index position (SRV fixed this     *
 *** Aug 23 2024                                                       *
************************************************************************/
    swap[1] = curr_position[city_id];
 /* let's get the cost */
 new_cost = curr_cost;
 new_cost = calc_new_cost(curr_tour, swap, new_cost);
 delta_e=new_cost - curr_cost; 

 return delta_e;

 }  /* end generate_move*/


/*** AcceptMove: sets new energy cost and tour in current cost and tour    * 
 *               for the next step and keeps track of the number of        *
 *               successful moves for acceptance statistics                *
 ***************************************************************************/

void AcceptMove(void)
{  /* begin accept_move*/
  int hold;                            /* temporary swap location */
 /* actually swap the elements in curr tour for real */
    /* new cost must be saved in curr cost */
     curr_cost = new_cost ; /* my way */
  /* elements must be swapped in curr tour and position */ 
      /* element of curr_tour is city ID */
      /* element of swap is index to curr_tour */
      /* element of curr_position is index to curr_tour */
      /* index to curr_postion is city ID -1 */
      hold = curr_tour [ swap [0]];
      curr_tour [swap[0] ] = curr_tour [swap [1]];
      curr_tour [swap[1] ] = hold;
       /* ** hold is a city id so must sub 1 to get index into curr_position** */
      curr_position [hold ] = swap [1];
      curr_position [curr_tour [swap[0]] ] = swap[0];
    /* count the acceptances we have */
       acc_tab.success++;
 if (tour_debug) 
{ /* begin debug */
   /* just for grins keep the best found so far */
      global_min ();
} /* end debug */


 }  /* end accept_move*/


/*******************************************
* Reject move does nothing right now.      *
********************************************/
 void RejectMove()
{  /* begin reject_move*/
 
  /* sorry, nothing to do. I leave original state intact until accept.   */
 }  /* end reject_move*/


/*** MOVE GENERATION - PART 2: A FUNC NEEDED IN MOVE.C (BUT NOT LSA.C) *****/

/*** UpdateControl: each 'interval' number of steps, acceptance stats are **
 *                  updated here; acceptance statistics are collected for: *
 *                  TSP:                                                   *
 *                  each tour generated is a step. Simplest case;          *
 *                  LJ:                                                    *
 *                  all particles in the same struct, since all parameters *
 *                  are of the same order of magnitude in this problem;    *
 *                  FLY gene circuit model:                                *
 *                  more complicated cost functions need move control for  *
 *                  each parameter individually;                           *
 *                  this function prints prolix stuff, if required (-p)    *
 ***************************************************************************/

void UpdateControl(uint16_t *m_success) 
{

  
/* Seb RV Nov 20 2025 */
/* Now UpdateControl has to manually check for itself.    * 
 *                                                         *
 * UpdateControl gets called every proc_tau moves, which   *
 * is a multiple of tau, so even if nsweeps=nhits this can *
 * can only be called during an UpdateStats step.           */
if (myid == 0){
}
  if (!(nsweeps % ap.interval) ){

  FILE       *prolixptr;                            /* prolix file pointer */

/* open prolix file for appending new move stats */
  if ( myid == 0 ) {
    if ( prolix ) {
      prolixptr = fopen(prolixfile, "a");
      if ( !prolixptr ) {
	perror("UpdateControl");
	exit(1);
      }
    }
  }

/* if parallel, pool the accpetance statistics */
/* Or how about we don't add two unneccessary, *
 * expensive MPI_Allreduces.                   */


/* calculate acceptance ratio (for all parameters) and adjust theta_bar    *
 * according to the gain (see King-Wai's thesis, p. 23)                    *
 *                                                                         *
 * Removed uneccessary acc_tab->acc_ratio computations                     *
 * I did this way earlier but I didn't date it - SRV NOV 20 2025           */

  /* tsp no log - we do not need x */
  acc_tab.theta_bar+= ap.gain_div_interval * (double)((int) *m_success-44);
  /* for all the trouble this stupid fucking function gave us
   * in trying to not make assumptions about sweeps, they
   * use this magic number that assumes we have 1 sweep every 100
   * moves. Whatever. SRV Nov 19 2025*/

    if ( acc_tab.theta_bar > neighbor_max ) {
      acc_tab.theta_bar = neighbor_max;
    }
    else if ( acc_tab.theta_bar < THETA_MIN ) {
      acc_tab.theta_bar = THETA_MIN;
    }
 
/* if -p: root node prints prolix information to prolix file */

  if ( prolix ) {
    if ( myid == 0 ) { 
      fprintf(prolixptr, "nsteps = %8d bar = %10.8e hits = %6d ",
	      nhits, acc_tab.theta_bar, acc_tab.hits ); 
      fprintf(prolixptr, "success = %6d acc_ratio = %5.2f\n", 
	      acc_tab.success, (double)acc_tab.success/(double)acc_tab.hits);
    }
  }

/* reset acceptance stats for next 'interval' */

  *m_success=0; 

/* close prolix file, if necessary */

  if ( myid == 0 )
    if ( prolix )
      fclose(prolixptr);
  
}
}



/*** FUNCTIONS FOR FILE I/O ***********************************************/

/*** WriteResults: writes the final states and annealing output to the *****
 *                 output file; this is an exception to the rule that I/O  *
 *                 functions should be in lj_sa.c, since it needs to know  *
 *                 about the tour arrays and ap parameters                 * 
 *                 which are static to move.c                              *
 ***************************************************************************/

void WriteResults(FILE *outptr, int precision)
{
  int i;

  /* ap.xx are static to move.c  */
  fprintf(outptr, "$annealing_parameters:\n");
  fprintf(outptr,"seed = %ld\n",ap.seed);
  fprintf (outptr,"initial temp = %f\n", ap.start_tempr);
  fprintf (outptr,"gain = %f\n",ap.gain_div_interval*(double)ap.interval);
  fprintf (outptr,"interval = %d\n",ap.interval);
  fprintf(outptr, "$$\n\n");

  /* state->tune.xxx are static globals declared in sa.h */
  fprintf(outptr, "$tune_parameters:\n"); 
  fprintf(outptr,"lambda = %f, lambda_mem_length_u =  %f\n",
		state->tune.lambda,state->tune.lambda_mem_length_u); 
  fprintf(outptr,"lambda_mem_length_v =  %f, control = %f\n",
	   state->tune.lambda_mem_length_v,state->tune.control); 
  fprintf(outptr,"initial_moves = %d, tau = %d\n",
	   state->tune.initial_moves, state->tune.tau);
  fprintf(outptr,"freeze_count = %d, update_S_skip = %d\n",
	   state->tune.freeze_count, state->tune.update_S_skip);
  fprintf (outptr,"criterion = %g\n", state->tune.criterion);

#ifdef MPI          
  fprintf (outptr,"mix_count = %d\n",state->tune.mix_interval);
#endif

  fprintf(outptr, "$$\n\n"); 


  /* DistP.xx are static from distributions.h  */
  fprintf(outptr, "$distribution_parameters:\n");
  fprintf (outptr,"distribution type=%d q=%lf\n",DistP.distribution,DistP.q);
  fprintf(outptr, "$$\n\n");

 
         /****************************************************************/
         /* equil_param.xxx are static to lsa.c ChuParam from sa.h       */ 
	 /* would need to write a new lsa.c routine to pass these        *
          *  back here, in case we need to output these                  */
         /* for now, I am not going to do this. LG 08:02                 */
         /*fprintf(outptr, "$equilibrate_parameters:\n");               */
	 /*fprintf (outptr, "end_T= %f, fix_T_skip= %d, fix_T_step= %d\n",
	  equil_param.end_T,equil_param.fix_T_skip,equil_param.fix_T_step); */
	 /*fprintf(outptr, "$$\n\n");                                    */ 
         /****************************************************************/
  
  fprintf(outptr, "$final_state:\n");
  fprintf(outptr,"\n annealing minimum cost is: %f", curr_cost); 
  fprintf(outptr," obtained in %d steps \n", nhits); 
  fprintf (outptr,"min tour is:\n");  
   for ( i=0; i<ncities; i++)
     {
          fprintf(outptr, "%d\t", curr_tour[i]); 
     }  /*end for print */
  fflush(outptr);
  fprintf(outptr, "\n$$\n\n");
  
 if (tour_debug) 
{ /* begin debug */
  fprintf(outptr, "$global_state:\n\n");
  fprintf(outptr,"\n For fun.absolute minimum cost is: %f\n",min_cost);  
  fprintf (outptr,"absolute min tour is:\n");  
  for ( i=0; i<ncities; i++)
    {
          fprintf(outptr, "%d\t", min_tour[i]); 
    }  /*end for print */
  fflush(outptr);
  fprintf(outptr, "$$\n\n");
} /* end debug */
 
  fprintf(outptr, "$annealing_output:\n");
  fprintf(outptr, "final_energy:\n");
  fprintf(outptr, "%.*f\n", precision, curr_cost);

 if (tour_debug) 
{ /* begin debug */
  fprintf(outptr, "global_energy:\n");
  fprintf(outptr, "%.*f\n", precision, min_cost);
} /* end debug */
fprintf(outptr, "max_count:\n"); /* this is number of iterations */
  fprintf(outptr, "%d\n", nhits);
  fprintf(outptr, "$$\n");
}

/*********************************************
* Print the current tour and its cost.       *
* commented out print to outfile.            *
**********************************************/
 void print_curr_tour(FILE* outfile)
{  /* begin print current tour */
int i;
fprintf(outfile,"\n current cost is: %f\n", curr_cost);  
fprintf (outfile,"current tour is:\n");  
 for ( i=0; i<ncities; i++)
  {
          fprintf(outfile, "%d\t", curr_tour[i]);  
  }  /*end for print */
 fflush(outfile);
}  /* end print current tour */



/***************************************************
* Print the annealing minimum tour & cost, and     *
* just for fun print the minimum tour and its cost *
****************************************************/
 void print_min_tour(FILE* outfile)
{  /* begin print min tour */
int i;
fprintf(outfile,"\n annealing minimum cost is: %f", curr_cost); 
fprintf(outfile," obtained in %d steps \n", nhits); 
fprintf (outfile,"min tour is:\n");  
 for ( i=0; i<ncities; i++)
  {
          fprintf(outfile, "%d\t", curr_tour[i]); 
  }  /*end for print */
fflush(outfile);
 if (tour_debug) 
{ /* begin debug */
fprintf(outfile,"\n For fun.absolute minimum cost is: %f\n",min_cost);  
fprintf (outfile,"absolute min tour is:\n");  
 for ( i=0; i<ncities; i++)
  {
          fprintf(outfile, "%d\t", min_tour[i]); 
  }  /*end for print */
fflush(outfile);
} /* end debug */
 }  /* end print min tour */



/*** FUNCTIONS THAT COMMUNICATE WITH OTHER SOURCE FILES ********************/

/*** SetProlix: sets flag for printing prolix output on acceptance stats ***
 *              and initializes the prolix file if required                *
 ***************************************************************************/

void SetProlix(int value, char *file, int init_flag)
{
  FILE       *prolixptr;

  const char *suffix = ".prolix";

/* sets the prolix file name static to move.c */

    prolixfile = (char *)calloc(MAX_RECORD, sizeof(char));
  prolixfile = strcpy(prolixfile, file);
  prolixfile = strcat(prolixfile, suffix);

/* this deletes an old .prolix file if present */

  if ( init_flag ) { 
    prolixptr = fopen(prolixfile, "w");
    if ( !prolixptr ) {
      perror("SetProlix");
      exit(1);
    }
    fclose(prolixptr);
  }

  prolix = value;
}


#ifdef MPI

/*** FUNCTIONS FOR SENDING MOVE STATES UPON MIXING (PARALLEL CODE ONLY) ****
 *   prototypes for these are in MPI.h; there's an extensive comment on    *
 *   how move state communication should be done at the beginning of lsa.c */

/*** MakeStateMsg: function to prepare a move state message which is then **
 *                 passed to other nodes via MPI; lsa.c doesn't know about *
 *                 the structs we use for acceptance statistics in         *
 *                 move(s).c, but we can safely assume that what we have   *
 *                 to send can be communicated as longs and doubles; thus, *
 *                 we split the move state message in two arrays, one for  *
 *                 the longs and one for the doubles; then we return the   *
 *                 arrays and their sizes to lsa.c                         *
 * to use for TSP need to pass:                                            * 
 *           curr_tour and curr_position arrays as well as curr_cost       *
 *           in addition to the nhits, nsweeps and acc_tab stuff           *
 ***************************************************************************/

void MakeStateMsg(unsigned char **buf, MPI_Aint padding, MPI_Aint *size)
{
  /*int    i;       don't need this anymore SRV 2025-11-15 */


/* calculate buffer size, compare with move parameters below */
/* Don't need to do that anymore - Seb RV 2025 November 14 */

/* allocate buffer */
  
    MPI_Alloc_mem(nbytes+padding, MPI_INFO_NULL, buf);

/* pack longs into longbuf */
/*Nope. */
  unsigned char *buf_pointer = *buf;
  /* equivalent to unsigned char *buf_pointer; buf_pointer=buf; */

  /*Since MPI_ALLOC_MEM gives us byte alignment for free,          *
   * we should use this in order to pack the massive arrays        *
   * first since this potentially will lead to more cache          *
   * hits, and more importantly extremely efficient vectorisation  *
   * of memcpy operations.                                         */
  memcpy(buf_pointer, curr_tour_and_pos, size_arr);
  buf_pointer += size_arr;
  memcpy(buf_pointer, &acc_tab.hits,sizeof(unsigned char));
  buf_pointer += sizeof(unsigned char);
  memcpy(buf_pointer, &acc_tab.success, sizeof(unsigned char));
  buf_pointer += sizeof(unsigned char);
  memcpy(buf_pointer, &ncities, sizeof(unsigned short));
  buf_pointer += sizeof(unsigned short);
  memcpy(buf_pointer, &nhits, sizeof(unsigned int));
  buf_pointer += sizeof(unsigned int);
  memcpy(buf_pointer,&nsweeps, sizeof(unsigned int));
  buf_pointer += sizeof(unsigned int);
  memcpy(buf_pointer,&curr_cost, sizeof(double));
  buf_pointer += sizeof(double);
  memcpy(buf_pointer,&acc_tab.theta_bar, sizeof(double));
  buf_pointer += sizeof(double);
  
  *size = nbytes; /*size now tells Lam where the move state ends
  so that it can append its information at the end */

}


/*** AcceptMsg: gets the move state message from lsa.c and reinstalls acc- *
 *              eptance statistics into move.c; see the comment for Make-  *
 *              StateMsg above for the rationale behind the two arrays     *
 *              that are passed                                            *
 *              DoMix in lsa.c does decision making for accepting this as  *
 *              the state for this node.                                   *
 ***************************************************************************/

void AcceptStateMsg(unsigned char **buf)
{
  int i;
  unsigned char *buf_pointer = *buf;
  /* equivalent to unsigned char *buf_pointer; buf_pointer=buf; */

  /*Since MPI_ALLOC_MEM gives us byte alignment for free,          *
   * we should use this in order to pack the massive arrays        *
   * first since this potentially will lead to more cache          *
   * hits, and more importantly extremely efficient vectorisation  *
   * of memcpy operations.                                         */
  memcpy( curr_tour_and_pos, buf_pointer,size_arr);
  buf_pointer += size_arr;
  memcpy(&acc_tab.hits,buf_pointer,sizeof(unsigned char));
  buf_pointer += sizeof(unsigned char);
  memcpy(&acc_tab.success, buf_pointer,sizeof(unsigned char));
  buf_pointer += sizeof(unsigned char);
  memcpy(&ncities, buf_pointer,sizeof(unsigned short));
  buf_pointer += sizeof(unsigned short);
  memcpy(&nhits, buf_pointer,sizeof(unsigned int));
  buf_pointer += sizeof(unsigned int);
  memcpy(&nsweeps, buf_pointer,sizeof(unsigned int));
  buf_pointer += sizeof(unsigned int);
  memcpy(&curr_cost, buf_pointer,sizeof(double));
  buf_pointer += sizeof(double);
  memcpy(&acc_tab.theta_bar, buf_pointer,sizeof(double));
  buf_pointer += sizeof(double);
/* unpack longs */
  }
#endif


/*** Utility cost functions for TSP tour stuff **************************/

/**********************************************
* Calc_new_cost calculates the cost of the    *
* tour using only city pairs that are proposed*
* to be swapped.  The proposed new cost is    *
* returned by this module                     *
***********************************************/
 double calc_new_cost(unsigned short *tour, unsigned short *swap, 
 double original_cost)

{/* begin update_cost*/
  double add_cost=0;/* new edge costs */
  double remove_cost=0;/* old edge costs */
  /* need to use modulo arithmetic to make nice with (i.e. generalize) the endpoints */
  /* in modulo arithmetic can't subtract a number.  You must add (dimension-number you want to subtract) */
  /*begin updating edges */
  if (((swap [0] == 0) && (swap[1] == ncities-1)) || ((swap [1] == 0) && (swap[0] == ncities-1)))
    {/* special treatment when both are endpoints */
       remove_cost += GetDistance(node_coords[tour[ncities-2]], node_coords[tour[ncities-1]]) ;
       remove_cost += GetDistance(node_coords[tour[0]], node_coords[tour[1]]) ;
       add_cost += GetDistance(node_coords[tour[1]], node_coords[tour [ncities-1]]) ;
       add_cost += GetDistance(node_coords[tour[0]], node_coords[tour [ncities-2]]) ;
     } /* both are endpoints */
  
  else if (((swap [1] - swap [0]) == -1)|| (swap [1] == 0 && swap[0] == 1) || (swap [1] == (ncities-2) && swap[0] == (ncities -1)))
    { /* Neighbors are special too.  Neighbors are swapping with swap[0] > swap[1]  */
      remove_cost += GetDistance(node_coords[tour[swap[0]]], node_coords[tour [(swap[0]+1)%ncities]]) ;
      remove_cost += GetDistance(node_coords[tour[swap[1]]], node_coords[tour [(swap[1]+ncities-1)%ncities]]) ;
      add_cost += GetDistance(node_coords[tour[swap[1]]], node_coords[tour [(swap[0]+1)%ncities]]) ;
      add_cost += GetDistance(node_coords[tour[swap[0]]], node_coords[tour [(swap[1]+ncities-1)%ncities]]) ;
    } /* neighbors are swapping with swap[0] > swap[1]  */
   else if (((swap [1] - swap [0]) == 1) || (swap [0] == 0 && swap[1] == 1) || (swap [0] == (ncities-2) && swap[1] == (ncities -1)))
    { /* Neighbors are special too.  Neighbors are swapping with swap[1] > swap[0]  */
      remove_cost += GetDistance(node_coords[tour[swap[0]]], node_coords[tour [(swap[0]+ncities-1)%ncities]]) ;
      remove_cost += GetDistance(node_coords[tour[swap[1]]], node_coords[tour [(swap[1]+1)%ncities]]) ;
      add_cost += GetDistance(node_coords[tour[swap[1]]], node_coords[tour [(swap[0]+ncities-1)%ncities]]) ;
      add_cost += GetDistance(node_coords[tour[swap[0]]], node_coords[tour [(swap[1]+1)%ncities]]) ;
    } /* neighbors are swapping with swap[1] > swap[0]  */
  else
    {  /* most likely not neighbors, not both endpoints.  This is the general case */
      remove_cost += GetDistance(node_coords[tour[swap[0]]], node_coords[tour [(swap[0]+ncities-1)%ncities]]) ;
      remove_cost += GetDistance(node_coords[tour[swap[0]]], node_coords[tour [(swap[0]+1)%ncities]]) ;
      remove_cost += GetDistance(node_coords[tour[swap[1]]], node_coords[tour [(swap[1]+ncities-1)%ncities]]) ;
      remove_cost += GetDistance(node_coords[tour[swap[1]]], node_coords[tour [(swap[1]+1)%ncities]]) ;
      add_cost += GetDistance(node_coords[tour[swap[1]]], node_coords[tour [(swap[0]+ncities-1)%ncities]]) ;
      add_cost += GetDistance(node_coords[tour[swap[1]]], node_coords[tour [(swap[0]+1)%ncities]]) ;
      add_cost += GetDistance(node_coords[tour[swap[0]]], node_coords[tour [(swap[1]+ncities-1)%ncities]]) ;
      add_cost += GetDistance(node_coords[tour[swap[0]]], node_coords[tour [(swap[1]+1)%ncities]]) ;
     } /* most likely not neighbors, not both endpoints */
 
 /*end for updating  edges */
  return original_cost + add_cost - remove_cost;
 } /* end update_cost*/




/*************************************************************
* compare new cost to global min cost; stores lowest in      * 
* min_tour and min_cost only for fun, so not in annealing run*
*************************************************************/
void global_min(void)
{  /* begin global_min */
int i; /*loop counter */
/* cost and tour must be saved */
if (curr_cost < min_cost)
   { /* we have a new minimum */
        min_cost = curr_cost ;
         for ( i=0; i<ncities; i++)
            {/*begin for */
                min_tour [i] = curr_tour [i];
            } /*end for */
 }/* we have a new minimum */
}  /* end global_min */

/* Needed for reading and writing state */
/* State stuff not working rn 2025-11-16 SRV */
//MoveState *MoveSave(void){
//  int i;
//  MoveState *move_stuff;
//  move_stuff = (MoveState *)malloc(sizeof(MoveState));
////
//  move_stuff->curr_tour     = (unsigned short*)calloc(ncities, sizeof(unsigned short));
//  move_stuff->curr_position = (unsigned short*)calloc(ncities, sizeof(unsigned short));
//  move_stuff->curr_cost     = curr_cost;
//  move_stuff->acc_tab_ptr   = acc_tab;
//  move_stuff->nhits         = nhits;
//  move_stuff->nsweeps       = nsweeps;
//  move_stuff->ncities       = ncities;
//
//  for (i=0; i<ncities; i++){
//     move_stuff->curr_tour[i]     = curr_tour[i];
//     move_stuff->curr_position[i] = curr_position[i];
//  }
//  return move_stuff;
//}
void FreeMoveState(MoveState *move_stuff){
  free(move_stuff->curr_tour);
  free(move_stuff->curr_position);
  free(move_stuff);
}
 
void FreeDistances(){
  free(distances);
}
