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
#include "distributions.h"   /* problem independent distributions */
#include "initialize.h"

static int tour_debug=0;  /*debug flag for tours and edges */

/* Global Variables */
int n_points;
int n_edges;
int degree;
spin_weight *weights_and_spin_states;
double curr_cost = 0.0;  /*cost of the current tour*/ 

/*** STATIC VARIABLES ******************************************************/
static int *rand_cluster;
static int *shuffled_list;
static bool *used;

static AParms    ap;                /* static copy of annealing parameters */

/* SRV 2025-11-16 made acc_tab no longer a pointer and made it into
 * the actual struct itself */
static AccStats  acc_tab;   /* struct to accumulate acceptance statistics */

static unsigned  int    nhits;    /* number of moves since the start of execution */
static unsigned int    nsweeps;     /* number of sweeps since start of execution */


static int       prolix;           /* flag for printing stat info (prolix) */
                                   /* value of zero means no trace */
static char      *prolixfile;                   /* filename of prolix file */


/*** TSP TOUR VARIABLES ******************************************************/

/* Edited by seb: Changed shorts to ints because some of the test problems are massive */

static double delta_e = 0.0;  /* cost of the proposed new tour*/ 
static double min_cost = 0.0;  /* cost of the minimum tour found so far 
                                  for fun-not really SA methodology*/
#define SET_USED_EDGE(i,j,neighbour) do { \
  used[i*degree +j ]=true; \
  int i_rel; \
  for(i_rel = 0 ; VERTEX_EDGE(neighbour,i_rel).neighbour != i; ++i_rel); \
  used[neighbour*degree + i_rel]=true;\
}while(0)

#define USED_EDGE(i,j) used[i*degree+j]

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

int InitConfig(int _n_points, int _n_edges){
  n_points = _n_points;
  n_edges  = _n_edges ;
  degree   = 2*(n_edges)/n_points;
  weights_and_spin_states = 
    malloc(sizeof(*weights_and_spin_states)*(n_edges*2+n_points));
  shuffled_list = malloc(sizeof(*shuffled_list)*(n_points));
  rand_cluster = malloc(sizeof(*rand_cluster)*(n_points/2+1));
  for (int i = 0; i < n_points; ++i){shuffled_list[i]=i;}
  used = malloc(sizeof(*used)*n_points*degree);
  return 0;
}



double InitFromStartPoints(FILE *start_tour_file){
  ReadStartingPoints(start_tour_file);
  /*todo*/
  return 0;
}


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

  seedval  = ap.seed;

  xsubj[0] = LOWBITS;

  middle16 = (unsigned short)seedval;
  xsubj[1] = middle16;

  left     = seedval >> (BYTESIZE * sizeof(unsigned short));
  left16   = (unsigned short)left;
  xsubj[2] = left16;

  InitERand(xsubj);        /* makes the xsubj array static to random.c */

/* acc_tab is for statistics like acceptance ratio etc. */

  acc_tab.theta_bar = 1;
  /* THETA_INIT varies by tsp instance so cannot be a constant */
  /* so do not use: acc_tab->theta_bar = THETA_INIT; */
  acc_tab.hits      = 0;
  acc_tab.success   = 0;
 
  nhits   = 0;
  nsweeps = 0;

  /* Initialize some byte info that we will need later 
   * for state messages */
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
  nhits     = moveptr->nhits;
  nsweeps   = moveptr->nsweeps;
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

  acc_tab = *(moveptr->acc_tab_ptr);
  free(moveptr->acc_tab_ptr);                                                   

  free(moveptr);
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
* deallocates dynamic memory for the spin arrays   *
* called by final move in ising_sa.c               *
****************************************************/
void ConfigDeallocate(){
  free(weights_and_spin_states);
  free(used);
  free(rand_cluster);
  free(shuffled_list);
}


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
  double theta; /* move control variable to pick neighbors*/
/* make a move, get energy and return delta_e */



/* increase counters */
/* SRV Nov 19 2025 - switching nhits to be after nsweeps *
 * in order to force UpdateControl to only occur after   *
 * communication steps.                                  */
  acc_tab.hits++;

  nhits++;
/* update statistics if interval passed & at least one sweep completed */

/* Seb RV August 1st 2024: Changed this to happen after nsweeps is 
 * calculated. It is what happens in the CDR code so, idk I don't really
 * think it makes a difference. */

 /* Removed all of this logic entirely. UpdateControl can only occur 
 during an update_stats step. Seb RV Nov 19 2025 */

  /*pick first tour index randomly */
  /* uniform dist i=[0, prob_dimension-1] */

  /* control the neighbor pick the lam way */
  theta = generate_dev(acc_tab.theta_bar, DistP.distribution, DistP.q); 
  int theta_int;
  if (theta+1 > (double)n_points/2)
    theta_int = n_points/2;
  else
    theta_int = (int)theta+1;

  if (DistP.distribution == 7) /* tsp needs a positive value for theta */
   {theta = fabs(theta);}

  /*****************************************************************/
  /* round robin method is another way tried once upon a time      */
  /*      if (theta > neighbor_max)                                */
  /*         { theta  = theta % neighbor_max; }  * modulo division */
  /*****************************************************************/
  RandomCluster(theta_int);
  delta_e = CalcDiff(curr_cost);


 return delta_e;

}  /* end generate_move*/


/*** AcceptMove: sets new energy cost and tour in current cost and tour    * 
 *               for the next step and keeps track of the number of        *
 *               successful moves for acceptance statistics                *
 ***************************************************************************/

void AcceptMove(void)
{  /* begin accept_move*/
  for (int i = 0 ; rand_cluster[i]!=-1; ++i){
    VERTEX_SPIN(rand_cluster[i])*=-1;
  }
  curr_cost += delta_e;
}

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

void UpdateControl(int *m_success) 
{

  
/* Seb RV Nov 20 2025 */
/* Now UpdateControl has to manually check for itself.    * 
 *                                                         *
 * UpdateControl gets called every proc_tau moves, which   *
 * is a multiple of tau, so even if nsweeps=nhits this can *
 * can only be called during an UpdateStats step.           */
  if (!(nhits % ap.interval) ){

  FILE       *prolixptr;                            /* prolix file pointer */

/* open prolix file for appending new move stats */
  if ( prolix ) {
    prolixptr = fopen(prolixfile, "a");
      if ( !prolixptr ) {
	perror("UpdateControl");
	exit(1);
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
  acc_tab.theta_bar+= ap.gain_div_interval * (double)(*m_success-44);
  /* for all the trouble this stupid fucking function gave us
   * in trying to not make assumptions about sweeps, they
   * use this magic number that assumes we have 1 sweep every 100
   * moves. Whatever. SRV Nov 19 2025*/

    if ( acc_tab.theta_bar > (double)n_points/2 ) {
      acc_tab.theta_bar = (double)n_points/2;
    }
    else if ( acc_tab.theta_bar < 1 ) {
      acc_tab.theta_bar = 1;
    }
 
/* if -p: root node prints prolix information to prolix file */

  if ( prolix ) {
      fprintf(prolixptr, "nsteps = %8d bar = %10.8e hits = %6d ",
	      nhits, acc_tab.theta_bar, acc_tab.hits ); 
      fprintf(prolixptr, "success = %6d acc_ratio = %5.2f\n", 
	      acc_tab.success, (double)acc_tab.success/(double)acc_tab.hits);
    }

/* reset acceptance stats for next 'interval' */

  *m_success=0; 

/* close prolix file, if necessary */

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
  fprintf(outptr,"\n annealing minimum cost is: %f", CalcCost()); 
  fprintf(outptr," obtained in %d steps \n", nhits); 
  fprintf (outptr,"spin configuration is:\n");  
   for ( i=0; i<n_points; i++)
     {
          fprintf(outptr, "%d\t", VERTEX_SPIN(i)); 
     }  /*end for print */
  fflush(outptr);
  fprintf(outptr, "\n$$\n\n");
  
 if (tour_debug) 
{ /* begin debug */
} /* end debug */
 
  fprintf(outptr, "$annealing_output:\n");
  fprintf(outptr, "final_energy:\n");
  fprintf(outptr, "%.*f\n", precision, CalcCost());

 if (tour_debug) 
{ /* begin debug */
  fprintf(outptr, "global_energy:\n");
  fprintf(outptr, "%.*f\n", precision, min_cost);
} /* end debug */
fprintf(outptr, "max_count:\n"); /* this is number of iterations */
  fprintf(outptr, "%d\n", nhits);
  fprintf(outptr, "$$\n");
}


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


/*** Utility cost functions for TSP tour stuff **************************/

/**********************************************
* Calc_new_cost calculates the cost of the    *
* tour using only city pairs that are proposed*
* to be swapped.  The proposed new cost is    *
* returned by this module                     *
***********************************************/
 double CalcDiff(double original_cost)
{ 
  memset(used,0,n_points*degree);
  double diff = 0;
  int j = 0;
  while (rand_cluster[j] != -1){
    int curr_vert = rand_cluster[j]; 
    int curr_spin = VERTEX_SPIN(curr_vert);
    for (int i = 0; i < degree; ++i){
      int neighbour = VERTEX_EDGE(curr_vert, i).neighbour;
      double weight = VERTEX_EDGE(curr_vert, i).weight;
      int neighbour_spin = VERTEX_SPIN(neighbour);
      if (USED_EDGE(curr_vert, i)){
        diff -= (curr_spin == neighbour_spin) ? - 2 * weight : + 2 * weight;
      } else {
        diff += (curr_spin == neighbour_spin) ? - 2 * weight : + 2 * weight;
        SET_USED_EDGE(curr_vert,i,neighbour);
      }
    }
    j++;
  }
  return diff;
} /* end update_cost*/




/*************************************************************
* compare new cost to global min cost; stores lowest in      * 
* min_tour and min_cost only for fun, so not in annealing run*
*************************************************************/
/* Needed for reading and writing state */
/* State stuff not working rn 2025-11-16 SRV */
MoveState *MoveSave(void){
  MoveState *move_stuff;
  move_stuff = (MoveState *)malloc(sizeof(MoveState));

  move_stuff->acc_tab_ptr   = &acc_tab;
  move_stuff->nhits         = nhits;
  move_stuff->nsweeps       = nsweeps;

  return move_stuff;
}
void SaveInitialState(void **p){
  /*TODO*/
}

void RestoreInitialState(void *p){
  /*TODO*/
}

void Wolff(int base_vertex, double theta){
  memset(rand_cluster, -1, (n_points/2+1)*sizeof(*rand_cluster));
  memset(used, 0, n_points);
  int cluster_index = 1;
  used[base_vertex] = true;
  rand_cluster[0]  = base_vertex;
  int i = 0;
  
   while(rand_cluster[i] != -1 && cluster_index < n_points/2){
    int curr_vert = rand_cluster[i];
    for (int j = 0; j < degree; ++j){
      int neighbour = VERTEX_EDGE(curr_vert, j).neighbour;
      double weight = VERTEX_EDGE(curr_vert, j).weight;
      if (used[neighbour]) continue;
      double rand_unif = RandomReal();
      /* If weights are large, then spins are more likely *
       * to be aligned (stronger coupling)                */
      if (rand_unif < 1-exp(-theta*weight)){
        rand_cluster[cluster_index++] = neighbour;
        if (cluster_index >= n_points/2) break;
      }
      used[neighbour] = true;
    }
    i++;
  }
    /* VERTEX_OFFSET, VERTEX_SPIN, VERTEX_EDGE */
}

void RandomCluster(int theta){
  int temp;
  memset(rand_cluster, -1, (n_points/2+1)*sizeof(*rand_cluster));
  for(int i = 0 ; i < n_points ; ++i){
    int j = RandomInt(i+1);
    temp = shuffled_list[i];
    shuffled_list[i] = shuffled_list[j];
    shuffled_list[j] = temp;
  }
  for(int i = 0; i < theta ; ++i){
    rand_cluster[i] = shuffled_list[i];
  }
}

double  CalcCost(){
  int neighbour; 
  double weight;
  double cost = 0;
  for (int i = 0 ; i < n_points ; ++i){
    for (int j = 0 ; j < degree; ++j){
      neighbour = VERTEX_EDGE(i,j).neighbour;
      weight = VERTEX_EDGE(i,j).weight;
      if(neighbour > i) continue;
      cost += VERTEX_SPIN(i)*VERTEX_SPIN(neighbour) * weight;
    }
  }
  return cost;
}

double CalcCostMod(){
  for (int i = 0 ; rand_cluster[i]!=-1; ++i){
      VERTEX_SPIN(rand_cluster[i])*=-1;
    }
  double cost = CalcCost();
  for (int i = 0 ; rand_cluster[i]!=-1; ++i){
      VERTEX_SPIN(rand_cluster[i])*=-1;
    }
  return cost;
}

void FreeMoveState(MoveState *move_stuff){
  free(move_stuff);
  return;
}

void DEBUG_print_state(){
  for (int i = 0 ; i < n_points; ++i){
    printf("%d\t", VERTEX_SPIN(i));
    }
  printf("\n");
  }
