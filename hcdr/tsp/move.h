/*****************************************************************
 *                                                               *
 *   move.h                                                      *
 *                                                               *
 *****************************************************************
 *                                                               *
 *   written by John Reinitz                                     *
 *   modified by King-Wai Chu and Johannes Jaeger 
 *               and Lorraine Greenwald                          *
 *                                                               *
 *****************************************************************
 *                                                               *
 * move.h contains problem-specific stuff for the TSP            *
 * problem; the functions and structs here deal with evaluation  *
 * of the TSP cost function and the move generation for Lam SA   *
 *                                                               *
 *****************************************************************/
/******************************************************************
 * my name is now move.h and now I am consistant with 9.3.1 lsa.c *
 * LG 03-03: use unsigned shorts for tour arrays                  *
 *           and floats rather than doubles for edges             * 
 *           for  BIG tsp instances memory crunch                 *
 * LG 08-02: intergrate 9.3.1 lsa.c with tsp                      *
 * there are all types of application specific stuff in here.     *
 * LG 03-02: need q for gen visiting distribution input file      *
 ******************************************************************/
#ifndef MOVE_INCLUDED
#define MOVE_INCLUDED
#endif

/* this def needed for func. defs that refer to (* FILE) */
#ifndef _STDIO_INCLUDED
#include <stdio.h>    
#endif

/* following for structures & consts used thruout */
#ifndef GLOBAL_INCLUDED
#include "global.h"
#endif

/* following for StopStyle enum style and SAType funcs below */
#ifndef SA_INCLUDED
#include "sa.h"
#endif


/*** CONSTANTS *************************************************************/

#define THETA_MIN  1.            /* minimum value of theta_bar (move size) */
/* THETA_INIT not used as a constant for tsp because dependent *
 * on problem size. Set to problem size/2 in InitMoves         */ 
/* #define THETA_INIT 5.0       initial value for all theta_bar (move size) */

#define LOWBITS    0x330E                /* following two are for drand to */
#define BYTESIZE   8                                   /* erand conversion */


/*** A GLOBAL **************************************************************/

int random_tweak;        /* flag for tweaking at random (not sequentially) */
                         /*** not used in TSP  **********/




/*** STRUCTS ***************************************************************/

/* range struct for limits and such */

typedef struct Range {       
  double min;    
  double max;    
} Range;


/***************************************************************************
 * The following are annealing parameters that are not specific to the Lam *
 * algorithm. In general they should be used in moves.c or fly_sa.c but    *
 * *not* in lsa.c. In the data file, the members of the struct labeled RO  *
 * are read from the $annealing_input section. They are used for initial   *
 * conditions of the annealer and do not change during a run. Members      *
 * labeled OUT are written to the $annealing_output section upon comple-   *
 * tion of a run.                                                          */

/***************************************************************************
 * The variable of this type was named tsp_AP and now is either ap         *
 * or l_aparams for reading in from input file                             *
 **************************************************************************/ 

typedef struct {
  long   seed;                      /* seed for random number generator RO */
  double start_tempr;          /* the initial equilibration temperature RO */
  double gain_div_interval;            /* gain for proportional control of move size RO */
  /* Seb may as well precompute gain/interval */
  double stop_energy;                /* the final energy of the answer OUT */
  int    max_count;                      /* total number of iterations OUT */
  int    interval;       /* number of sweeps between updating theta_bar RO */
/***************************************************************************/
/* int distribution;        move generation distribution type       RO     */
                         /* 1 - uniform; 2 - exp; 3 - normal; 4 - lorentz  */
                         /* formerly dist_type in king's lj code           */
                         /* need q for gen visiting distribution input file*/
                         /* LG: 05-02 moved to distributions.h DistP       */
/***************************************************************************/

} AParms;

/* Struct for acceptance statistics, which are used for keeping the accep- *
 * tance ratio for each parameter as close as possible to .44; there's one *
 * such struct for each parameter to be tweaked                            */

typedef struct {
/*double acc_ratio; Seb RV 2025 November 14 - removed
 * acceptance ratio for parameter */
/*Changed hits and success to unsigned chars */
  double theta_bar;              /* theta bar is proportional to move size */
  unsigned char hits; /* number of moves since last call to UpdateControl() */
  unsigned char success; /* number of these moves that were accepted */
} AccStats;
                         /* determines when to call UpdateControl routine */

/* contains copies of the static variables of moves.c together with the 
 * current tour */
typedef struct {
  unsigned short *curr_tour;
  unsigned short *curr_position;
  double curr_cost;
  AccStats *acc_tab_ptr;
  unsigned int nhits;
  unsigned int nsweeps; 
  unsigned short ncities;
} MoveState;

/* Opts struct is used to save command line options in savestate.c */
typedef struct {
  char *inname; /* filename of input file */
  char *outname; /* filename of output file */
  int diff_outfile;
  char *argv;   /* original command line */
  StopStyle stop_flag;  /*stop criterion */
  int prolix_flag;
  int landscape_flag;
  int log_flag;
  int time_flag;
  long state_write;
  long print_freq;
  long captions;
  int precision;
  int quenchit; /* flag for quenchit mode */
  int equil;
#ifdef MPI
  int tuning;
  int covar_index;
  int write_tune_stat;
  int auto_stop_tune;
#endif
} Opts;

/*** FUNCTION PROTOTYPES ***************************************************/

/* tsp_sa.c: I/O functions for miscellaneous stuff modified from Yogi's stuff*/


/*** InitEquilibrate: reads the equilibrate section of the data file, ******
 *                    which is needed for equilibration runs               *
 ***************************************************************************/

void InitEquilibrate(FILE *fp);

/*** ReadTune: reads the tune_parameters section in a data file and ********
 *             turns a SAType structure to the caller                      *
 ***************************************************************************/

SAType ReadTune(FILE *fp);

/*** ReadAParameters: reads the AParm struct from an annealing_input sec- **
 *                    tion; these are the annealing parameters that are    *
 *                    not Lam-specific (and should NOT go into lsa.c)      *
 ***************************************************************************/

AParms ReadAParameters(FILE *fp);


/*** PrintTimes: prints user and wallclock times to the .times file ********
 ***************************************************************************/

void PrintTimes(FILE *fp, double *times);

/*** PrintEquil: writes an 'equilibrate_variance' section with 'title' *****
 *               to the stream specified by fp                             *
 ***************************************************************************/

void PrintEquil(FILE *fp, double *equil_var, char *title);

/* move.c: functions for move generation */

/* a function for finalizing a run */

/*** GetFinalInfo: collects stop energy and final count and returns them ***
 *                 to the caller                                           *
 ***************************************************************************/

AParms GetFinalInfo(void);


/* initializing functions */

/*** InitTSP: initializes the search space and the TSP cost function       *
 *            directly in move.c No longer returns starting energy 3-3-03  *
 ***************************************************************************/

void InitTSP(FILE *infile);


/***************************************************************************/
/*** InitMoves: initializes the following move.c-specific stuff: ***********
 *              - static annealing parameter struct (ap)                   *
 *              - initializes random number generator in lsa.c             *
 *              - initializes acc_tab for acceptance statistics            *
 *              it then returns the initial temperature to the caller      *
 ***************************************************************************/

double InitMoves(FILE *fp);

/* a move generation function used in move.c, but not lsa.c                */

/***************************************************************************/
/*** UpdateControl: each 'interval' number of steps, acceptance stats are **
 *                  updated here; note that each problem may have specific *
 *                  needs.  TSP is straight forward collection.            *
 *                  For LJ problem, we collect                             *
 *                  acceptance statistics for all particles in the same    *
 *                  struct, since all parameters are of the same order of  *
 *                  magnitude in this problem; more complicated cost func- *
 *                  tions like the fly gene circuit model needs move con-  *
 *                  trol for each parameter individually; this function    *
 *                  also prints prolix stuff, if required (-p)             *
 ***************************************************************************/
void UpdateControl(void);

/***************************************************************************/
/* miscellaneous functions */
/***************************************************************************/

/***************************************************************************/
/*** WriteResults: writes the final states and annealing output to the *****
 *                 output file                                             *
 ***************************************************************************/

void WriteResults(FILE *outptr, int precision);

/***************************************************************************/
/*** SetProlix: sets flag for printing prolix output on acceptance stats ***
 *              and initializes the prolix file if required                *
 *              enables and disables trace called from lsa.c as an option  *
 ***************************************************************************/

void SetProlix(int value, char *file, int init_flag);

 

/****************************************************/
/* tsp prototypes for functions in move.c           */
/****************************************************/
/****************************************
* init_cost initializes tour costs to 0 *
*****************************************/
void init_cost(void);

/**********************************************
* Calc_new_cost calculates the cost of the    *
* tour using only city pairs that are proposed*
* to be swapped.  The proposed new cost is    *
* returned by this module                     *
***********************************************/
double calc_new_cost(unsigned short *tour, unsigned short* swap, double cost);

/***************************************************
* allocates dynamic memory for the tsp tour arrays *
***************************************************/
int tour_allocate(void);

/*****************************************************
* deallocates dymanic memory for the tsp tour arrays * 
******************************************************/
void tour_deallocate(void);

/*************************************************************
* compare new cost to global min cost; stores lowest in      * 
* min_tour and min_cost only for fun, so not in annealing run*
*************************************************************/
void global_min(void);

/***********************************************************
* this module prints the not annealed minimum tour for fun *
??????????????????may want to get rid of this??????????
***********************************************************/
void print_min_tour(FILE* outfile);

/***************************************
* this module prints the current tour  *
****************************************/
void print_curr_tour(FILE* outfile);

/*******************************************
* generates starting tour and sets min tour *
********************************************/
double StartTour( void );

/*******************************************
 * State File functions                    *
 ******************************************/

void StateRead(char *statefile, Opts *options, MoveState *moveptr, double *state, 
    unsigned short *rand, double *delta);

void RestoreProlix();

MoveState *MoveSave(void);

void StateRm();

void RestoreMoves(MoveState *moveptr);

void FreeMoveState(MoveState *move_stuff);

Opts *GetOptions(void);

void RestoreOptions(Opts *options);

