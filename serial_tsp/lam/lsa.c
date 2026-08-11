/*******************************************************************
 *
 *                                                                 *
 *   lsa.c                                                         *
 *                                                                 *
 *******************************************************************
 *                                                                 *
 *   Credits:                                                      *
 *                                                                 *
 *   Originally written by Jimmy Lam and Dan Greening              *
 *   Adaptation for continuous problems and original implemen-     *
 *   tation of the parallel algorithm by John Reinitz              *
 *   Tuning, equilibration and quenchit mode by King-Wai Chu       *
 *   Partly rewritten, extended & documented by Johannes Jaeger    *
 *                                                                 *
 *******************************************************************
 *                                                                 *
 *   see ../doc/prog_man.ps for details of how this works          *
 *                                                                 *
 *******************************************************************
 *                                                                 *
 * IMPORTANT: IF YOU EVER CHANGE ANYTHING IN THIS FILE, LET ALL    *
 *            YOUR FELLOW PROGRAMMERS KNOW WELL IN ADVANCE AND     *
 *            CONSULT WITH THEM IF THEY AGREE ON YOUR CHANGES!!    *
 *                                                                 *
 *******************************************************************
 *                                                                 *
 * Copyright (C) 1989-2003 John Reinitz                            *
 *                                                                 *
 * This program is free software; you can redistribute it and/or   *
 * modify it under the terms of the GNU General Public License as  *
 * published by the Free Software Foundation; either version 2 of  *
 * the License, or (at your option) any later version.             *
 *                                                                 *
 * This program is distributed in the hope that it will be useful, *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of  *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the   *
 * GNU General Public License for more details.                    *
 *                                                                 *
 * You should have received a copy of the GNU General Public       *
 * License along with this program; if not, write to the           *
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330,     *
 * Boston, MA  02111-1307, U.S.A.                                  *
 *                                                                 *
 *******************************************************************/


#ifdef ICC
#include <mathimf.h>
#else
#include <math.h>
#endif
#include <float.h>                                    /* for double limits */
#include <limits.h>                                  /* for integer limits */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>                        /* these two are for times() */
#include <sys/times.h>
#include <time.h>                                    /* this is for time() */
#include <unistd.h>          /* for command line option stuff and access() */

/* NOTE: do not ever dare to include moves.h in here or in any of the hea- */
/*       ders below; lsa.c must remain truly problem independent           */

#include <sa.h>
#include <error.h>
#include <random.h>


/* STATIC VARIABLES ********************************************************/

/* log and input/output file related variables *****************************/

static char   *inputfile;                        /* name of the input file */
static char   *outputfile;                      /* name of the output file */
static char   *statefile;                        /* name of the state file */
static char   *logfile;                    /* name of the global .log file */
static char   *landscapefile;                /* filename of landscape file */

/* some energy-related variables *******************************************/

static double energy;                                    /* current energy */
static double old_energy;
static double energies[1000];


static double S;                                 /* current inverse energy */
static double dS;                      /* delta S: change in S during move */
static double S_0;                      /* the initial inverse temperature */

static double exp_arg; /* the exponent of the Metropolis criterion (-dE/T) */

/* Lam stats stuff: estimators, stats and acceptance ratios ****************/

static double mean;          /* mean energy, collected from tau last steps */
static double vari;      /* energy variance, collected from tau last steps */

static double estimate_mean;              /* Lam estimator for mean energy */
static double estimate_sd;  /* Lam estimator for energy standard deviation */

static int    success;                       /* number of successful moves */

static double alpha;         /* the third term of the Lam schedule formula */
static double acc_ratio;    /* average acceptance ratio for all parameters */



/* Lam stats stuff: variables for calculating the estimators  **************/

/* mean estimator */

static double w_a;                       /* w_a is the weight for the mean */

static double usyy;         /* these parameters store intermediate results */
static double usxy;               /* for the updating formulas for A and B */
static double usxx;                       /* see Lam & Delosme, 1988b, p10 */
static double usx;
static double usy;
static double usum;

static double A;            /* A and B are the parameters for the rational */
static double B;                /* function for the estimation of the mean */


/* sd estimator */

static double w_b;         /* w_b is the weight for the standard deviation */

static double vsyy;         /* these parameters store intermediate results */
static double vsxy;               /* for the updating formulas for D and E */
static double vsxx;                       /* see Lam & Delosme, 1988b, p10 */
static double vsx;
static double vsy;
static double vsum;

static double D;            /* D and E are the parameters for the rational */
static double E;  /* function for the estimation of the standard deviation */


/* Lam stats stuff: variables related to tau *******************************/

static double Tau;     /* double version of tau to calculate mean and vari */
static int    proc_tau;  /* proc_tau = tau                     in serial   */
                         /* proc_tau = tau / (# of processors) in parallel */
static long   count_tau;  /* how many times we did tau (or proc_tau) moves */

/* the actual number of moves for collecting initial statistics ************/

static int    proc_init;                        /* number of initial moves */

/* stuff used by Frozen ****************************************************/

static double old_mean;                    /* old mean as stored by Frozen */
static int    counter;                           /* counter used by Frozen */

/* skip tells the annealer how often it should update S ********************
 * implemented by Lam for increasing performance, i.e. currently we only   *
 * update the temperature every 10 moves (possibly obsolete now)           */

static int    skip = -1;

/* flag used by Landscape generation ****************************************/
/*      Set by InitLandscape called from xxx_sa.c****************************/
static int    landscape = 0;  
                            
/* vars used for equilibration runs ****************************************/
/* note: these need to be static to be passed on to the file that writes   */
/*       them to wherever they need to be written to (see GetEquil())      */

static ChuParam equil_param;             /* equilibration parameter struct */

static double fix_T_avg   = 0.0;   /* overall energy average at fixed temp */
static double fix_T_var   = 0.0;   /* energy variance at fixed temperature */
static double r_est=0;
static double r_n_est=0;
  
/* variables used for timing */

/* these variables are used to evaluate real time using time() */
static double      start;                     /* wallclock time before run */
static double      finish;                     /* wallclock time after run */

/* these structs are used to evaluate user time using times() or MPI_Wtime */
static struct tms *cpu_start;                      /* user time before run */
static struct tms *cpu_finish;                      /* user time after run */


/* MAIN HERE ***************************************************************/
/* This should be pretty self-explanatory.                                 */

int main(int argc, char **argv )
{
  double *delta;                            /* used to store elapsed times */


/* code for timing: wallclock and user times */
                               
  cpu_start  = (struct tms *)malloc(sizeof(struct tms));      /* user time */
  cpu_finish = (struct tms *)malloc(sizeof(struct tms));

  times(cpu_start);

  start = time(NULL);     /* returns wallclock time since EPOCH (1/1/1970) */


/* initialize cost function and move state, do initial moves (or restore   */
/* annealing state if restart                                              */
   
  Initialize( argc, argv ); 
/* the following is for non-equlibration runs and equilibration runs that  */
/* have not yet settled to their equilibrium temperature                   */

  if ( (bench != 1) && ((equil != 1) || (1.0/S > equil_param.end_T)) )
    Loop();
/* there's an alternative Loop for equlibration runs at stable temperature */

  if ( equil == 1 )
    FixTLoop();

/* code for timing */

  if ( time_flag ){
    delta = GetTimes();                  /* calculates times to be printed */
      WriteTimes(delta);                            /* then write them out */
    free(delta);
  }
  return 0;

}






/*** INITIALIZING FUNCTIONS ************************************************/

/*** Initialize: calls ParseCommandLine first; then does either initial ****
 *               randomization and collecting Lam stats or restores state  *
 *               of the annealer as saved in the state file.               *
 ***************************************************************************/

void Initialize(int argc, char **argv)
{
  int    opt_index;         /* pointer to current argument of command line */
  int    stateflag = 0;                              /* state file or not? */
  double initial_temp;                 /* initial temperature for annealer */

/* allocate memory for static file names */

  inputfile = (char *)calloc(MAX_RECORD, sizeof(char));
  statefile = (char *)calloc(MAX_RECORD, sizeof(char));

/* allocate memory for static Lam parameters and dance partners */
  state = (NucStateType *)malloc(sizeof(NucStateType));

                        
/* parse the command line and return index to input file name; then we     *
 * initialize the static file names for inputfile and state file; these    *
 * two are needed here, since they are required for reading the state file *
 * whereas all other file names depend on the output file name and need to * 
 * be initialized after calling RestoreState() below                       */

  opt_index = ParseCommandLine(argc, argv);
  inputfile = strcpy(inputfile, argv[opt_index]);

/* state files: used for the case that a run terminates or crashes unex-   *
 * pectedly; we can then restore the state of the run *precisely* as it    * 
 * was before the crash by restarting it from the state file               *
 * editing it so that shmem procs get their own statefile also.            */

    sprintf(statefile, "%s.state", inputfile);
/* check if a state file exists (access() is in unistd.h) */

  if ( 0 == access(statefile, F_OK) ){
    stateflag = 1;
  }
 /* first get Lam parameters, initial temp and energy and initialize S_0 */
/* if we restore a run from a state file: call RestoreState() */

/* Seb, 12/12/2023, fixed a bug where else was missing causing the */
/* program to attempt to restore state even if stateflag = 0. */
  if ( !stateflag ) {
    initial_temp = InitialMove(argc, argv, opt_index, state, &energy );
    S_0 = 1./initial_temp;
  } else{
    initial_temp=InitialMove(argc, argv, opt_index, state, &energy);
    RestoreState(statefile, state, &energy);
  }
/* initialize those static file names that depend on the output file name */
  InitFilenames();
  proc_tau  = state->tune.tau;                       /* static copy to tau */
  proc_init = state->tune.initial_moves;             /* # of initial moves */
  Tau = (double)state->tune.tau;                  /* double version to tau */
                                             /* for calculating estimators */

/* if we're not restarting: do the initial moves for randomizing and ga-   *
 * thering initial statistics                                              */
  if ( !stateflag )
    InitialLoop();
                                      
/* write first .log entry and write first statefile right after init; note *
 * that equilibration runs are short and therefore don't need state files  *
 * which would be rather complicated because of all the stats collected    * 
 * during equilibration; for the exact same reasons we don't write state   *
 * files for tuning either; in case of a restart, we just restore the .log */

  if ( !equil && !bench && !nofile_flag) {
    if ( !stateflag ) {
      WriteLog();
	StateWrite(statefile);
    } else 
      RestoreLog();
  }
  
}




/*** InitFilenames: initializes static file names that depend on the out- **
 *                  put file name (i.e. this *must* be called after we     *
 *                  have called RestoreState())                            *
 ***************************************************************************/

void InitFilenames(void)
{

/* allocate memory for static file names */

  logfile    = (char *)calloc(MAX_RECORD, sizeof(char));

/* mixlogfile added by SRV on Jul 31 2023. See DoMix and WriteMixLog for   *
 * more details                                                            */

/* if we're not using -w: output = inputfile */

  if ( !outputfile )
    outputfile = inputfile;

/* the global .log file: stores iterations, temperature, change in tempe-  *
 * rature, global means and standard deviation, Lam estimators for mean    *
 * and sd as well as global acceptance ratios over an annealing run        */

    sprintf(logfile, "%s.log", outputfile);

                                               
}



/*** InitialLoop: performs the two sets of initial moves: ******************
 *                   1. randomizing moves (not parallelized)               *
 *                   2. loop for initial collection of statistics          *
 ***************************************************************************/

void InitialLoop(void)
{
  int         i;                                     /* local loop counter */

  double      energy_change;               /* change of energy during move */




/* randomize initial state; throw out results; DO NOT PARALLELIZE! */
  for (i=0; i<state->tune.initial_moves; i++) {
	/* make a move: will either return the energy change or FORBIDDEN_MOVE */

    energy_change = GenerateMove();
   
/* Metropolis stuff here; we usually want FORBIDDEN_MOVE to be very large  *
 * that's why we want to prevent overflows here (hence the 'if')           */

    if ( energy_change != FORBIDDEN_MOVE )
      exp_arg = -S_0 * energy_change;  

    /* MIN_DELTA provides a min. probability with which any move is accepted */

    if ( exp_arg <= MIN_DELTA )
      exp_arg = MIN_DELTA;  
      
    /* below, we apply the Metropolis criterion to accept or reject a move */

    if ( energy_change == FORBIDDEN_MOVE )
      RejectMove();
    else if ( (energy_change <= 0.0) || (exp(exp_arg) > RandomReal()) ) {
      energy += energy_change;
      AcceptMove();
    } else
      RejectMove();
  }   /* end randomize initial state */
  
/* set all stats to zero, collection starts below */

  mean      = 0.0;
  vari      = 0.0;
  success   = 0;

/* loop to collect initial statistics; this one is parallelized */

  for ( i=0; i<proc_init; i++ ) {
  /* make a move: will either return the energy change or FORBIDDEN_MOVE */

    energy_change = GenerateMove();

/* Metropolis stuff here; we usually want FORBIDDEN_MOVE to be very large  *
 * that's why we want to prevent overflows here (hence the 'if')           */
    if ( energy_change != FORBIDDEN_MOVE )
      exp_arg = -S_0 * energy_change; 

    /* MIN_DELTA provides a min. probability with which any move is accepted */

    if( exp_arg <= MIN_DELTA )
      exp_arg = MIN_DELTA;

    /* below, we apply the Metropolis criterion to accept or reject a move */

    if ( energy_change == FORBIDDEN_MOVE )  
      RejectMove();
    else if ( (energy_change <= 0.0) || (exp(exp_arg) > RandomReal()) ) {
      energy += energy_change;
      AcceptMove();
      success++;
    } else
      RejectMove();

/* collect stats */

    mean += energy;        
    vari += energy * energy;

  }
 
    
/* global stats are calculated here */

  mean     /= (double)state->tune.initial_moves;
  vari      = vari / ((double)state->tune.initial_moves) - mean * mean;
  acc_ratio = ((double)success) / ((double)state->tune.initial_moves); 

/* initialize Lam parameters used for calculating Lam estimators */
  
  InitializeParameter();
  if (r_flag)
    memset(energies, 0, sizeof(energies));
}



/*** InitializeParameter: initializes variables for Lam annealing: this is *
 *                        executed after doing the initial steps;          *
 *                        the local parameters only need to be set for tu- *
 *                        ning code                                        *
 ***************************************************************************/

void InitializeParameter(void)
{
  double   d;             /* d is used to store intermediate results below */



/* 1. set global parameters (serial and parallel code) *********************/

/* set estimators to stats collected during initializing phase of run */

  estimate_sd   = sqrt(vari);
  estimate_mean = mean;

/* initialize A,B,D,E according to Lam & Delosme, 1988b, p10 */

  A = estimate_sd * estimate_sd / (estimate_mean * estimate_mean);
  B = (1.0 / estimate_mean) - (A * S_0);
  D = estimate_sd / estimate_mean;
  E = (1.0 / estimate_sd) - (D * S_0); 

/* initialize these intermediate variables for updating funcs for A,B,D,E */

  usum = vsum = 1.0;
  usxy = usxx = usx = 0.0;
  usy  = 1.0 / estimate_mean;
  usyy = usy * usy;
  vsxy = vsxx = vsx = 0.0;
  vsy  = 1.0 / estimate_sd;
  vsyy = vsy * vsy;

/* set the initial temperature and the initial delta S */

  S  = S_0; 
  dS = 0.5 / estimate_sd;            /* keep--may not need--based on s_0=0 */

/* alpha is the third term of the main Lam schedule formula */

  d     = (1.0 - acc_ratio) / (2.0 - acc_ratio);
  alpha = 4.0 * acc_ratio * d * d;


/* weights determine how the estimators are sampled for times before tau */

  InitializeWeights();

}



/*** InitializeWeights: initialize weights a and b for calculating Lam *****
 *                      estimators; these weights are computed from the    *
 *                      lambda memory length products                      *
 ***************************************************************************/

void InitializeWeights(void)
{
  FILE  *logptr;

/* w_a is the weight for the mean */

  w_a = state->tune.lambda_mem_length_u / state->tune.lambda;
  w_a = 1.0 - state->tune.tau / w_a;
  if (w_a < 0.0) 
    w_a = 0.0;          

/* w_b is the weight for the standard deviation */

  w_b = state->tune.lambda_mem_length_v / state->tune.lambda;
  w_b = 1.0 - state->tune.tau / w_b;
  if (w_b < 0.0) 
    w_b = 0.0;            


  if ( !equil && !bench && !nofile_flag ) {
      logptr = fopen(logfile, "w");
      if ( !logptr )
	file_error("InitializeWeights");
      fprintf(logptr, "InitializeWeights:  w_a = %g w_b = %g\n", w_a, w_b );
      fclose(logptr);
      if ( log_flag )
	printf("InitializeWeights:  w_a = %g w_b = %g\n", w_a, w_b );
  }
}




/*** MAIN LOOP AND UPDATE FUNCTIONS ****************************************/

/*** Loop: loops (making moves, updating stats etc.) until the system is ***
 *         considered frozen according to the stop criterion               *
 ***************************************************************************/

void Loop(void)
{

  int    i;                                          /* local loop counter */
  double energy_change;                                   /* local Delta E */
  double d;                /* difference between energy and estimated mean */
  
/* quenchit mode: set temperature to (approximately) zero immediately */

  if ( quenchit )
    S = DBL_MAX;

/* loop till the end of the universe (or till the stop criterion applies) */

  while (1) {
     
/* reset statistics */

    mean    = 0.0;  
    vari    = 0.0;
    success = 0;
    
    
/* do proc_tau moves here */
    
    for (i=0; i<proc_tau; i++) {    
        if (r_flag){
          old_energy=energy;
          }
/* make a move: will either return the energy change or FORBIDDEN_MOVE */

          energy_change = GenerateMove();

      if (energy+energy_change<0)
        printf("%f %f\n",energy_change,energy);    
/* Metropolis stuff here; we usually want FORBIDDEN_MOVE to be very large  *
 * that's why we want to prevent overflows here (hence the 'if'); we also  *
 * need to avoid overflows with quenchit (where S is (almost) infinite!)   */
      
      if ( !quenchit && (energy_change != FORBIDDEN_MOVE ) )
	exp_arg = -S * energy_change;

/* MIN_DELTA provides a min. probability with which any move is accepted */

      if ( (exp_arg <= MIN_DELTA) )
	exp_arg = MIN_DELTA;         
      
/* below, we apply the Metropolis criterion to accept or reject a move; in *
 * quenchit mode, only lower energies are accepted                         */
            if ( energy_change == FORBIDDEN_MOVE ) {
	RejectMove();
      } else if ( (energy_change <= 0.0) || 
		  ( (!quenchit) && (exp(exp_arg) > RandomReal()) ) ) {	
	energy += energy_change;

	AcceptMove();
	success++;

	
      } else {
	RejectMove();
      }

/* update statistics */

      mean    += energy;
      d        = energy - estimate_mean;
      vari    += d * d;


/* update temperature every 'skip' steps; this was put here by Jimmy Lam,  *
 * probably to save computation time on his old Spark; I guess it's obso-  *
 * lete by now, but since it doesn't seem to do any harm and saves us some *
 * time, we left it in here                */
        if ( !quenchit ) 
	if ( --skip <= 0 )
	  UpdateS();      
  if (r_flag){
    r_est+=(energy-estimate_mean)/(old_energy-estimate_mean);
    r_n_est+=(energy-estimate_mean)/(energies[100*((count_tau-9)%10)+i]-estimate_mean);
    energies[100*(count_tau%10)+i]=energy;
  }
    
    }                 /* this is the end of the proc_tau loop */            
/* have done tau moves here: update the 'tau' counter */

    count_tau++;  
/* calculate mean, variance and acc_ratio for the last tau steps; i is     *
 * passed as an argument for checking if all local moves add up to Tau     */
    UpdateStats(i);

/* check if the stop criterion applies: annealing and tuning runs (that    *
 * aren't stopped by the tuning stop criterion) leave the loop here; equi- *
 * libration runs exit below                                               */


  
/* update Lam stats: estimators for mean, sd and alpha from acc_ratio (we  *
 * don't need this in quenchit mode since the temperature is fixed to 0)   */
    if ( !quenchit ) 
      UpdateParameter(); 


/* write the log every print_freq * tau (not proc_tau!) */

    if ( (count_tau % print_freq == 0) && !equil && !nofile_flag ){
      if(r_flag){
        r_est/=print_freq*proc_tau;
        r_n_est/=print_freq*proc_tau;
      }
      WriteLog();                     
      if (r_flag){
        r_est=0;
        r_n_est=0;
      }
    }
/* the state file gets written here every state_write * tau */
    if ( (count_tau % state_write == 0) && !equil && !nofile_flag )
      StateWrite(inputfile); 
  if (Frozen()){
    FinalMove();
    return;
  }

  }                                /* this is the end of the while(1) loop */

}                          /* this is the end of Loop */



/*** UpdateS: update inverse temperature S at every Sskip step *************
 ***************************************************************************/

void UpdateS(void)
{
  register double    d;              /* used to store intermediate results */

  S += dS;                   /* here, inverse temperature is updated by dS */

/* we need to update Lam parameters here, since S has changed; A, B, C and */
/* D get updated in UpdateParameters                                       */

  estimate_mean = 1.0 / (A*S + B);    /* for temperature updating formulas */
  estimate_sd   = 1.0 / (D*S + E);        /* see Lam & Delosme, 1988b, p10 */

  

  d   = S * estimate_sd;             /* intermediate for the specific heat */

/* following lines implement the main Lam schedule formula */

  dS  = state->tune.lambda * alpha / (d*d * estimate_sd);
  dS *= state->tune.update_S_skip;       /* ... we have to muliply by skip */



/* reset skip */

  skip = state->tune.update_S_skip;

}



/*** UpdateStats: updates mean, variance and acc_ratio after tau moves *****
 ***************************************************************************/

void UpdateStats(int i)
{
  mean /= Tau;                                  /* collect some statistics */
  vari /= Tau;


  acc_ratio = ((double)success) / Tau ;         /* update acceptance ratio */
  }


/*** UpdateParameter: update parameters A, B, D and E and the estimators ***
 *                    for mean and standard deviation for the current S    *
 ***************************************************************************/

void UpdateParameter(void)
{
  register double   d;    /* d is used to store intermediate results below */


/* this part of the code updates the estimator for the mean */

  d     = 1.0 / mean; 

/* first: multiply all intermediate vars by weights */

  usyy *= w_a;
  usxy *= w_a;
  usy  *= w_a;
  usx  *= w_a;
  usxx *= w_a;
  usum *= w_a;

/* then: update all intermediate vars */

  usyy += d*d;                      
  usxy += S*d;
  usy  += d;
  usx  += S;
  usxx += S*S;
  usum += 1.0;
                        
/* ... and use intermediate vars to update A and B ... */
/* Seb: Added a statement to check for division by 0. */
  if (usum*usxx-usx*usx)
    A = (usum * usxy - usx * usy) / (usum * usxx - usx * usx);
  else
    A = (usum * usxy - usx * usy) / DBL_MIN;

  B = (usy - A * usx)/usum;

/* ... which are then used to update the estimator for the mean */

  estimate_mean = 1.0 / (A * S + B);


/* this part of the code updates the estimator for the standard deviation */

  if (vari > 0.0) {

    d     = 1.0 / sqrt(vari);

/* first: multiply all intermediate vars by weights */

    vsyy *= w_b;       
    vsxy *= w_b;
    vsy  *= w_b;
    vsx  *= w_b;
    vsxx *= w_b;
    vsum *= w_b;

/* then: update all intermediate vars */

    vsyy += d*d;                     
    vsxy += S*d;
    vsy  += d;
    vsx  += S;
    vsxx += S*S;
    vsum += 1.0;
                        
/* ... and use intermediate vars to update D and E ... */

    D = (vsum * vsxy - vsx * vsy) / (vsum * vsxx - vsx * vsx);
    E = (vsy - D * vsx) / vsum;
  }

/* ... which are then used to update the estimator for the std dev */

  estimate_sd = 1.0 / (D*S + E);

/* alpha corresponds to the third term in the main Lam schedule formula,   *
 * which is a measure of how efficiently the space state is samples; this  *
 * term is at a maximum for acc_ratio = 0.44, see Lam & Delosme, 1988b, p1 */

  d = (1.0 - acc_ratio) / (2.0 - acc_ratio);   
  alpha = 4.0 * acc_ratio * d * d;

}                                 

/*** Frozen: returns TRUE if frozen, FALSE otherwise; 'frozen' is defined **
 *           by 'freeze_count', 'stop_flag' and 'criterion' (see also sa.h *
 *           for a more extensive comment on this)                         *
 ***************************************************************************/

int Frozen(void)
{
  double  delta;
  if(stop_flag == proportional_freeze)
    delta = (mean - old_mean)/mean;

  else if (stop_flag == absolute_freeze)
    delta = mean - old_mean;
  else if (stop_flag == absolute_energy)
    delta = mean;

  if (delta <= 0.)
    delta = -delta;

  if (delta <= state->tune.criterion) {
    counter++;
  }
  else {
    counter  = 0;
    old_mean = mean;
  }
  return(counter >= state->tune.freeze_count);
}






/*** FUNCTIONS REQUIRED FOR EQUILIBRATION RUNS ONLY ************************/

/*** FixTLoop: loops (making moves, updating stats etc.) keeping the temp- *
 *             erature fixed for an equilibration run                      *
 ***************************************************************************/

void FixTLoop(void)
{
  int    i;                                               /* loop counters */
  int    j;

  char   *varfile;                            /* global variance file name */
  FILE   *varptr;                          /* global variance file pointer */

  char   *acfile;                             /* autocorrelation file name */
  FILE   *acptr;                           /* autocorrelation file pointer */

  double energy_change;                                   /* local Delta E */

/* array(s) to store energies */

  double *energy_storage;     /* array used to store intermediate energies */

/* below are vars that are used to calculate statistics at equilibrium */

/* first some variables to calculate average energy and variances */

  double tmp1        = 0.0;                   /* three temporary variables */
  double tmp2        = 0.0;                   /* for calculating variances */
  double var_sum     = 0.0;      

  double instant_avg;     /* used for calculating the inst. average energy */

/* below some variables used to calculate autocorrelation */


  double *gamma;                                  /* array for covariances */
  double *rho;                                /* array of autocorrelations */

                    /* array of intervals for calculating autocorrelations */
  int    h[56] = {  0,
     1,      2,      3,      4,      5,      6,      7,      8,      9,
     10,     20,     30,     40,     50,     60,     70,     80,     90,
     100,    200,    300,    400,    500,    600,    700,    800,    900,
     1000,   2000,   3000,   4000,   5000,   6000,   7000,   8000,   9000,
     10000,  20000,  30000,  40000,  50000,  60000,  70000,  80000,  90000,
     100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000,
     1000000 }; 

  int    num_of_h = 56;                   /* number of elements in h array */
  int    ntemp;                     /* temporary var for max. step minus h */



/* allocate serial-specific file names */

  varfile        = (char *)calloc(MAX_RECORD, sizeof(char));
  acfile         = (char *)calloc(MAX_RECORD, sizeof(char));

/* allocate the energy array */

  energy_storage = (double *)calloc(equil_param.fix_T_step+1, 
				    sizeof(double)); 

/* allocate autocorrelation arrays */

  gamma          = (double *)calloc(num_of_h, sizeof(double));
  rho            = (double *)calloc(num_of_h, sizeof(double));

  sprintf(varfile, "%s.var", outputfile);
  sprintf(acfile,  "%s.ac" , outputfile);


/* Equilibrate and throw away statistics */

  for (i=0; i < equil_param.fix_T_skip; i++) {

/* randomize initial state; throw out results; DO NOT PARALLELIZE! */

    energy_change = GenerateMove();
    
/* Metropolis stuff here; we usually want FORBIDDEN_MOVE to be very large  *
 * that's why we want to prevent overflows here (hence the 'if')           */

    if ( energy_change != FORBIDDEN_MOVE )
      exp_arg = -S * energy_change; 

/* MIN_DELTA provides a min. probability with which any move is accepted */

    if( exp_arg <= MIN_DELTA )
      exp_arg = MIN_DELTA;    
    
/* below, we apply the Metropolis criterion to accept or reject a move */

    if ( energy_change == FORBIDDEN_MOVE ) { 
      RejectMove();
    } else if ((energy_change <= 0.0) || (exp(exp_arg) > RandomReal()) ) {
      energy += energy_change;
      AcceptMove();
    } else {
      RejectMove();
    }

  }

/* Equilibrate and gather statistics */

  for (i=0; i <= equil_param.fix_T_step; i++) {

/* make a move: will either return the energy change or FORBIDDEN_MOVE */

    energy_change = GenerateMove();

/* Metropolis stuff here; we usually want FORBIDDEN_MOVE to be very large  *
 * that's why we want to prevent overflows here (hence the 'if')           */

     if ( energy_change != FORBIDDEN_MOVE )
      exp_arg = -S * energy_change;      

/* MIN_DELTA provides a min. probability with which any move is accepted */

    if(exp_arg <= MIN_DELTA ) 
      exp_arg = MIN_DELTA;    

/* below, we apply the Metropolis criterion to accept or reject a move */

    if (energy_change == FORBIDDEN_MOVE)
      RejectMove();
    else if ((energy_change <= 0.0) || (exp(exp_arg) > RandomReal())) {
      energy += energy_change;
      AcceptMove();
      if (landscape) {
       WriteLandscape(landscapefile,i,energy_change);}
    }
    else
      RejectMove();


/* serial code: collect energies at every step */

    energy_storage[i] = energy;


  }                                    /* end of loop to gather statistics */

/* open files to write into */

  varptr = fopen(varfile, "w");                   /* this is the .var file */
  if ( !varptr )
    file_error("FixTLoop");
  
  acptr = fopen(acfile, "w");                      /* this is the .ac file */
  if ( !acptr )
    file_error("FixTLoop");

/* first calculate local average and variances */

/* calculate overall energy average */

  for (i=0; i <= equil_param.fix_T_step; i++)
    fix_T_avg += energy_storage[i];
  fix_T_avg /= equil_param.fix_T_step;

/* initialize variables */

  tmp1 = (energy_storage[0]-fix_T_avg)*(energy_storage[0]-fix_T_avg);
  tmp2 = energy_storage[0]-fix_T_avg;
  instant_avg = energy_storage[0];

/* print captions */
  fprintf(varptr, "# nsteps     variance     ");
  fprintf(varptr, "inst. avg.   overall avg.\n\n");
  
/* calculate and print local variances here */

  for(i=1; i <= equil_param.fix_T_step; i++) {

    instant_avg += energy_storage[i];

    tmp1    += (energy_storage[i]-fix_T_avg)*(energy_storage[i]-fix_T_avg);
    tmp2    += (energy_storage[i]-fix_T_avg);
    var_sum += (tmp1 - tmp2*tmp2/(i+1)) / i;

    fix_T_var = var_sum / (i+1);

/* print variances and averages every 1000 steps or when finished */

    if ( !(i % 1000) || (i == equil_param.fix_T_step) ) {
      fprintf(varptr,"%8d %12.5E   %12.5E   %12.5E\n",
              i, fix_T_var, instant_avg/(i+1), fix_T_avg);
      fflush(varptr);
    }

  }

/* calculate auto-correlation below: this is important for evaluating the  */
/* appropriate number of initial steps for a problem (see Kingwai's thesis */
/* section 2.5.3 (pp. 23 - 24)                                             */
/* gamma is the covariance for a certain interval (taken from the h array) */
/* rho is the auto-correlation for the same interval                       */

  fprintf(acptr, "#      h           rho          gamma\n\n");

  for(i=0; i<num_of_h; i++) {     /* calculate autocorrelation for all h's */

    ntemp = equil_param.fix_T_step - h[i];  

/* note: if h[i] is bigger than fix_T_step, gamma and rho will simply be 0 */

    gamma[i] = 0.0;
    for (j=0; j<=ntemp; j++)
      gamma[i] += 
	(energy_storage[j+h[i]]-fix_T_avg)*(energy_storage[j]-fix_T_avg);
    gamma[i] /= equil_param.fix_T_step;
 
    rho[i] = gamma[i] / gamma[0];
    fprintf(acptr, "%8d   % 10.8f   %12.5E\n", h[i], rho[i], gamma[i]);

  }

/* clean up and go home ... */

  fclose(varptr);
  free(varfile);
  fclose(acptr);
  free(energy_storage);
  free(acfile);
  free(gamma);
  free(rho);

  FinalMove();
  return;

} 



/*** SetEquilibrate: simply makes the equil_param struct static to lsa.c ***
 ***************************************************************************/

void SetEquilibrate(ChuParam ep)
{
  equil_param = ep;
}

/*** GetEquil: returns the results of an equilibration run *****************
 ***************************************************************************/

void GetEquil(double *equil_var)
{

}






/*** FUNCTIONS USED TO COMMUNICATE WITH OTHER FILES  ***********************/

/****************************************************************************  
 *** InitLandscape: sets flag for printing landscape output and acceptance****
 *                 landscape and initializes the landscape file names      ***
 *                 called from xxx_sa.c to make filenames static and       ***
 *                 set landscape flag                                      ***
 ****************************************************************************/

void InitLandscape(int value, char *file)
{

  const char *suffix = ".landscape"; /* landscape file in equilibrate */

/* sets the landscape and acceptance landscape file names static to lsa.c */

  landscapefile = (char *)calloc(MAX_RECORD, sizeof(char));
  landscapefile = strcpy(landscapefile, file);
  landscapefile = strcat(landscapefile, suffix);

  landscape = value;


}
/*** GetLamstats: returns Lam statistics in an array of doubles; used to ***
 *                store Lam statistics in a state file                     *
 ***************************************************************************/

double *GetLamstats(void)
{
  double *stats;

  stats = (double *)calloc(31, sizeof(double));

  stats[0] = (double)counter;

  stats[1]  = old_mean;
  stats[2]  = energy;

  stats[3]  = mean;
  stats[4]  = vari;

  stats[5]  = estimate_mean;
  stats[6]  = estimate_sd;

  stats[7]  = S;
  stats[8]  = dS;
  stats[9]  = S_0;

  stats[10] = alpha;
  stats[11] = acc_ratio;

  stats[12] = w_b;
  stats[13] = vsyy;
  stats[14] = vsxy;
  stats[15] = vsxx;
  stats[16] = vsx;
  stats[17] = vsy;
  stats[18] = vsum;
  stats[19] = D;
  stats[20] = E;

  stats[21] = w_a;
  stats[22] = usyy;
  stats[23] = usxy;
  stats[24] = usxx;
  stats[25] = usx;
  stats[26] = usy;
  stats[27] = usum;
  stats[28] = A;
  stats[29] = B;

  stats[30] = (double)count_tau;

  return(stats);
}



/*** GetTimes: returns a two-element array with the current wallclock and **
 *             user time to be saved in the state file; for parallel code  *
 *             we average the times for all processes                      *
 ***************************************************************************/

double *GetTimes(void)
{
  double *delta;
  double clk_tck = (double)CLOCKS_PER_SEC; /* system-specific clock tick duration */


  delta = (double *)calloc(2, sizeof(double));
                                                      /* measure user time */
  times(cpu_finish);
                                                    /* then wallclock time */
  finish = time(NULL);

  delta[0] = finish - start;
  delta[1] = (cpu_finish->tms_utime - cpu_start->tms_utime)/clk_tck;

  return delta;
}



/*** SetOutname: sets the output filename in lsa.c; this is necessary to ***
 *               have the diverse .log and .ac and .mb etc files have the  *
 *               name of the output file if -w is chosen                   *
 ***************************************************************************/

void SetOutname(char *outname)
{
  outputfile = outname;
}





/*** FUNCTIONS TO RESTORE THINGS IN LSA.C AFTER A RESTART ******************/

/*** RestoreLamstats: restores static Lam statistics in lsa.c from an ******
 *                    array of doubles; used to restore runs from a state  *
 *                    file.                                                *
 ***************************************************************************/

void RestoreLamstats(double *stats)
{
  counter = (int)rint(stats[0]);

  old_mean      = stats[1];
  energy        = stats[2];
  mean          = stats[3];
  vari          = stats[4];

  estimate_mean = stats[5];
  estimate_sd   = stats[6];

  S             = stats[7];
  dS            = stats[8];
  S_0           = stats[9];
  alpha         = stats[10];
  acc_ratio     = stats[11];

  w_b           = stats[12];
  vsyy          = stats[13];
  vsxy          = stats[14];
  vsxx          = stats[15];
  vsx           = stats[16];
  vsy           = stats[17];
  vsum          = stats[18];
  D             = stats[19];
  E             = stats[20];

  w_a           = stats[21];
  usyy          = stats[22];
  usxy          = stats[23];
  usxx          = stats[24];
  usx           = stats[25];
  usy           = stats[26];
  usum          = stats[27];
  A             = stats[28];
  B             = stats[29];

  count_tau = lrint(stats[30]);

  free(stats);
}



/*** RestoreTimes: restores the wallclock and user times if -t is used *****
 ***************************************************************************/

void RestoreTimes(double *delta)
{
  double clk_tck = (double)CLOCKS_PER_SEC; /* system-specific clock tick duration */

  start = time(NULL);
  start -= delta[0];

  times(cpu_start);
  cpu_start->tms_utime -= (delta[1] * clk_tck);
}



/*** RestoreLog: restores the .log (and the .llog files) upon restart ******
 ***************************************************************************/

void RestoreLog(void)
{
  char   *shell_cmd;                             /* used by 'system' below */
  char   *outfile;                           /* temporary output file name */

  FILE   *logptr;                                     /* .log file pointer */
  FILE   *outptr;

  char   *logline;                                /* array of read buffers */
  long   saved_count_tau;          /* count_tau as read from the .log file */
  long   max_saved_count;    /* last count_tau that was saved in .log file */
  long   i;                                                /* loop counter */


/* this is the last line we've written into the .log file */

  max_saved_count = state->tune.initial_moves+proc_init+count_tau*proc_tau;

  logline   = (char *)calloc(MAX_RECORD, sizeof(char));
  shell_cmd = (char *)calloc(MAX_RECORD, sizeof(char));

/* get temporary file name for output file */


    outfile   = (char *)calloc(MAX_RECORD, sizeof(char));
    outfile   = strcpy(outfile,"logXXXXXX");      /* required by mkstemp() */
    if ( mkstemp(outfile) == -1 )         /* get unique name for temp file */
      error("RestoreLog: error creating temporary (log) file");
   
/* restore the global .log file */

    saved_count_tau = -1;                             /* reset the counter */

/* open the log file and a temporary output file */

    logptr = fopen(logfile, "r");
    if ( !logptr )
      file_error("RestoreLog (at open log file for reading)");

    outptr = fopen(outfile, "w");
    if ( !outptr )
      file_error("RestoreLog (at open temp log file for writing)");

/* read and write the first few title and caption lines */

    for (i=0; i<4; i++) {
      if ( NULL == fgets(logline, MAX_RECORD, logptr ) )
	error("RestoreLog: error reading log file captions");
      fprintf(outptr, "%s", logline);
    }
    
/* read and write the actual log lines till we are at current time */

    while ( (saved_count_tau < max_saved_count) && 
	    (NULL != fgets(logline, MAX_RECORD, logptr)) ) {
      if ( 1 != sscanf(logline, "%ld", &saved_count_tau) )
	error("RestoreLog: error reading saved_count_tau (after %d)", 
	      saved_count_tau);
      fprintf(outptr, "%s", logline);
    } 

    fclose(logptr);
    fclose(outptr);

/* rename tmpfile into new file */

    sprintf(shell_cmd, "cp -f %s %s", outfile, logfile);
    
    if ( -1 == system(shell_cmd) )
      error("RestoreLog: error renaming temp file %s", outfile);
    
    if ( remove(outfile) )
      warning("RestoreLog: temp file %s could not be deleted", 
	      outfile);


}





/*** FUNCTIONS WHICH WRITE LOG FILES ***************************************/

/*** WriteLandscape: writes iterations, temperature, dS/S, energy, delta_energy,  **
 *              mean, std deviation, estimate_mean, estimate_sd, acceptance ratio*
 *              to look at the landscape of the problem Will be either  **
 *              the entire landscape or only the accepted landscape *
 ***************************************************************************/

void WriteLandscape(char *landfile, int iteration, double delta_energy)
{
  const char *format =
    "  %9d %14.6f  %10.6e %16.6f %16.6f %5.2f \n";

  FILE *landptr;

    landptr = fopen(landfile, "a");   /* first  open appropriate landscape file */
    if ( !landptr ) 
      file_error("WriteLandscape");
    fprintf(landptr, format, 
	    iteration, 
	    1.0/S, 1.0/dS, energy, delta_energy,
	    acc_ratio);
    fclose( landptr );

}


/*** WriteLog: writes things like mean and variation, Lam estimators, dS, **
 *             alpha and acceptance ratio to the log files and to stdout   *
 *             (if -l is chosen).                                          *
 ***************************************************************************/

void WriteLog(void)
{
  FILE   *logptr;                      /* file pointer for global log file */

    logptr = fopen(logfile, "a");   /* first write to the global .log file */
    if ( !logptr ) 
      file_error("WriteLog");
    PrintLog(logptr, 0);
    fclose( logptr );
  
    if ( log_flag ) {                        /* display log to the screen? */
      PrintLog(stdout, 0);
      fflush( stdout );
    }
    
     
 
}


/*** PrintLog: actually prints the log to wherever it needs to be printed **
 ***************************************************************************/

void PrintLog(FILE *outptr, int local_flag)
{
  if (r_flag){
    const char *format="  %10ld %14.6f  %10.6e %16.6f %16.6f %16.6f %16.6f %5.2f %8.5f %1.8f %1.8f\n";
    if ( count_tau % (print_freq * captions) == 0 ) {
      fprintf(outptr, 
        "\n iterations              T          dS/S            meanE");
      fprintf(outptr, 
        "              sdE         (e)meanE           (e)sdE");
      fprintf(outptr, 
        "   acc    alpha  r_estimate r_n_estimate\n\n");
    }
                                                               /* print data */
      fprintf(outptr, format, 
        (state->tune.initial_moves+proc_init+count_tau*proc_tau), 
        1.0/S, dS/S, 
        mean, sqrt(vari), estimate_mean, estimate_sd, 
        acc_ratio, alpha, r_est, r_n_est);
  }


  else{
    const char *format =
      "  %10ld %14.6f  %10.6e %16.6f %16.6f %16.6f %16.6f %5.2f %8.5f\n";

    if ( count_tau % (print_freq * captions) == 0 ) {
      fprintf(outptr, 
        "\n iterations              T          dS/S            meanE");
      fprintf(outptr, 
        "              sdE         (e)meanE           (e)sdE");
      fprintf(outptr, 
        "   acc    alpha\n\n");
    }
                                                               /* print data */
      fprintf(outptr, format, 
        (state->tune.initial_moves+proc_init+count_tau*proc_tau), 
        1.0/S, dS/S, 
        mean, sqrt(vari), estimate_mean, estimate_sd, 
        acc_ratio, alpha);
  }
}

/*Mix Log */
/* writes a detailed log giving information about each nodes energy, 
 * probability of being chosen, and dance_partner at each mix interval.       */


void WriteMixLog(double *node_prob, int* dance_partner)
{


}

void shift_left(char **buff, int length_buff, int length_loop, int *set_sizes, int start_point){
  int displacements = 0;
  int corrected_start_point = start_point;
  int i;
  int items_shifted=0;
  /*first we shift the start point to the first '\0' */
  while ((buff[corrected_start_point][0]!='\0') && (items_shifted<length_buff)){
    ++corrected_start_point;
    ++items_shifted;
  }
  i=corrected_start_point;
  while(items_shifted<length_buff && i<length_loop){
    if (strcmp(buff[i], "\0")){
      set_sizes[i-displacements]=set_sizes[i];
      free(buff[i-displacements]);
      buff[i-displacements]=(char*)malloc(set_sizes[i]*sizeof(char));
      strncpy(buff[i-displacements], buff[i], set_sizes[i]);
    }
    else
      ++displacements;
    ++i;
  }
}

void copy_dyn_int_array_to_const_int_array(int dest[], int *src, int size){
  int i;
  for (i=0; i<size; ++i)
    dest[i]=src[i];
    
}


char** extract_string_set_from_char_list(char* list, int length_list, int *length_set, int *sizes, int *displacements, int *set_sizes){
  char **buff_list;
  char **set;
  char *curr_string;
  int  BuffNotEmpty = 1;
  int  length_buff=length_list;
  int  loop_length=length_buff;
  int i;
  *length_set=0;
  /*first, buff is a copy of the list and set_sizes is a copy of sizes*/
  buff_list=(char **)malloc(length_list*sizeof(char*));
  for (i=0; i<length_list; ++i){
    buff_list[i]=(char *)malloc(sizes[i]*sizeof(char));
    strncpy(buff_list[i], list+displacements[i], sizes[i]);
    set_sizes[i]=sizes[i];
  }
  /*Next, we extract the first element from the buffer, increment the length of the final set,
   * and erase the first element of the buffer. */
  curr_string=buff_list[0];
  *length_set=*length_set + 1;
  while (BuffNotEmpty){
    for (i=*length_set; i<loop_length; ++i){
      /*erase all elements that are the same as curr_string */
      if (strcmp(curr_string, buff_list[i])==0){
        length_buff--;
        buff_list[i][0]='\0';
      }
    }
    if (length_buff==*length_set)
      BuffNotEmpty=0;
    else
      shift_left(buff_list, length_buff, loop_length, set_sizes, *length_set);

    loop_length=length_buff;
    curr_string=buff_list[*length_set];
    *length_set=*length_set + 1;
    }
  /*The while loop increments *length_set 1 more time than it should. */
  *length_set=*length_set-1;
  set=(char**)malloc(sizeof(char*)*(*length_set));
  for (i=0; i<*length_set; ++i){
   set[i]=(char*)malloc(sizeof(char)*set_sizes[i]);
   strncpy(set[i],buff_list[i],set_sizes[i]);
 }
  for (i=1; i<length_list; ++i)
    free(buff_list[i]);
  free(buff_list);
  return set;
}


