/******************************************************************
 *                                                                *
 *   tsp_sa.c                                                     *
 *                                                                *
 ******************************************************************
 *                                                                *
 * This is the file which specifies the user interface of the TSP * 
 * code. It contains functions for file I/O and  command line     *
 * parsing a la Yogi for the tsp_sa program                       *
 *                                                                *
 ******************************************************************/
/**********************************************************************
* File Name: tsp_sa.c Created 03-00 by: Lorraine Greenwald            *
*   -D landscape option by Lorraine Greenwald Oct 2002                *
* LG 03-03: use unsigned shorts for tour arrays                       *
*           and floats rather than doubles for edges                  * 
*           for  BIG tsp instances memory crunch                      *
* LG 08-02: incorporate lsa.c and sa.h version 9.3.1                  *
* LG 08-01: use distribution type from input file                     *
*           do not print problem instance in output file              *
* LG 03-02: need q for gen visiting distribution input file           *
* LG 05-02: set factors only dependent on q for general visiting      *
* distribution input file here: call qgt2_init qlt2_init              *
*                               from distributions.c                  *
*                                                                     *
* Problem Specific 'main' routines for tsp                            *
* i.e. these are routines that would be in a Main progam              *
* opens files; initializes (InitialMove); and finishes (FinalMove)    *
*                                                                     *
* This file has routines that:                                        *
* read LAM tune parameters from input file (ReadAParameters)          *
* read problem specific parms from file (InitTSP move.c)              *
* Calls routines that:                                                *
* init stats (InitMove from move.c)                                   *
* set up starting tour (start_tour from move.c)                       *
*       calculate cost for entire tour (tour_cost from move.c)        *
*       sets initial global min tour and cost                         *
* allocates tour memory (called from InitTSP);                        *
* deallocates tour memory (called from FinalMove)                     *
***********************************************************************/
#ifdef ICC
#include <mathimf.h>
#else
#include <math.h>
#endif
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>


#include "error.h"
#include "edge_wt.h" /* prototypes for edge handling routines */
#include "move.h"   /* prototypes for move routines and AP*/
#include "distributions.h"   /* DistP.variables and prototypes */
#include "sa.h"   /* generic lam parms and problem specific state variables */
#include "initialize.h"
#include "random.h"
/* also includes generic prototypes for app specific move routines */

#ifndef MPI_INCLUDED
#include <mpi.h>
#include "MPI.h"
#endif

#define  OPTS       ":b:c:C:e:Ef:hlLNpQrStTvw:W:y:"
                                             /* command line option string */
                                             /* D will be debug, like fly */
                     /* must start with :, option with argument must have a : following */
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/


/* Help, usage and version messages */

#ifdef MPI
static const char usage[]    =
"Usage: tsp_sa.mpi [-C <covar_ind>]  [-e <freeze_crit>]\n"
"                 [-E] [-f <param_prec>] [-h] [-l] [-L] [-N] [-p] [-r]\n"
"                 [-S] [-t] [-T] [-v] [-w <outfile> ] [-W <tune_stat>]\n"
"                 [-y <log_freq> ] <infile> \n";
#else
static const char usage[]    =
"Usage: tsp_sa [-e <freeze_crit>] [-E ] [-f <param_prec>]\n"
"             [-h] [-l] [-p] [-Q] [-N ] [-r] [-t] [-v] [-w <outfile>]\n"
"             [-y <log_freq>] <infile> \n";
#endif

static const char help[]     =

#ifdef MPI
"Usage: tsp_sa.mpi [options] <infile>\n\n"
#else
"Usage: tsp_sa [options] <infile>\n\n"
#endif

"Arguments:\n"
"  <infile>            input data file\n\n"

"Options:\n"
#ifdef MPI
"  -C <covar_ind>      set covar sample interval to <covar_ind> * tau\n"
#endif
"  -e <freeze_crit>    set annealing freeze criterion to <freeze_crit>\n"
"  -E                  run in equilibration mode\n"
"  -f <param_prec>     float precision of parameters is <param_prec>\n"
"  -h                  prints this help message\n"
"  -l                  echo log to the terminal\n"
#ifdef MPI
"  -L                  write local logs (llog files)\n"
#endif
"  -N                  generates landscape to .landscape file in equilibrate mode\n"
"  -p                  prints move acceptance stats to .prolix file\n"
"  -r                  tweak coordinates in random order\n"
#ifndef MPI
"  -Q                  quenchit mode, T is lowered immediately to zero\n"
#else
"  -S                  disable tuning stop flag\n"
#endif
"  -t                  write timing information to .times file\n"
#ifdef MPI
"  -T                  run in tuning mode\n"
#endif
"  -v                  print version and compilation date\n"
"  -w                  print output and log file to x <outfile> and <outfile>.log\n"
#ifdef MPI
"  -W <tune_stat>      tuning stats written <tune_stat> times per interval\n"
#endif
"  -y <log_freq>       write log every <log_freq> * tau moves\n\n"

"Please report bugs to <lorraine@odd.bio.sunysb.edu>. Thank you!\n";

static char version[MAX_RECORD];                 /* version gets set below */

static const char verstring[] = 

"compiled by:      %s\n"
"         on:      %s\n"
"      using:      %s\n"
"       date:      %s at %s\n";

/* other static variables */

static char   *inname;                   /* filename of input file */
static char   *outname;                  /* filename of output file */
static char   *param_inname;             /* filename of input parameter file */
static char   *instance_inname;          /* filename of tsp instance file */


static int    argcsave;                 /* static variable for saving argc */
static char   *argvsave;               /* static variable for saving argv */
 
static int    precision   = 6;                    /* precision for eqparms */
static int    prolix_flag = 0;               /* to prolix or not to prolix */
static int    landscape_flag = 0;            /* generate energy landscape data */
static int    diff_outfile = 0;  /* set to 1 if the output file has a diff na- *
                                  * me than input file when using -w           */
           /* set the landscape flag (and the landscape filename) in lsa.c */


/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

static const char prog_buf[]="@(#)tsp_sa version3.1 send";
/* enables us determine without a doubt which version was used */ 
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/


/*** FUNCTIONS *************************************************************/

/*** COMMAND LINE OPTS ARE ALL PARSED HERE *********************************/

/*** ParseCommandLine: well, parses the command line and returns an index **
 *                     to the 1st argument after the command line options  *
 ***************************************************************************/
 
int ParseCommandLine(int argc, char **argv)
{
  int           c;                   /* used to parse command line options */
  
/* external declarations for command line option parsing */

  extern char *optarg;                     /* command line option argument */
  extern int  optind;                /* pointer to current element of argv */
  extern int  optopt;               /* contain option character upon error */
          /* the above three variables are declared in #include <unistd.h> */

/* set the version string */
#ifdef MPI    
      /*sprintf(version, "tsp_sa version %s parallel", VERS);*/
#else
        /* sprintf(version, "tsp_sa version %s serial", VERS); */
#endif


/* following part sets default values for command line options */
/* Edited by Seb, Aug 30 2024, to have more reasonable defaults for
 * print_freq and state_write */
  captions        = 100000000;   /* default freq for writing captions (off) */
  print_freq      = 100;            /* default freq for writing to log file */
  state_write     = 1000;         /* default freq for writing to state file */

  stop_flag       = absolute_freeze;             /* type of stop criterion */
  time_flag       = 0;                                  /* flag for timing */
  log_flag        = 0;              /* flag for writing logs to the screen */

  random_tweak    = 1;                      /* random tweak is the default */

#ifdef MPI
  tuning          = 0;                    /* tuning mode is off by default */
  covar_index     = 1;      /* covariance sample will be covar_index * tau */
  write_tune_stat = 1;         /* how many times do we write tuning stats? */
  auto_stop_tune  = 1;               /* auto stop tuning runs? default: on */
  write_llog      = 0; /* write local llog files when tuning; default: off */
#endif

/* following part parses command line for options and their arguments      */
  
  optarg = NULL;
  while( (c = getopt(argc, argv, OPTS)) != -1) {
    switch(c) {
    case 'b':            /* -b sets backup frequency (to write state file) */
/*      state_write = strtol(optarg, NULL, 0);
      if ( state_write < 1 )
	error("tsp_sa: max. backup frequency is every tau steps i.e. -b 1");
      if ( state_write == LONG_MAX )
        error("tsp_sa: argument for -b too large"); */
      warning("tsp_sa: state files not supported!");
      break;
    case 'c':               /* -c sets the frequency for printing captions */
      error("tsp_sa: -c is not supported anymore, captions are off for good");
/* if you want to be able to insert captions into .log files, uncomment    *
 * the following lines; CAUTION: make sure that RestoreLog() works proper- *
 * ly with captions before you do this                                     */
/*      captions = strtol(optarg, NULL, 0);
      if ( captions < 1 )
        error("tsp_sa: can't print captions more than every line (-c 1)");
      if ( captions == LONG_MAX )
        error("tsp_sa: argument for -c too large");                         */
      break;
    case 'C':             /* parallel code: -C sets the covar sample index */
#ifdef MPI
/* -C is currently disabled since I (Yogi) could not get it to work and I  *
 * got tired of debugging and I didn't really need it at that time and so  *
 * on and so forth; try running on a smaller number of nodes or get it to  *
 * work yourself; that's what I say. Grmbl.                                */
/*    covar_index = atoi(optarg); 
      if ( covar_index < 1 )
        error("tsp_sa: covariation sample index must be >= 1");            */
      error("tsp_sa: -C does not work yet, try to run on fewer nodes");
#else
      error("tsp_sa: can't use -C in serial, tuning only in parallel");
#endif
      break;
    case 'e':                            /* -e sets the stopping criterion */
      if( !(strcmp(optarg,"pfreeze")) ) 
	stop_flag = proportional_freeze;  
      else if( !(strcmp(optarg,"afreeze")) )
	stop_flag = absolute_freeze;
      else if( !(strcmp(optarg,"abs")) )
	stop_flag = absolute_energy;
      else
	error("tsp_sa: valid stopping criteria are pfreeze, afreeze, abs");
      break;
    case 'E':                                /* -E does equilibration runs */
      equil = 1;
      break;
    case 'f':
      precision = atoi(optarg);           /* -f determines float precision */
      if ( precision < 0 )
	error("tsp_sa: what exactly would a negative precision be???");
      if ( precision > MAX_PRECISION )
	error("tsp_sa: max. float precision is %d!", MAX_PRECISION);
      break;
    case 'h':                                            /* -h help option */
      PrintMsg(help, 0);
      break;
    case 'l':                         /* -l displays the log to the screen */
      log_flag = 1;
      break;
    case 'L':                   /* -L writes local .llog files when tuning */
#ifdef MPI
      write_llog = 1;
#else
      error("tsp_sa: can't use -L in serial, tuning only in parallel");
#endif
      break;
    case 'N':                  /* -N sets laNdscape flag and Equilibrate mode */
       equil = 1;
       landscape_flag = 1;
      break;        
    case 'p':                                   /* -p sets prolix printing */
      prolix_flag = 1;
      break;
    case 'Q':                  /* -Q sets quenchit mode (serial code only) */
#ifdef MPI
      error("tsp_sa: can't use -Q for parallel code");
#else
      quenchit = 1;
#endif
      break;
    case 'r':
      random_tweak = 0;
      break;
    case 'S':                         /* -S unsets the auto_stop_tune flag */
#ifdef MPI
      auto_stop_tune = 0;
#else
      error("tsp_sa: can't use -S in serial, tuning only in parallel");
#endif   
      break;
    case 't':                   /* -t saves times in data and .times files */
      time_flag = 1;
      break;
    case 'T':                       /* -T does tuning (parallel code only) */
#ifdef MPI
      tuning = 1;
#else
      error("tsp_sa: can't use -T in serial, tuning only in parallel");
#endif
      break;
    case 'v':                                  /* -v prints version message */
     /* fprintf(stderr, "%s\n", version);
      fprintf(stderr, verstring, USR, MACHINE, COMPILER, __DATE__, __TIME__); 
      exit(0);*/
      break;
    case 'w': /* -w changes the output file */
      diff_outfile=1;
      outname=(char *)calloc(MAX_RECORD + 1, sizeof(char));
      strcpy(outname, optarg);
      break;
    case 'W':       /* -W determines the frequency of writing tuning stats */
#ifdef MPI
      write_tune_stat = atoi(optarg);
      if ( write_tune_stat < 1 )
	error("tsp_sa: frequency of writing tune stats must be >= 1");
#else
      error("tsp_sa: can't use -W with serial code, tuning only in parallel");
#endif
      break;
    case 'y':                             /* -y set frequency to print log */
      print_freq = strtol(optarg, NULL, 0);
      if ( print_freq < 1 )
	error("tsp_sa: can't print status more than every tau steps (-y 1)");
      if ( print_freq < 1000 )
#ifdef MPI
	if ( myid == 0 )
	  warning("tsp_sa: writing log often (< 1000) slows down your code!");
#endif
      if ( print_freq == LONG_MAX )
	error("tsp_sa: argument for -y too large");
      break;
    case ':':
      error("tsp_sa: need an argument for option -%c", optopt);
      break;
    case '?':
    default:
      error("tsp_sa: unrecognized option -%c", optopt);
    }
  }

/* error checking here */

#ifdef MPI
  if ( (tuning == 1) && (equil == 1) )
    error("tsp_sa: can't combine -E with -T");
  if ( write_llog && !tuning )
    error("tsp_sa: -L only makes sense when tuning");
#else
  if ( (quenchit == 1) && (equil == 1) )
    error("tsp_sa: can't combine -E with -Q");
#endif
  /* Modified this to allow inname=outname */
  if ( ((argc - (optind - 1)) != 2) ) 
    PrintMsg(usage, 1);

/* set output file name here so that lsa.c stores the correct .log file */
  if(!diff_outfile){
    outname = (char *)calloc(MAX_RECORD + 1, sizeof(char));   /* output file */
    sprintf(outname, "%s.output", argv[optind]); 
    SetOutname(inname);                    /* communicates outname to lsa.c */
  }
  else{
    SetOutname(outname);
  }
  /*come back later */
  
  return optind;
} /* end ParseCommandLine*/



/*** FUNCTIONS THAT INITIALIZE/RESTORE THINGS ******************************/

/**********************************************************************
 * This routine reads the distribution parameters  from the input file*
 * and stores them in DistP.xx from distributions.h                   *
 * LG 03-02: need q for gen visiting distribution input file          *
 * LG 05-02: set factors only dependent on q for general visiting     *
 * distribution by calling qgt2_init or glt2_init from distributions.c*
 **********************************************************************/

void InitDistribution(FILE *fp)
{ /* begin init distribution */

   /* read in the data from file */
   /* read in header info then output it*/
   /* At the end of file exit from the loop */


  fp = FindSection(fp,"distribution_parameters"); /* find distribution section */
  if( !fp )
    error("ReadTune: could not locate distribution_parameters\n");

  fscanf(fp,"%*s\n");    /* advance past title line no blanks allowed! */

                                                  /* read distribution stuff */
   if ( 2!= (fscanf(fp,"%d %lf\n",  &(DistP.distribution), &(DistP.q) )) )
      error ("ReadTune: error reading distribution stuff"); 
 

	if ((DistP.distribution > 11) || (DistP.distribution < 1))
	  {error ("tsp_sa: distribution must be int from 1 to 11 \n");}
	else if 
	  ((DistP.distribution == 5)||(DistP.distribution == 8)||(DistP.distribution == 10))
	  {error ("tsp_sa: chose a distribution that returns positive values \n");
	  }
	else if (DistP.distribution == 7 )
          {  
	    if ((DistP.q >= 3.0) || (DistP.q <= 1.0))
	      {error ("tsp_sa: q must be between 1 and 3 \n");}
	    else if (DistP.q == 2.0)
	      {DistP.distribution = 4;
	      /* tsp needs abs lorentz, fly should be 5*/
	      printf ("tsp_sa: q=2 is lorentz--setting distribution to 4\n");}
	    else if (DistP.q > 2.0)
	      {qgt2_init(); }
	    else 
	      {qlt2_init(); }
            }
	    /***************LG 05-02***************/ 
	    /* calculate q dependent factors that */
	    /* do not change for entire run.      */
	    /* these live in distributions.h      */
	    /***************LG 05-02***************/ 
}  /* end InitDistribution */




/***************************************************************************
 * LG 03-03 use two files params and instance. Needed for BIG problems     *
 * LG 08-02: initial_move becomes InitialMove for 9.3.1.compatibility      *
 *                - various static things for filenames and such           *
 *                - Lam parameter struct for use in lsa.c (from tune sect) *
 *                - cost function in move.c (InitTSP))                     *
 *                - move generation in move.c (InitMoves())                *
 *                - sets initial energy and checks validity of initial     *
 *                  parameters according to limit ranges                   *
 *                then it returns the initial temperature to the caller    *
 * LG 08-01: Killed envp                                                   *
* initial move is the Lam plug-in that:                                    * 
* opens the input file,                                                    *
* *** should open output file, but initialize(lsa.c) does it  ***          * 
* calls InitTSP(move.c) - reads problem specific data, and gets intial tour*
* calls InitMoves (move.c) - stats=0, gets aparms & randomizes xsubj       *
* calculates initial energy and returns it as a parameter                  *
* sets global minimum variables;                                           *
* returns initial temperature (ap.start_tempr)                             *
***************************************************************************/

double InitialMove (int argc, char ** argv, int opt_index,
                    NucStatePtr state_ptr, double *p_chisq)
{ /* BEGIN InitialMove */

  FILE    *param_infile;
  FILE    *instance_infile;
  
  double  i_temp=0;

  char    *p;
  int     i;
  
  SAType  in_tune;          /* temporary struct to read in tune parameters */

/* save some things in static storage for FinalMove */
  inname  = (char *)calloc(MAX_RECORD + 1, sizeof(char)); /* input file name*/
  param_inname  = (char *)calloc(MAX_RECORD + 1, sizeof(char)); /* parameter file */
  instance_inname  = (char *)calloc(MAX_RECORD + 1, sizeof(char));  /* instance file */

  inname  = strcpy(inname, argv[opt_index]);

  sprintf(param_inname, "%s.params", inname);
  sprintf(instance_inname, "%s.instance", inname);

  argvsave = (char *)calloc(MAX_RECORD, sizeof(char));
  for (i=0; i < argc; i++){
    if(i>0)
      argvsave=strcat(argvsave, " ");
    argvsave=strcat(argvsave, argv[i]);
  }

/* set the prolix flag (and the .prolix filename) in move.c, if necessary */

  if ( prolix_flag )
    SetProlix(prolix_flag, outname, 1);       /* 1 means: new .prolix file */

/* set the landscape flag (and the landscape filename) in lsa.c */

  if ( landscape_flag )
    InitLandscape(landscape_flag, outname);       /*  lives in lsa.c */
  param_infile = fopen(param_inname, "r");
   if ( !param_infile ) {
     perror("tsp_sa");
     exit(1);
   }

  in_tune  = ReadTune(param_infile);   /* read tune_parameter section */  
/* initialize Lam parameters (see sa.h for further detail) */
    state_ptr->tune.lambda              = in_tune.lambda;
    state_ptr->tune.lambda_mem_length_u = in_tune.lambda_mem_length_u;
    state_ptr->tune.lambda_mem_length_v = in_tune.lambda_mem_length_v;
    state_ptr->tune.control             = in_tune.control;
    state_ptr->tune.initial_moves       = in_tune.initial_moves;
    state_ptr->tune.tau                 = in_tune.tau;
    state_ptr->tune.freeze_count        = in_tune.freeze_count;
    state_ptr->tune.update_S_skip       = in_tune.update_S_skip;
    state_ptr->tune.criterion           = in_tune.criterion;
#ifdef MPI
    state_ptr->tune.mix_interval        = in_tune.mix_interval;
    state_ptr->tune.glob_interval       = in_tune.glob_interval;
    state_ptr->tune.ngroups             = in_tune.ngroups;
    state_ptr->tune.score_method        = in_tune.score_method;
    ngroups=in_tune.ngroups;
    score_method=in_tune.score_method;
    glob_interval=in_tune.glob_interval;
#endif
    AssignGroups();
  
/* initialize some Lam/Greening structures */
  p = state_ptr->tune.progname;     /* tune.progname contains program name */
  /*p = strcpy(p, version);*/

  state_ptr->tune.debuglevel = 0;          /* following stuff not used now */
  p = state_ptr->tune.tunename;
  p = strcpy(p, "Squelch the Weasel");    /* Ween tune, what else? */ 
                                          /*Oye Vey! et tu, Yogi?*/
  state_ptr->tune.tunefile = NULL;
/* read data file and initialize different things for cost function and    */
/* move generation (in move.c)                                             */

 
  instance_infile = fopen(instance_inname, "r");
  if ( !instance_infile ) {
    perror("tsp_sa");
    exit(1);
  }

                    /* read input from tsp lib file */
	           /* generate initial tour and get initial energy */
        InitTSP(instance_infile); /* this routine is in move.c */
        *p_chisq = StartTour(); /* this routine is in move.c */

 /* set initial temperature and initialize random number generator and *
  * annealing parameters                                               */
	i_temp   = InitMoves(param_infile, state_ptr->tune.tau); 
        InitDistribution(param_infile);   /* initialize distribution stuff */
        fclose(param_infile);
        fclose(instance_infile);

  return i_temp;
}

/***************************************************************************
 * functions that communicate with savestate.c.                            *
 * *************************************************************************/

/* GetOptions: returns command line options to savestate.c. For detailed 
 * meanings of all these options see ParseCommandLine().
 * Opts struct defined in moves.h */

Opts *GetOptions(void)
{
  Opts *options;
  options = (Opts *)malloc(sizeof(Opts));
  options->diff_outfile = diff_outfile;
  options->inname = inname;
  options->outname = outname;
  options->argv = argvsave;
  options->stop_flag = stop_flag;
  options->prolix_flag = prolix_flag;
  options->landscape_flag = landscape_flag;
  options->time_flag = time_flag;
  options->log_flag = log_flag;
  options->print_freq = print_freq;
  options->captions = captions;
  options->precision = precision;
  options->quenchit = quenchit;
  options->equil = equil;
  options->state_write=state_write;
#ifdef MPI
  options->tuning = tuning;
  options->covar_index = covar_index;
  options->write_tune_stat=write_tune_stat;
  options->auto_stop_tune=auto_stop_tune;
#endif
  return options;
}
/*** RestoreOptions: restores the values of the command line opt variables *
 *                   from the Opts struct (used for restoring a run)       *
 ***************************************************************************/

void RestoreOptions(Opts *options)
{

/* restore input/output file names and the full command line string; note  *
 * that the output file name needs to be communicated to lsa.c, since we   * 
 * need it there for setting up log and tuning file names                  */
  inname = options->inname;
  if (!options->diff_outfile){
    sprintf(outname, "%s.output", inname); 
    SetOutname(inname);
  }
  else if (strcmp(outname, "(null)")==0){ /*i.e., outname == "(null)" */
    sprintf(outname, "%s.output", inname); 
    SetOutname(inname);
  }
  else{
    outname = options->outname;
    SetOutname(outname);
  }
  argvsave = options->argv;

  stop_flag   = options->stop_flag;

  prolix_flag = options->prolix_flag;
  if ( prolix_flag )
    SetProlix(prolix_flag, outname, 0);

  landscape_flag = options->landscape_flag;
  if ( landscape_flag )
    error("RestoreOptions: cannot restore an equilibration (lanDscape) run");

  log_flag    = options->log_flag;
  time_flag   = options->time_flag;

  state_write = options->state_write;
  print_freq  = options->print_freq;
  captions    = options->captions;

  precision   = options->precision;

  quenchit    = options->quenchit;
  equil       = options->equil;
  if ( equil )
    error("RestoreOptions: cannot restore an equilibration run");

#ifdef MPI
  tuning          = options->tuning;
  if ( tuning )
    error("RestoreOptions: cannot restore a tuning run");
  covar_index     = options->covar_index;
  write_tune_stat = options->write_tune_stat;
  auto_stop_tune  = options->auto_stop_tune;
#endif

  free(options);
}





/*** RestoreState: called when an interrupted run is restored; does the ****
 *                 following:                                              *
 *                 - stores various static things for filenames and such   *
 *                 - initializes Lam parameters for lsa.c                  *
 *                 - initializes model and scoring funcs & solver stepsize *
 *                 - initializes move generation in moves.c                *
 *                 - restores state at which previous run was interrupted  *
 *                   according to state file                               *
 *                                                                         *
 * Comment by JR:  RestoreState was originally going to be implemented with*
 * branches in InitialMove. That won't work because when when this func.   *
 * returns, control should go right to Loop(), skipping all the initiali-  *
 * zation stuff in Initialize. Hence most of the code in InitialMove is    *
 * just repeated here.                                                     *
 ***************************************************************************/


void RestoreState(char *statefile, NucStatePtr state_ptr, double *p_chisq)
{
  FILE           *infile;                /* pointer to original input file */
  FILE           *instance_infile, *param_infile;
  char           *instance_inname, *param_inname;
  char           *p;                                   /* temporary string */
  double         i_temp;
  SAType         in_tune;          /* temporary struct for tune parameters */

  Opts           *options;         /* used to restore command line options */
  MoveState      *move_ptr;                       /* used to restore moves */
  double         *stats;                      /* used to restore Lam stats */
  unsigned short *rand;                         /* used to restore ERand48 */
  double         delta[2];                        /* used to restore times */
  double         _temp;

/* allocate memory for structures that will be returned *
 * No, stats does not get allocated in StatRead.        */
  options = (Opts *)malloc(sizeof(Opts));
  options->inname    = (char *)calloc(MAX_RECORD, sizeof(char));
  options->outname   = (char *)calloc(MAX_RECORD, sizeof(char));
  options->argv      = (char *)calloc(MAX_RECORD, sizeof(char));
  instance_inname    = (char *)calloc(MAX_RECORD, sizeof(char));
  param_inname       = (char *)calloc(MAX_RECORD, sizeof(char));

  stats    = (double *)calloc(31, sizeof(double));
  move_ptr = (MoveState *)malloc(sizeof(MoveState));
  rand     = (unsigned short *)calloc(3, sizeof(unsigned short));

  StateRead(statefile, options, move_ptr, stats, rand, delta);
  
  RestoreOptions(options);
  sprintf(param_inname, "%s.params", inname);
  sprintf(instance_inname, "%s.instance", inname);
  
  p = state_ptr->tune.progname;     /* tune.progname contains program name */
  /*p = strcpy(p, version);*/

  state_ptr->tune.debuglevel = 0;          /* following stuff not used now */
  p = state_ptr->tune.tunename;
  p = strcpy(p, "Squelch the Weasel");    /* Ween tune, what else? */ 
                                          /*Oye Vey! et tu, Yogi?*/
  state_ptr->tune.tunefile = NULL;

/* read data file and initialize different things for cost function and    */
/* move generation (in move.c)                                             */

  param_infile = fopen(param_inname, "r");
  if ( !param_infile ) {
    perror("tsp_sa");
    exit(1);
  }
  instance_infile = fopen(instance_inname, "r");
  if ( !instance_infile ) {
    perror("tsp_sa");
    exit(1);
  }
 InitTSP(instance_infile);
 in_tune = ReadTune(param_infile);
 InitDistribution(param_infile);
 state_ptr->tune.lambda              = in_tune.lambda;
  state_ptr->tune.lambda_mem_length_u = in_tune.lambda_mem_length_u;
  state_ptr->tune.lambda_mem_length_v = in_tune.lambda_mem_length_v;
  state_ptr->tune.control             = in_tune.control;
  state_ptr->tune.initial_moves       = in_tune.initial_moves;
  state_ptr->tune.tau                 = in_tune.tau;
  state_ptr->tune.freeze_count        = in_tune.freeze_count;
  state_ptr->tune.update_S_skip       = in_tune.update_S_skip;
  state_ptr->tune.criterion           = in_tune.criterion;
#ifdef MPI
  state_ptr->tune.mix_interval        = in_tune.mix_interval;
#endif

 _temp   = InitMoves(param_infile, state_ptr->tune.tau);
 fclose(param_infile);
 fclose(instance_infile);
     RestoreMoves(move_ptr);
    RestoreLamstats(stats);
  if (time_flag)
    RestoreTimes(delta);
    InitERand(rand);
  if( prolix_flag )
    RestoreProlix();
  free(param_inname);
  free(instance_inname);
}


/*** THE FINAL MOVE FUNCTION ***********************************************/


/*****************************************************
* Final move is problem specific finish routine.     *
* Calculates and  prints the final tour cost;        *
* deallocates tour memory.                           *
******************************************************/
void FinalMove( void )

{  /* begin final_move*/
  FILE     *outfile;                             /* pointer to output file */
  int      i = 0;                           /* loop counter */
  AParms   ap;           /* annealing parameter stuct for annealing output */
  double   equil_var[2];         /* array for results of equilibration run */


#ifdef MPI
  int      winner = 0;                           /* id of the winning node */
  double   minyet = DBL_MAX;         /* minimum score, used to find winner */

  double   *final_e;               /* array of final energies of all nodes */
  
  final_e = (double *)calloc(nnodes, sizeof(double));
#endif

/* get the final energies and iterations from move.c */
  ap = GetFinalInfo();

/* for equilibration runs: get the final equilibration results */
  if ( equil )
    GetEquil(equil_var);
#ifdef MPI
/* Seb: Might change this to use MPI_Reduce instead of this interesting choice *
 * We shall see. */

/* parallel code: find the node with the lowest energy */

/* Seb RV 12 Aug 2024: the code below is theoretically useless since calloc  *
 * initializes to 0 for us. Additionally, MPI_Allgather overwrites all of    *
 * this anyway                                                               */
//  for(i=0; i<nnodes; i++)                   /* initialize the energy array */
//    final_e[i] = 0;  
  
                                /* collect the final scores from all nodes */
  MPI_Allgather(&ap.stop_energy, 1, MPI_DOUBLE, final_e, 1, MPI_DOUBLE, 
		MPI_COMM_WORLD);

  for(i=0; i<nnodes; i++) {     /* evaluate the node with the lowest score */
    if ( final_e[i] <= minyet ) {
      minyet = final_e[i];
      winner = i;
    }
  }

/* write the answer to the output file */

  if ( myid == winner ) {
#endif           

    outfile = fopen(outname, "w");
    if ( !outfile ) {
      perror("tsp_sa");
      exit(1);
    }
    
    fprintf(outfile, "$version\n");
    fprintf(outfile, "%s\n", "WIP");
    for ( i=0; i<argcsave; i++ )
      fprintf(outfile, "%s ", argvsave);
    fprintf(outfile, "\n$$\n");
    
    if ( equil ) 
      PrintEquil(outfile, equil_var, "equilibrate_variance"); 
    else
      WriteResults(outfile, precision);
   
    fclose(outfile);
#ifdef MPI
  }
#endif
  if ( !equil && !nofile_flag )
#ifdef MPI
    if ( ! tuning )
#endif
      /*States Not Implemented at the moment - SRV 2025-11-16 */
      //StateRm();

/* free all memory */
	tour_deallocate();
  free(node_coords);
  
}  /* end FinalMove*/






/*** FILE FUNCTIONS FOR READING AND WRITING MISCELLANEOUS STUFF ************/


/*** InitEquilibrate: reads the equilibrate section of the data file, ******
 *                    which is needed for equilibration runs; then puts    *
 *                    the parameters into a static struct in lsa.c         *
 ***************************************************************************/

void InitEquilibrate(FILE *fp)
{
  ChuParam   l_equil_param;                    /* local equil_param struct */

  fp = FindSection(fp, "equilibrate");                /* find tune section */
  if( !fp )
    error("ReadEquilibrate: could not locate equilibrate section");
  
  fscanf(fp,"%*s\n");                         /* advance past title line 1 */
  
  if ( 1 != fscanf(fp, "%lf\n", &(l_equil_param.end_T)) )
    error("ReadEquilibrate: error reading equilibration temperature");

  fscanf(fp,"%*s\n");                         /* advance past title line 2 */

  if ( 1 != fscanf(fp, "%d\n", &(l_equil_param.fix_T_skip)) )
    error("ReadEquilibrate: error reading fixed temperature skip");

  fscanf(fp,"%*s\n");                         /* advance past title line 3 */

  if ( 1 != fscanf(fp, "%d\n", &(l_equil_param.fix_T_step)) )
    error("ReadEquilibrate: error reading fixed temperature step");

  SetEquilibrate(l_equil_param);
}



/*** ReadTune: reads the tune_parameters section in a data file and ********
 *             turns a SAType structure to the caller                      *
 ***************************************************************************/

void ReadGroupInfo(FILE *fp){
  char buff[200];
  fp = FindSection(fp, "tune_parameters");
  if (!fp)
    error("ReadTune: could not locate tune_parameters\n");
  for (int i=0; i<7; i++)
    fgets(buff, 200, fp);
  fscanf(fp,"%*d %*d %d", &ngroups);
}

SAType ReadTune(FILE *fp)
{
  SAType in_tune;              /* this is the SAType struct that we return */

  int intbuf1;         /* following four are used as temp. buffer for ints */
  int intbuf2;
  int intbuf3;
  int intbuf4;

  double dbuf1;     /* following four are used as temp. buffer for doubles */
  double dbuf2;
  double dbuf3;
  double dbuf4;

  unsigned short sbuf;



  fp = FindSection(fp, "tune_parameters");            /* find tune section */
  if( !fp )
    error("ReadTune: could not locate tune_parameters\n");

  fscanf(fp,"%*s\n");                         /* advance past title line 1 */

                                                  /* read Lam lambda stuff */
  if ( 4!= (fscanf(fp,"%lg %lg %lg %lg\n", &dbuf1, &dbuf2, &dbuf3, &dbuf4)) )
    error ("ReadTune: error reading Lam lambda stuff");

  in_tune.lambda              = dbuf1;
  in_tune.lambda_mem_length_u = dbuf2;
  in_tune.lambda_mem_length_v = dbuf3;
  in_tune.control             = dbuf4;

  fscanf(fp,"%*s\n");                         /* advance past title line 2 */

                        /* read initial moves, tau, freeze count and Sskip */
  if ( 4 !=(fscanf(fp,"%d %d %d %d\n",&intbuf1,&intbuf2,&intbuf3,&intbuf4)))
    error ("ReadTune: error reading line 2 of tune parameters");

  in_tune.initial_moves = intbuf1;
  in_tune.tau           = intbuf2;
  in_tune.freeze_count  = intbuf3;
  in_tune.update_S_skip = intbuf4;

  fscanf(fp,"%*s\n");                         /* advance past title line 3 */

                                                    /* read stop criterion */
  if ( 1 != (fscanf(fp,"%lg\n",&dbuf1)) )
    error ("ReadTune: error reading stop criterion");

  in_tune.criterion = dbuf1;

#ifdef MPI
  fscanf(fp,"%*s\n");                         /* advance past title line 4 */

  if ( 3 != fscanf(fp, "%d %d %d\n", &intbuf1, &intbuf2, &intbuf3) )
    error("ReadTune: error reading mixing interval");

  in_tune.mix_interval  = intbuf1;
  in_tune.glob_interval = intbuf2;
  in_tune.ngroups       = intbuf3;
 fscanf(fp, "%*s\n");
  if ( 1!= fscanf(fp, "%hu\n", &sbuf) )
    error("ReadTune: error reading score method");
  if (sbuf != 1){
    if (sbuf != 2)
      error("ReadTune: score method value was not recognised.");
  }
  in_tune.score_method=sbuf;
#endif

  in_tune.tunefile  = NULL;                  /* We're not using this stuff */
  in_tune.debuglevel = 0;

  return (in_tune);
}




/*** ReadAParameters: reads the AParm struct from an annealing_input sec- **
 *                    tion; these are the annealing parameters that are    *
 *                    not Lam-specific (and should NOT go into lsa.c       *
 ***************************************************************************/

AParms ReadAParameters(FILE *fp)
{
  AParms     l_aparms;                              /* local Aparms struct */

  fp = FindSection(fp, "annealing_input");            /* find tune section */
  if( !fp )
    error("ReadAParameters: could not locate annealing_parameters");
  
  fscanf(fp,"%*s\n");                         /* advance past title line 1 */
  
  if ( 1 != fscanf(fp, "%li\n", &(l_aparms.seed)) )
    error("ReadAParameters: error reading random number generator seed");

  fscanf(fp,"%*s\n");                         /* advance past title line 2 */

  if ( 1 != fscanf(fp, "%lf\n", &(l_aparms.start_tempr)) )
    error("ReadAParameters: error reading start temperature");

  fscanf(fp,"%*s\n");                         /* advance past title line 3 */

  if ( 1 != fscanf(fp, "%lf\n", &(l_aparms.gain_div_interval)) )
    error("ReadAParameters: error reading gain");

  fscanf(fp,"%*s\n");                         /* advance past title line 4 */

  /* Minor code change: interval is going to have to 
   * be a multiple of 100 now. Sorry. */
  if ( 1 != fscanf(fp, "%d\n", &(l_aparms.interval)) )
    error("ReadAParameters: error reading interval");
  /* precompute gain/interval - Seb RV 2025 Nov 14*/
  l_aparms.gain_div_interval/=(double)l_aparms.interval;
  return l_aparms;
}




/*** WriteTimes: writes the time-structure to a .times file ****************
 ***************************************************************************/

void WriteTimes(double *times)
{
  char   *timefile;                             /* name of the .times file */
  FILE   *timeptr;                         /* file pointer for .times file */
  
           /* create time file name by appending .times to input file name */
  timefile = (char *)calloc(MAX_RECORD, sizeof(char));
  timefile = strcpy(timefile, outname);
  timefile = strcat(timefile, ".times");

  timeptr = fopen(timefile, "w");
  if ( !timeptr ) {
    perror("main");
    exit(1);
  } 

  PrintTimes(timeptr, times);                /* write times to .times file */

  fclose(timeptr);                                             /* clean up */
  free(timefile);
}



/*** PrintTimes: writes two (parallel: three) times sections ***************
 ***************************************************************************/

void PrintTimes(FILE *fp, double *times)
{
  fprintf(fp, "wallclock: %.3f\n", times[0]);
  fprintf(fp, "user:      %.3f\n", times[1]);
}



/*** PrintEquil: writes an 'equilibrate_variance' section with 'title' *****
 *               to the stream specified by fp                             *
 ***************************************************************************/

void PrintEquil(FILE *fp, double *equil_var, char *title)
{
  fprintf(fp, "$%s\n", title);
  fprintf(fp, "variance:\n");
  fprintf(fp, "%g\n", equil_var[0]);
  fprintf(fp, "equilibrate_final_energy:\n");
  fprintf(fp, "%g\n", equil_var[1]);
  fprintf(fp, "$$\n");
}




