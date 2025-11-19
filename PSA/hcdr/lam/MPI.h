/*****************************************************************
 *                                                               *
 *   MPI.h                                                       *
 *                                                               *
 *****************************************************************
 *                                                               *
 *   written by John Reinitz                                     *
 *   modified by King-Wai Chu and Johannes Jaeger                *
 *                                                               *
 *   last modified, 06/13/2002                                   *
 *                                                               *
 *****************************************************************
 *                                                               *
 * IMPORTANT: IF YOU EVER CHANGE ANYTHING IN THIS FILE, LET ALL  *
 *            YOUR FELLOW PROGRAMMERS KNOW WELL IN ADVANCE AND   *
 *            CONSULT WITH THEM IF THEY AGREE ON YOUR CHANGES!!  *
 *                                                               *
 *****************************************************************
 *                                                               *
 * MPI.h contains structs and constants that are specific to     *
 * parallel annealing code using MPI.                            *
 *                                                               *
 * This includes prototypes for all the tuning functions below.  *
 *                                                               *
 * There are two problem-specific functions declared below that  *
 * need to be defined in move(s).c.                              *
 *                                                               *
 *****************************************************************
 *                                                               *
 * NOTE: this header only contains prototypes for functions used *
 *       for parallel annealing code only; all prototypes of     *
 *       functions that include serial code need to go into sa.h *
 *                                                               *
 *****************************************************************/

#ifndef MPI_INCLUDED
#define MPI_INCLUDED
#include <mpi.h>


/*** CONSTANTS *************************************************************/

#define MAX_MIX         10000   /* max number of mixes during a tuning run */
                    /* set this to a lower number if you run out of memory */
#define GROUP_SIZE         10    /* group size for calculating upper bound */

#define STOP_TUNE_CNT    10000000                       /* stop tune count */
#define STOP_TUNE_CRIT   0.05                     /* tuning stop criterion */

#define LSTAT_LENGTH        1    /* length of Lam msg array when annealing */
#define LSTAT_LENGTH_TUNE  28       /* length of Lam msg array when tuning */
#define GSTAT_LENGTH       20 /* length of global Lam msg array when       *
                                 annealing                                 */

#define MAX_P 200

/*** PARALLEL GLOBALS ******************************************************/

int myid;                                  /* id of local node (processor) */
int nnodes;               /* number of nodes (processors), also known as P */
int tuning;                           /* flag for switching on tuning mode */
int covar_index;    /* covariance sample index for tuning (in 'tau' units) */
int write_tune_stat;               /* how often to write tuning statistics */
int auto_stop_tune;       /* auto stop tune flag to stop tuning runs early */
int write_llog;                        /* flag for writing local log files */ 
int logging_mix;       /* flag for writing detailed logs on mixing process */
int lam_group_size;  /* The size of an MPI group. In future, may change to *
                      * array if we want multiple groups of different size */
int ngroups;        /* number of groups */
int *root_ids;              /* array for root ids necessary for prolix and *
                             * equilibration runs                          */

unsigned short score_method;
int glob_interval;

MPI_Group world;
MPI_Comm *my_comm;
MPI_Group *my_group;
MPI_Group *local_groups;
MPI_Comm *local_comms;
MPI_Group root_group;
MPI_Comm root_comm;

/* logging_mix added by Seb on 31 Jul 2023 for debugging purposes. For     *
 * more info, see comments under DoMix and WriteMixLog in lsa.c            */

/*** FUNCTION PROTOTYPES ***************************************************/

/* lsa.c: parallel non-tuning funcs: update func for local Lam parameters */

/*** UpdateLParameter: update local parameters l_A, l_B, l_D and l_E and ***
 *                    the local estimators for mean and standard deviation *
 *                    for both upper and lower bounds of M_Opt for the     *
 *                    current S                                            *
 ***************************************************************************/

void UpdateLParameter(void);



/* lsa.c: parallel non-tuning funcs: mixing functions */

/*** DoMix: does the mixing; sends move state and local Lam stats to the ***
 *          dance partner(s)                                               *
 *          Local mix occurs only inside local group. Global occurs        *
 *          between all communicators.                                     *
 ***************************************************************************/

void DoMix(void);

void DoGlobalMix(void);

void DoLocalMix(void);

void AssignDancePartner(int nodesInMix, MPI_Comm comm, double score);
/*** DoFixMix: for equilibration run, we only need to pass the move state **
 *             since Lam stats are not needed at constant temperature      *
 ***************************************************************************/

void DoFixMix(void);

/*** MakeLamMsg: packages local Lam stats into send buffer *****************
 ***************************************************************************/

void MakeLamMsg(unsigned char**sendbuf, MPI_Aint size);

/*** AcceptLamMsg: receives new energy and Lam stats upon mixing ***********
 ***************************************************************************/

void AcceptLamMsg(unsigned char**recvbuf, MPI_Aint size);



/* lsa.c: tuning functions */

/*** InitTuning: sets up/restores structs and variables for tuning runs ****
 ***************************************************************************/

void InitTuning(void);

/*** DoTuning: calculates the cross-correlation (for lower bound) and the **
 *             variance of local means (for upper bound) for a sub_tune_   *
 *             interval; the results are added to arrays, which span the   *
 *             whole tuning_interval; when these values get written to the *
 *             'bound' files, they will be divided by the number of mixes  *
 *             we have done already to average them out (see also King-Wai *
 *             Chu's thesis, p. 63)                                        *
 ***************************************************************************/

void DoTuning(void);

/*** WriteTuning: writes tuning stats to files ever tune_interval **********
 ***************************************************************************/

void WriteTuning(void);

/*** StopTuning: is to tuning runs what Frozen() is to a normal annealing **
 *               run: it basically checks if the tuning stop criterion     *
 *               applies and returns true if that's the case               *
 ***************************************************************************/

int StopTuning(void);




/* move(s).c: functions for communicating move state for mixing */

/*** MakeStateMsg: function to prepare a message which is then passed ******
 *                 to other nodes via MPI telling the other nodes about    *
 *                 move stats; since we don't know about the move state    *
 *                 structs, but can assume that we'll only have to send    *
 *                 longs and doubles, we split the message in two arrays,  *
 *                 one for the longs and one for the doubles; then we re-  *
 *                 turn the arrays and their sizes for lsa.c to send them  *
 *                 to the dance partners                                   *
 *                                                                         *
 *                 note that all the arguments need to get passed by refe- *
 *                 rence since we need to allocate the arrays according to *
 *                 their problem-specific size in move(s).c                *
 ***************************************************************************/

void MakeStateMsg(unsigned char **buf, MPI_Aint padding, MPI_Aint *size);

void MakeGlobalLamMsg(unsigned char **sendbuf, MPI_Aint size);
/*** AcceptMsg: communicates a message about move stats received via MPI ***
 *              to move(s).c; see the comment for MakeStateMsg for the ra- *
 *              tionale behind the two arrays that are passed              *
 ***************************************************************************/

void AcceptStateMsg(unsigned char **buf);

void AcceptGlobalLamMsg(unsigned char **recvbuf, MPI_Aint size);


/*** WriteMixLog: writes a detailed log giving information about each      *
 *                nodes energy, probability of being chosen, and dance     *
 *                partner at each mix interval.                            *
 *                                                                         *
 * NOTE: As of now, this will not be restored properly if a run is         *
 * restarted. This feature may be implemented at some point, but not for   *
 * now.                                                                    */

void WriteMixLog(double *node_prob, int* dance_partner);

/* WriteMixLog added by Seb on 31 Jul 2023 for debugging purposes. For     *
 * more info, see comments under DoMix and WriteMixLog in lsa.c            */

/*** AssignGroups: Writes array of ranks to put in each group.           ***/

void AssignGroups();

/***Some utility functions for AssignGroups:                             ***/
char** extract_string_set_from_char_list(char* list, int length_list, int *length_set, int *sizes, int *displacements, int *set_sizes);

void shift_left(char **buff, int length_buff, int length_loop, int *set_sizes, int start_point);

void copy_dyn_int_array_to_const_int_array(int dest[], int *src, int size);

/*** UpdateGlobalStats: Updates statistics in each group.                ***/
void UpdateGlobStats(double Inv_Sum, double *new_stats);

/* Create an inverse sum operation for MPI to use */
#endif
