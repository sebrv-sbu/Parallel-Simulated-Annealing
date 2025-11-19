/*****************************************************************
 *                                                               *
 *  savestate.c                                                  *
 *                                                               *
 *****************************************************************
 *                                                               *
 *  Written by Seb RV                                            *
 *                                                               *
 *****************************************************************
 *                                                               *
 * savestate.c contains three functions that read, write, and    *
 * remove the state file for an annealing run. The frequency     *
 * with which states are saved can be chosen by command line op- *
 * tion -b (for backup stepsize). The state file is very useful  *
 * for the case when long annealing runs have to be interrupted  *
 * or crash for some reason or another. The run can then be res- *
 * umed by indicating the state file as an additional argument   *
 * to tsp_sa. It will find them by default too.                  *
 *                                                               *
 *****************************************************************
 *                                                               *
 * Copyright (C) ??? Sebastian Ramirez Villarreal                *
 * the full GPL copyright notice can be found in lsa.c           *
 *                                                               *
 *****************************************************************/


#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>                                        /* getopt stuff */

#include <error.h>
#include "move.h"
#include <random.h>
#include <sa.h>

#ifdef MPI
#include <MPI.h>                                               /* for myid */
#endif

/*** A STATIC VARIABLE *******************************************************/

static char *filename; /* name of state file */

/*** FUNCTIION DEFINITIONS ***************************************************/

/*** StateRead: reads Lam statistics, move state, and erand state from a     *
 *              state file and restores the annealer's state to the same     *
 *              state it was in before it got interrupted.                   *
 *     CAUTION: InitMoves must be called before calling StateRead!           *
 *****************************************************************************/
 /* Not functional right now - SRV 2025-11-16 */

void StateRead(char *statefile, Opts *options, MoveState *moveptr, double *stats,
  unsigned short *rand, double *delta){
//  FILE *infile;
//  int i;
//
//  filename = (char *)calloc(MAX_RECORD, sizeof(char));
//  filename = strcpy(filename, statefile);
//
//  /* open the state file and read it */
//  infile = fopen(statefile, "r");
//  if( !infile )
//    file_error("StateRead");
//
//  options->argv = fgets(options->argv, MAX_RECORD, infile);
//
//  fscanf(infile, "%s\n",    options->inname);
//  fscanf(infile, "%s\n",    options->outname);
//  fscanf(infile, "%d\n",  &(options->stop_flag));
//  fscanf(infile, "%d\n",  &(options->prolix_flag));
//  fscanf(infile, "%d\n",  &(options->landscape_flag));
//  fscanf(infile, "%d\n",  &(options->log_flag));
//  fscanf(infile, "%d\n",  &(options->time_flag));
//  fscanf(infile, "%ld\n", &(options->state_write));
//  fscanf(infile, "%ld\n", &(options->print_freq));
//  fscanf(infile, "%ld\n", &(options->captions));
//  fscanf(infile, "%d\n",  &(options->precision));
//  fscanf(infile, "%d\n",  &(options->quenchit));
//  fscanf(infile, "%d\n",  &(options->equil));
//#ifdef MPI
//  fscanf(infile, "%d\n",  &(options->tuning));
//  fscanf(infile, "%d\n",  &(options->covar_index));
//  fscanf(infile, "%d\n",  &(options->write_tune_stat));
//  fscanf(infile, "%d\n",  &(options->auto_stop_tune));
//#endif
//
//  if ( options->time_flag ) {
//    fscanf(infile, "%lf\n", &(delta[0]));
//    fscanf(infile, "%lf\n", &(delta[1]));
//  }
//
//  fscanf(infile, "%d\n",  &(moveptr->ncities));
//  fscanf(infile, "%d\n",  &(moveptr->nhits));
//  fscanf(infile, "%d\n",  &(moveptr->nsweeps));
//
//  moveptr->acc_tab_ptr =
//    (AccStats *)malloc(sizeof(AccStats));
//  moveptr->curr_tour=(unsigned short*)calloc(moveptr->ncities, sizeof(unsigned short));
//  moveptr->curr_position=(unsigned short *)calloc(moveptr->ncities, sizeof(unsigned short));
//  for(i=0; i < moveptr->ncities-1; i++)
//    fscanf(infile, "%d", &(moveptr->curr_tour[i]));
//  fscanf(infile,"%d \n", &moveptr->curr_tour[moveptr->ncities-1]);
//  for (i=0; i < moveptr->ncities-1; i++)
//    fscanf(infile, "%d", &(moveptr->curr_position[i]));
//  fscanf(infile,"%d \n", &(moveptr->curr_tour[moveptr->ncities-1]));
//  fscanf(infile,"%lf\n", &(moveptr->curr_cost));
//
//  fscanf( infile, "%lg %d %d\n",
//   &(moveptr->acc_tab_ptr->theta_bar),
//   &(moveptr->acc_tab_ptr->hits),
//   &(moveptr->acc_tab_ptr->success) );
//
//  for(i=0; i<31; i++)
//    fscanf(infile, "%lg\n", &(stats[i]));
//
//  for(i=0; i<3; i++)
//    fscanf(infile, "%hd\n", &(rand[i]));
//
//
//
}

/*** StateWrite: collects command line options and arguemnts, Lam statis- **
 *               tics, move state and the state of the erand48 random num- *
 *               ber generator and writes all that into the state file,    *
 *               which can then be used to restore the run in case it gets *
 *               interrupted                                               *
 ***************************************************************************/
/* Currently Not Working SRV 2025-11-16 */

void StateWrite(char *statefile){
//  int i;
//  FILE *outfile;
//  Opts *options;
//  MoveState *move_status;
//  double *lamsave;
//  unsigned short *prand;
//  double *delta;
//
//  /* if StateWrite() called for the first time, make filename static. */
//  if (filename==NULL){
//    filename = (char *)calloc(MAX_RECORD, sizeof(char));
//    filename = strcpy(filename, statefile);
//  }
//  
//  options = GetOptions();
//  move_status = MoveSave();
//  lamsave = GetLamstats();
//  prand = GetERandState();
//  if (time_flag)
//    delta = GetTimes();
//  outfile = fopen(filename, "w");
//  if ( !outfile )
//    file_error("StateWrite");
//
//  fprintf(outfile, "%s\n",      options->argv);
//
//  fprintf(outfile, "%s\n",    options->inname);
//  fprintf(outfile, "%s\n",    options->outname);
//  fprintf(outfile, "%d\n",    options->stop_flag);
//  fprintf(outfile, "%d\n",    options->prolix_flag);
//  fprintf(outfile, "%d\n",    options->landscape_flag);
//  fprintf(outfile, "%d\n",    options->log_flag);
//  fprintf(outfile, "%d\n",    options->time_flag);
//  fprintf(outfile, "%ld\n",   options->state_write);
//  fprintf(outfile, "%ld\n",   options->print_freq);
//  fprintf(outfile, "%ld\n",   options->captions);
//  fprintf(outfile, "%d\n",    options->precision);
//  fprintf(outfile, "%d\n",    options->quenchit);
//  fprintf(outfile, "%d\n",    options->equil);
//#ifdef MPI
//  fprintf(outfile, "%d\n",    options->tuning);
//  fprintf(outfile, "%d\n",    options->covar_index);
//  fprintf(outfile, "%d\n",    options->write_tune_stat);
//  fprintf(outfile, "%d\n",    options->auto_stop_tune);
//#endif
// if ( time_flag ) {
//    fprintf(outfile, "%.3f\n", delta[0]);
//    fprintf(outfile, "%.3f\n", delta[1]);
//  }
//
//  fprintf(outfile, "%d\n",    move_status->ncities);
//  fprintf(outfile, "%d\n",    move_status->nhits);
//  fprintf(outfile, "%d\n",    move_status->nsweeps);
//
//  for(i=0; i < move_status->ncities; i++)
//    fprintf(outfile, "%d ", move_status->curr_tour[i]);
//  fprintf(outfile,"\n");
//  
//  for(i=0; i < move_status->ncities; i++)
//    fprintf(outfile, "%d ", move_status->curr_position[i]);
//  fprintf(outfile, "\n");
//  fprintf(outfile, "%lf\n", move_status->curr_cost);
//
//  fprintf(outfile, "%.16g %d %d\n",
//    move_status->acc_tab_ptr->theta_bar,
//    move_status->acc_tab_ptr->hits,
//    move_status->acc_tab_ptr->success);
//
//  for(i=0; i < 31; i++)
//    fprintf(outfile,"%.16g\n", lamsave[i]);
//
//  for(i=0; i < 3; i++)
//    fprintf(outfile,"%d\n", prand[i]);
//  fclose(outfile);
//
//  free(options);
//  FreeMoveState(move_status);
//  free(lamsave);
//  if ( time_flag )
//    free(delta);
}
//
void StateRm(void){
//  if ( remove(filename) )
//    warning("StateRm: could not delete %s", filename);
//  free(filename);
}
