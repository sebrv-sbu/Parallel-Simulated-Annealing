#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>                                          /* for getopt */

#include "move.h"
#include "error.h"
#include "initialize.h"

/*** Constants *************************************************************/

#define  OPTS      ":Dp:f:h:v"  /* command line option string */

static const char usage[]  =

"Usage: printscore [-D] [-p <float_prec>] [-f <output_file>] [-h][-v]\n";

static const char help[]  =
"Usage: printscore [options] <datafile>\n\n"

"Arguments:\n"
"  -D                  debugging mode, prints out all kinds of debugging info\n"
"  -p <float_prec>     float precision of output is <float_prec>\n"
"  -f <output_file>        set the output file name to <output_file>\n"
"  -h                  prints this help message\n"
"  -v                  print version and compilation date\n"

"Please fix any bugs that you find.\n";

static const char verstring[] =

"%s version %s\n"
"compiled by:      %s\n"
"         on:      %s\n"
"      using:      %s\n"
"      flags:      %s\n"
"       date:      %s at %s\n";

int main(int argc, char **argv){
  int  c;
  FILE *instance_file;
  FILE *output_file;
  char *instance_name;
  char *output_name;
  char *inname;
  short debug_flag=0, diff_output_file=0, prec_flag=0;
  
  int i;
  int precision;
  char *section_title;
  char *format;
  double orig_cost;
  double cost;
  char buff[MAX_RECORD];
  int neighbour;
  double weight;
  
  extern char *optarg;
  extern int optind;
  extern int optopt;

  int n_points;
  int n_edges;
  int site_index, neighbour_index;
  int num_used;
  spin_weight *weights_and_spin_states;

  int nsought;
  char *base;
  long looksite;
  int degree;

  section_title=(char*)calloc(MAX_RECORD, sizeof(char));
  section_title=strcpy(section_title,"final_state");
  optarg = NULL;
  while((c=getopt(argc, argv, OPTS))!= -1){
    switch (c) {
      case 'D':
        debug_flag=1;
        break;
      case 'p':
        precision=atoi(optarg);
        prec_flag=1;
        break;
      case 'f':
        diff_output_file=1;
        output_name=(char *)calloc(MAX_RECORD + 1, sizeof(char));
        strcpy(output_name, optarg);
        break;
      case 'h':
        PrintMsg(help,0);
        break;
      case 'v':
        printf("WIP\n");
        exit(0);
        break;
      case ':':
        error("printscore: need an argument for option -%c", optopt);
        break;
      case '?':
      default:
        error("printscore: unrecognised option -%c",optopt);
    }
  }
  if( (argc - (optind -1)) != 2)
    PrintMsg(usage, 1);
  if (prec_flag == 0)
    precision=6;
  inname = (char*)calloc(MAX_RECORD, sizeof(char));
  strcpy(inname, argv[optind]);
  instance_name=(char *)calloc(MAX_RECORD, sizeof(char));
  sprintf(instance_name, "%s.instance", inname);
  if ( diff_output_file == 0 ){
    output_name = (char*)calloc(MAX_RECORD, sizeof(char));
    sprintf(output_name, "%s.output", inname);
  }
  instance_file=fopen(instance_name, "r");


  fscanf(instance_file, "%d %d", &n_points, &n_edges);
  /* For now we assume regular toroidal */
  weights_and_spin_states = 
    malloc(sizeof(*weights_and_spin_states)*(n_edges*2+n_points));
  degree = 2*(n_edges)/n_points;
  bool *used = calloc(sizeof(*used), n_points*degree);
  if (used == NULL){ printf("Malloc error"); exit(1); }

  char line_data[180];
  line_data[0] = 1;
  while (fgets(line_data, sizeof(line_data), instance_file) != NULL){
    if ( sscanf(line_data, "%d %d %lf", 
      &site_index, &neighbour_index, &weight) != 3
      ) continue;
    site_index --;
    neighbour_index -- ;

    for (num_used = 0; used[site_index*degree + num_used]; ++num_used);

    VERTEX_EDGE(site_index, num_used) = 
      (edge_weight){.neighbour = neighbour_index, .weight = weight};

    used[site_index*degree+num_used]=true;
    
    for (num_used = 0; used[neighbour_index*degree + num_used]; ++num_used);

    VERTEX_EDGE(neighbour_index, num_used) =
      (edge_weight){.neighbour = site_index, .weight = weight};
    used[neighbour_index*degree+num_used]=true;
  }

  output_file = fopen(output_name, "r");

  rewind(output_file);                /* start looking at the beginning of the file */
  
  nsought = strlen("final_state");
  base = (char *)calloc(MAX_RECORD, sizeof(char));
  
/*** while loop goes through the file character by character... ************/
   while ( (c=getc(output_file)) != EOF) {               /* ...until if finds a '$' */
    
    if ( c == '$') {                  /* found a sectioning control string */
      looksite = ftell(output_file);                               /* where are we? */
      base = fgets(base,MAX_RECORD,output_file);  /* get full sect title (without $)*/
      if ( !(strncmp(base, "final_state", nsought)) ) {
        fseek(output_file, looksite, 0);        /* found the sought string: reposi- */
        fscanf(output_file, "%*s\n");              /* tion the pointer to '$', then */
        break;
	    free(base);
      }                                     /* section title and return it */
      else {                          /* didn't find it: skip this control */
        fseek(output_file, looksite, 0);                    /* record, keep looking */
        fscanf(output_file, "%*s");                 /* NOTE: "%*s" advances pointer */
      }                                              /* without assignment */
    }
  }
   
  free(base);


  fgets(buff, MAX_RECORD, output_file);
  sscanf(buff, "%*s\n");
  if (debug_flag){
    fgets(buff, MAX_RECORD, output_file);
    sscanf(buff, "annealing minimum cost is: %lf", &orig_cost);
  }
  fgets(buff, MAX_RECORD, output_file);
  for(i=0; i<n_points; i++){
    fscanf(output_file, "%d", &VERTEX_SPIN(i)); 
  }
  
  for (int i = 0 ; i < n_points ; ++i){
    for (int j = 0 ; j < degree; ++j){
      neighbour = VERTEX_EDGE(i,j).neighbour;
      weight = VERTEX_EDGE(i,j).weight;
      if(neighbour > i) continue;
      cost += VERTEX_SPIN(i)*VERTEX_SPIN(neighbour) * weight;
    }
  }
  format = calloc(MAX_RECORD,sizeof(*format));
  if (debug_flag){
    sprintf(format, "Cost from program: %%.%dlf\\n", precision);
    printf(format, orig_cost);
  }
  sprintf(format, "Tour Cost: %%.%dlf\n", precision);
  printf(format, cost);
  fclose(output_file);
  fclose(instance_file);
  free(format);
  free(section_title);

  free(instance_name);
  free(inname);
  free(output_name);
  free(weights_and_spin_states);
  free(used);
  /* FIX LATER */
  return 0;
}
