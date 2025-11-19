#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>                                          /* for getopt */

#include "move.h"
#include "edge_wt.h"
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
  
  unsigned short *curr_tour;
  extern char *optarg;
  extern int optind;
  extern int optopt;

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
  ReadTSP(instance_file);
  output_file = fopen(output_name, "r");
  output_file=FindSection(output_file, "final_state");
  fgets(buff, MAX_RECORD, output_file);
  sscanf(buff, "%*s\n");
  if (debug_flag){
    fgets(buff, MAX_RECORD, output_file);
    sscanf(buff, "annealing minimum cost is: %lf", &orig_cost);
  }
  fgets(buff, MAX_RECORD, output_file);
  curr_tour=(unsigned short*)calloc(ncities, sizeof(unsigned short));
  for(i=0; i<ncities; i++){
    fscanf(output_file, "%d", &curr_tour[i]); 
  }
  cost=tour_cost(curr_tour);
  format=(char*)calloc(MAX_RECORD, sizeof(char));
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
  free(curr_tour); 
  return 0;
}
