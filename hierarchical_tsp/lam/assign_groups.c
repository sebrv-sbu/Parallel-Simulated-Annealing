#include <MPI.h>
#include <assign_groups.h>
#include <float.h>                                    /* for double limits */
#include <limits.h>                                  /* for integer limits */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <sys/types.h>                        /* these two are for times() */
#include <time.h>                                    /* this is for time() */
#include <unistd.h>          /* for command line option stuff and access() */
#include <stddef.h>
#include <stdint.h>

/* NOTE: do not ever dare to include moves.h in here or in any of the hea- */
/*       ders below; lsa.c must remain truly problem independent           */

#include <sa.h>
#include <error.h>

#include <MPI.h>
#include <mpi.h>
#define PACK_BYTES(buf, arr, len) \
do {  \
  memcpy((buf),(arr), (len)); \
  (buf)+= (len); \
}while (0)
#define PACK_INTS(buf,arr,len)\
do{ \
  memcpy((buf),(arr),(len)*sizeof(int));\
  (buf)+=(len);\
}while (0)

#define GROUP_IS_FULL(g) ((g)[lam_group_size - 1] != 0)



void AssignGroups(){
  char machine[MPI_MAX_PROCESSOR_NAME];                  /* For machine name */
  int name_len;
  int  i;
  char *machines;
  char **machines_set;
  char *curr_string;
  int *displacements;
  int *sizes;
  int len_and_id[2];
  int **all_ranks;
  int local_ranks[MAX_M][MAX_P];
  int root_ranks[MAX_P];
  int unique_machines=0;
  int *group_size, *group_index;
  int max_group_id;
  int curr_rank;
  int max_size = 0;
  MPI_Group world;
  local_comms =(MPI_Comm*)malloc(ngroups*sizeof(MPI_Comm));
  local_groups=(MPI_Group*)malloc(ngroups*sizeof(MPI_Group));
  displacements = calloc(nnodes,sizeof(*displacements));
  sizes         = calloc(nnodes,sizeof(*sizes));
  root_ids      = calloc(ngroups,sizeof(*root_ids));
  MPI_Get_processor_name(machine, &name_len);

  MPI_Allgather(&name_len, 1, MPI_INT,sizes, 1, MPI_INT, MPI_COMM_WORLD);

  for (i=0; i<nnodes; i++){
    sizes[i]++;
    if (sizes[i]>max_size)
      max_size=sizes[i];
  }

  for ( i = 1; i < nnodes; i++)
    displacements[i]=sizes[i-1]+displacements[i-1];
   /* Displacements tells Gatherv how far away from the pointer to the start of the buffer to put the ith entry*/
 machines = (char*)calloc(displacements[nnodes-1]+sizes[nnodes-1], sizeof(*machines));
 MPI_Allgatherv(machine, name_len, MPI_CHAR, machines, sizes, displacements, MPI_CHAR, MPI_COMM_WORLD);

   machines_set=ExtractSetFromCharList(machines, nnodes, &unique_machines);

   /*need to keep track of group ids and number of nodes in group */
   group_index=malloc(unique_machines*sizeof(*group_index));
   group_size=calloc(unique_machines,sizeof(*group_size));

   for (int i=0; i<unique_machines; ++i){ group_index[i]=i; }
   /*this allows a machine to contain multiple groups */

   /*max_group_id tells us which group id to use next if a machine has more nodes than fit
    * on a group */

    if(nnodes % ngroups == 0)
      lam_group_size=nnodes/ngroups;
    else
      perror("tsp_sa: number of processors not divisible by number of groups");

   max_group_id=unique_machines;
   curr_string=(char*)malloc(sizeof(char)*max_size);
   all_ranks=calloc((ngroups+unique_machines),sizeof(*all_ranks));
   for (i=0; i<ngroups+unique_machines; ++i)
     all_ranks[i]=calloc(lam_group_size, sizeof(**all_ranks));

   /* We now need complicated logic to assign everyone a *
    * group. Fuck.                                       *
    * We are going to implement the following algorithm. *
    *                                                    *
    * Every time that we call this function, we can make *
    * initially the lists as before, except keep one     *
    * additional list for each machine:                  *
    * [P_0 P_1 ... P_n-1][P_n...P_2n-1]...[P_kn P_kn+1 .. P_kn+x 0 0 0 0]
    * We then merge these buffers and add them to a group. */
   int machine_idx;
   for (curr_rank=0; curr_rank<nnodes; ++curr_rank){
     machine_idx=0;
     memcpy(curr_string,machines+displacements[curr_rank],sizes[curr_rank]);
     curr_string[sizes[curr_rank]-1]='\0';
     while(strcmp(curr_string,machines_set[machine_idx])!=0){
       ++machine_idx;
     }
    if (group_size[machine_idx]==lam_group_size){
       group_index[machine_idx]=max_group_id++;
       group_size[machine_idx]=0;
     }
     all_ranks[group_index[machine_idx]][group_size[machine_idx]++]=curr_rank;
   }
     int *merge_buffer=malloc(unique_machines*lam_group_size*sizeof(*merge_buffer));
     int *merge_buffer_pos = merge_buffer;
     int all_ranks_index;
     int mismatches=0;
   /* Now we merge overflows */
   for (machine_idx = 0 ; machine_idx < unique_machines; ++ machine_idx){
     if (lam_group_size == group_size[machine_idx])
       continue;
     all_ranks_index = group_size[machine_idx];
     PACK_INTS(merge_buffer_pos,
     all_ranks[group_index[machine_idx]],
     (all_ranks_index));
     mismatches += (all_ranks_index);
     }
     /* Put them back into all_ranks */
     merge_buffer_pos=merge_buffer;
     machine_idx=0;
     while (mismatches){
       if (lam_group_size == group_size[machine_idx]){
         machine_idx++;
         continue;
       }
       memcpy(all_ranks[group_index[machine_idx]],
       merge_buffer_pos, lam_group_size*sizeof(*merge_buffer));
       group_size[machine_idx]=lam_group_size;
       merge_buffer_pos+=lam_group_size;
       mismatches -= lam_group_size;
       machine_idx++;
     }
     free(merge_buffer);
     /* root ids: */
   machine_idx=0;
   int roots_assigned=0;
   for (int g_index=0; roots_assigned < ngroups; g_index++){
     if (!GROUP_IS_FULL(all_ranks[g_index]))
       continue;
     root_ids[roots_assigned++] = all_ranks[g_index][0];
     }

   
   MPI_Comm_group(MPI_COMM_WORLD,&world);
   /*now we assign all the groups. This requires a workaround to get past
    * the fact that MPI_Comm stuff needs const int[ranks] */
   CpyAllRanksToStack(local_ranks, all_ranks, ngroups);
  for (int i = 0; i < ngroups; ++i){
    MPI_Group_incl(world, lam_group_size, local_ranks[i], &local_groups[i]);
    MPI_Comm_create_group(MPI_COMM_WORLD, local_groups[i], i, &local_comms[i]);
  if (MPI_COMM_NULL != local_comms[i]) {
    my_comm  = &local_comms[i];
    my_group = &local_groups[i];
    MPI_Comm_rank(local_comms[i], &my_group_id);
    }
  }
  copy_dyn_int_array_to_const_int_array(root_ranks, root_ids, ngroups);
  MPI_Group_incl(world, ngroups, root_ranks, &root_group);
  MPI_Comm_create_group(MPI_COMM_WORLD, root_group, i, &root_comm);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  MPI_Comm_size(MPI_COMM_WORLD,&nnodes);
  /* time to free stuff. */
  free(curr_string);
  free(displacements);
  free(sizes);
  free(machines);
  free(group_index);
  free(group_size);

  for (i=0; i<ngroups; ++i){
    free(all_ranks[i]);
  }
  free(all_ranks);

  for (i=0; i<unique_machines; ++i){
    free(machines_set[i]);
  }
  free(machines_set);

 }

void copy_dyn_int_array_to_const_int_array(int dest[], int *src, int size){
  int i;
  for (i=0; i<size; ++i)
    dest[i]=src[i];
}

void CpyAllRanksToStack(int dest[][MAX_P], int **src, size_t size){
  int groups_assigned=0;
  for (int all_ranks_idx=0; groups_assigned < size ; all_ranks_idx++){
    if (!GROUP_IS_FULL(src[all_ranks_idx]))
        continue;
       memcpy(dest[groups_assigned++], src[all_ranks_idx], lam_group_size*sizeof(**src));
    }
  }
char** ExtractSetFromCharList(char* list, int length_list, int *length_set){
  char **set;
  int *pos=calloc(length_list, sizeof(*pos));
  *length_set=1;
  char *list_ptr=list;
  uint8_t comp_flag=0;

  for (int i = 1; i < length_list; i++){
    while (*list_ptr != '\0'){ list_ptr++; }
    list_ptr++;
    for (int j = 0; j<*length_set; ++j){
    /* Recall strcmp returns 0 if strings are identical (for some reason) */
      if (0 == strcmp(list_ptr, list + pos[j])) {
        comp_flag=1;
        break;
      }
    }
    if (!comp_flag){
      pos[*length_set] = list_ptr-list;
      (*length_set)++;
    }
    comp_flag=0;

    }

  set=calloc(*length_set, sizeof(*set));
  int len;
  for ( int j =0 ; j < *length_set; ++j){
    len=strlen(list+pos[j]);
    set[j]=calloc(len+1, sizeof(**set));
    strcpy(set[j], list+pos[j]);
  }
  free(pos);
  return(set);
}

int compFun(const void *a, const void *b){
    int ia=*(const int *)a;
    int ib=*(const int *)b;
    return ( (ia > ib) - (ia<ib) );
  }
