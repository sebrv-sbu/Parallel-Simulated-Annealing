#include <MPI.h>
#include <mpi.h>

void copy_dyn_int_array_to_const_int_array(int dest[], int *src, int size);
void CpyAllRanksToStack(int dest[][MAX_P], int **src, size_t size);
char** ExtractSetFromCharList(char* list, int length_list, int *length_set);
int compFun(const void *a, const void *b);
