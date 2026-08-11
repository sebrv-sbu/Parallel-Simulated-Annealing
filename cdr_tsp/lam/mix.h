#include <mpi.h>
/* For Faster Mean and Variance computations */
typedef struct{
  double mean;
  double vari;
  uint8_t success ;
} mean_vari_succ;

MPI_Datatype MPI_Meanvarisucc;

extern MPI_Datatype MPI_Meanvarisucc;
extern MPI_Op MPI_Meanvarisucc_sum;

void Meanvarisucc_MPI_Init();
void Meanvarisucc_MPI_Free();
