int *sorted_lists;

typedef struct{
  double coord_x;
  double coord_y;
} coord;

typedef struct distance_from_city_struct
{
  double distance;
  int city;
}
dist_and_city;


int ncities;
coord *node_coords;

dist_and_city *distances;
dist_and_city *temp;

/*Static Variables Needed in move.c and others */
int GetJthNearestNeighbour(int home, int j);

double GetDistance(coord A, coord B);

void MergeSort(int home, int* final_result);

void MergeInc(int left, int increment);

void MergeLim(int left, int increment);

void FreeSortedLists();

void GenSortedLists();
