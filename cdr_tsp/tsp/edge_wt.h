typedef struct{
  double coord_x;
  double coord_y;
} coord;

typedef struct distance_from_city_struct
{
  double distance;
  unsigned short city;
}
dist_and_city;


unsigned short ncities;
coord *node_coords;

dist_and_city *distances;

/*Static Variables Needed in move.c and others */
int GetJthNearestNeighbour(unsigned short home, unsigned short j);

void FloydRivestSelect(unsigned short left, unsigned short right, unsigned
short k);

double GetDistance(coord A, coord B);

void FreeDistances();

__attribute__((always_inline))
void swap_inplace(unsigned short left, unsigned short right);
