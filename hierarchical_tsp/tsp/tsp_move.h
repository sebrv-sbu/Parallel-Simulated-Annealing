/* my name is tsp_move.h and I was created on 4-01 by Lorraine Greenwald
* There are tsp function prototypes from tsp_move.c here */ 

/***************************************************
* allocates dynamic memory for the tsp tour arrays *
***************************************************/
int tour_allocate(void);

/*****************************************************
* deallocates dymanic memory for the tsp tour arrays * 
******************************************************/
void tour_deallocate(void);

/***************************************
* this module prints the current tour  *
****************************************/
void print_curr_tour(void);

/*******************************************
* generates starting tour and sets min tour *
********************************************/
double start_tour( void );


