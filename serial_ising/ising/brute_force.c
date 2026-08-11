#include "move.h"
#include "brute_force.h"
#include "initialize.h"
#include <stdlib.h>
#include "sa.h"

double BruteForceOptimum() {
  state = (NucStateType *)malloc(sizeof(NucStateType));
    BestSolution best;
    best.cost = 1e100;
    // Try all 2^16 = 65,536 configurations
    for (unsigned int config = 0; config < (1u << n_points); ++config) {
        // Set spins according to bit pattern
        for (int i = 0; i < n_points; ++i) {
            VERTEX_SPIN(i) = (config & (1u << i)) ? 1 : -1;
        }
        
        double cost = CalcCost();
        
        if (cost < best.cost) {
            best.cost = cost;
            // Save configuration
            for (int i = 0; i < n_points; ++i) {
                best.config[i] = VERTEX_SPIN(i);
            }
        }
        
        // Progress indicator
        if (config % 10000 == 0) {
            printf("Checked %u / %u configurations...\r", config, 1u << n_points);
            fflush(stdout);
        }
    }
    
    printf("\nBrute force complete!\n");
    printf("Global minimum: %f\n", best.cost);
    printf("Optimal configuration: ");
    for (int i = 0; i < n_points; ++i) {
        printf("%+d ", best.config[i]);
    }
    printf("\n");

// Set to optimal config
  for (int i = 0; i< n_points; ++i){
    VERTEX_SPIN(i) = best.config[i];
  }

  return CalcCost()+SumWeights();
  }
