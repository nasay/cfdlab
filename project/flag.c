#include "LBDefinitions.h"
#include "helper.h"
#include "flag.h"
#include "computeCellValues.h"
#include <omp.h>

void makeAvgDistFn(float * collideField, int * flagField, int * n, int * cell) {
    /*
    A GAS cell that is promoted to INTERFACE needs an initial distribution function, which 
    is calculated via f_eq(rho_avg, v_avg), 
    where v_avg, rho_avg are averaged from the neighboring FLUID and INTERFACE cells.

    collideField  An array of DFs for each cell in the domain, excluding boundary cells
    n             The dimensions of the domain, including boundary cells
    cell          The coordinates of the cell in need of a DF
    */    

    // TODO: How does this conserve mass?
    // TODO: Do we have to set the mass field?
    // TODO: Can we make this faster?

    int i, neighbor[3], nNeighbors, flag;
    float density, density_avg, velocity[D], velocity_avg[D], * cellDF, * neighborDF;
    

    cellDF = getEl(collideField, cell, 0, n);

    nNeighbors = 0;
    density_avg = 0;
    velocity_avg[0] = 0;
    velocity_avg[1] = 0;
    velocity_avg[2] = 0;

    // for each i <- lattice direction
    for (i = 0; i < Q; i++) {

        // Retrieve coordinates of neighbor in direction i
        neighbor[0] = cell[0] + LATTICEVELOCITIES[i][0];
        neighbor[1] = cell[1] + LATTICEVELOCITIES[i][1];
        neighbor[2] = cell[2] + LATTICEVELOCITIES[i][2];

        // Do not overstep boundaries
        if (neighbor[0] < 1 || neighbor[0] > n[0]-2 || neighbor[1] == 0 || 
            neighbor[1] == n[1] || neighbor[2] == 0 || neighbor[2] == n[2]) {
            continue;
        }

        flag = *getFlag(flagField, neighbor, n);

        if (flag != FLUID && flag != INTERFACE) {
            continue;
        }

        // Retrieve distribution function of that neighbor
        neighborDF = getEl(collideField, neighbor, 0, n);

        // Extract density, velocity from that neighbor
        computeDensity(neighborDF, &density);
        computeVelocity(neighborDF, &density, velocity);

        nNeighbors++;
        density_avg += density;
        velocity_avg[0] += velocity[0];
        velocity_avg[1] += velocity[1];
        velocity_avg[2] += velocity[2];
    }

    density_avg /= nNeighbors;
    velocity_avg[0] /= nNeighbors;
    velocity_avg[1] /= nNeighbors;
    velocity_avg[2] /= nNeighbors;
    computeFeq(&density_avg, velocity_avg, cellDF);
}


// TODO Use a real data structure
void removeFromEmptyList(int ** emptiedCells, int * nEmptied, int * targetCell) {

    //printf("REMOVING: %d,%d,%d\n", targetCell[0], targetCell[1], targetCell[2]);
    // printf("BEFORE: List contains: ");
    // for (int k=0; k < *nEmptied; k++) {
    //     printf("(%d,%d,%d), ", emptiedCells[k][0],emptiedCells[k][1],emptiedCells[k][2]);

    // }
    // printf("\n");
    int j = 0; 
    while (j < *nEmptied && !(
        emptiedCells[j][0] == targetCell[0] && 
        emptiedCells[j][1] == targetCell[1] &&
        emptiedCells[j][2] == targetCell[2])) {
        j++;
    }
    
    emptiedCells[j][0] = -1;
    emptiedCells[j][0] = -1;
    emptiedCells[j][0] = -1;
  
    // j is either nEmptied or index of match
//    for (int k=j; k < (*nEmptied)-1; k++) {
//        emptiedCells[k][0] = emptiedCells[k+1][0];
//        emptiedCells[k][1] = emptiedCells[k+1][1];
//        emptiedCells[k][2] = emptiedCells[k+1][2];
//    }
//    // Decrement list length
//    if (j != *nEmptied) {
//        // printf("*************** REMOVED SOMETHING ******************\n");
//        (*nEmptied)--;
//    }

    // printf("After list contains: ");
    // for (int k=0; k < *nEmptied; k++) {
    //     printf("(%d,%d,%d), ", emptiedCells[k][0],emptiedCells[k][1],emptiedCells[k][2]);
    // }
    // printf("\n");
}


void performFill(float * collideField, int * flagField, int * n, int ** filledCells, int nFilled, int ** emptiedCells, int * nEmptied) {
    /*
    For collections of interface cells that get emptied or filled, examine the neighboring cells 
    and update their flags to maintain the integrity of the interface layer. For example, if a cell 
    is updated to FLUID, any neighboring GAS cells become INTERFACE.

    collideField  An array of DFs for each cell in the domain
    flagField     An array of flags <- {FLUID, INTERFACE, GAS, ...} for each cell in domain
    n             The dimensions of flagField
    filledCells   An ?x3 array containing coordinates of cells which have just been filled
    nFilled       The length of emptiedCells
    emptiedCells  An ?x3 array containing coordinates of cells which have just been emptied
    nEmptied      The length of emptiedCells
    */

    int i, k, neighbor[3], *flag;

    // for each k <- cell that has been updated
#pragma omp parallel for schedule(dynamic) private(i, neighbor, flag)
    for (k = 0; k < nFilled; k++) {         

        // Update the cell's own flag
        *getFlag(flagField, filledCells[k], n) = FLUID;

        // Update each neighbor to ensure mesh integrity
        for (i = 0; i < Q; i++) {            // for each i <- lattice direction

            // Retrieve coordinates of neighbor in direction i
            neighbor[0] = filledCells[k][0] + LATTICEVELOCITIES[i][0];
            neighbor[1] = filledCells[k][1] + LATTICEVELOCITIES[i][1];
            neighbor[2] = filledCells[k][2] + LATTICEVELOCITIES[i][2];


            // Check if neighbor is on the domain boundary, in which case we ignore it
            // TODO: See if this actually makes things faster, if not, delete it. 
            //       Unless someone sets the FLUID flag as a domain boundary condition, 
            //       which would probably cause lots of problems, this won't do anything.
            if (neighbor[0] == 0 || neighbor[0] == n[0] || neighbor[1] == 0 || 
                neighbor[1] == n[1] || neighbor[2] == 0 || neighbor[2] == n[2]) {
                continue;
            }

            // Retrieve the flag corresponding to neighbor
            flag = getFlag(flagField, neighbor, n);

            // If neighbor needs to be updated
            if (*flag == GAS) {

                *flag = INTERFACE;

                // update distribution function from average of neighbors
                makeAvgDistFn(collideField, flagField, n, neighbor);

                // Remove this neighbor from 'empty' list
                removeFromEmptyList(emptiedCells, nEmptied, neighbor);
            }
        }
    }
}

void performEmpty(float * collideField, int * flagField, int * n, int ** updatedCells, int nUpdated) {
    /*
    For collections of interface cells that get emptied or filled, examine the neighboring cells 
    and update their flags to maintain the integrity of the interface layer. For example, if a cell 
    is updated to FLUID, any neighboring GAS cells become INTERFACE.

    collideField  An array of DFs for each cell in the domain
    flagField     An array of flags <- {FLUID, INTERFACE, GAS, ...} for each cell in domain
    n             The dimensions of flagField
    updatedCells  An ?x3 array containing coordinates of cells which have just been updated
    nUpdated      The length of updatedCells
    */

    int i, k, neighbor[3], * flag;

#pragma omp parallel for schedule(dynamic) private(i, neighbor, flag)
    // for each k <- cell that has been updated
    for (k = 0; k < nUpdated; k++) {         
        if (updatedCells[k][0] == -1) continue;
        
        * getFlag(flagField, updatedCells[k], n) = GAS;

        // Heal interface by updating neighbors
        for (i = 0; i < Q; i++) {            // for each i <- lattice direction

            // Retrieve coordinates of neighbor in direction i
            neighbor[0] = updatedCells[k][0] + LATTICEVELOCITIES[i][0];
            neighbor[1] = updatedCells[k][1] + LATTICEVELOCITIES[i][1];
            neighbor[2] = updatedCells[k][2] + LATTICEVELOCITIES[i][2];

            // Check if neighbor is on the domain boundary, in which case we ignore it
            // TODO: See if this actually makes things faster, if not, delete it. 
            //       Unless someone sets the FLUID flag as a domain boundary condition, 
            //       which would probably cause lots of problems, this won't do anything.
            if (neighbor[0] == 0 || neighbor[0] == n[0] || neighbor[1] == 0 || 
                neighbor[1] == n[1] || neighbor[2] == 0 || neighbor[2] == n[2]) {
                continue;
            }

            // Retrieve the flag corresponding to neighbor
            flag = getFlag(flagField, neighbor, n);

            // Swap out flag
            if (*flag == FLUID) {
                *flag = INTERFACE;
            }
        }
    }
}

void updateFlagField(float * collideField, int * flagField, float * fractionField, int * length) {
    int x, y, z, i, flag, nFilled = 0, nEmptied = 0;
    int node[3];
    float fraction, eps = 1e-3;
    int n[3] = { length[0] + 2, length[1] + 2, length[2] + 2 };
    
    // int filledCells[n[0] * n[1] * n[2]][3];
    // int emptiedCells[n[0] * n[1] * n[2]][3];
    int **filledCells = (int **) malloc((size_t)( n[0] * n[1] * n[2] * sizeof( int * )));
    int **emptiedCells = (int **) malloc((size_t)( n[0] * n[1] * n[2] * sizeof( int * )));
    for (i = 0; i < (n[0] * n[1] * n[2]); ++i){
        filledCells[i] = (int *) malloc((size_t)( 3 * sizeof( int )));
        emptiedCells[i] = (int *) malloc((size_t)( 3 * sizeof( int )));
    }

    /*
      Updating flags for INTERFACE cells:
         if fraction > 1 set flag to FLUID;
         if fraction < 0 set flag to GAS.
      Saving all emptied and filled cells to emptiedCells and filledCells arrays.
    */
    /* TODO fill filledCells and emptiedCells on collide step */
    for (z = 1; z <= length[2]; z++) {
        node[2] = z;
        for (y = 1; y <= length[1]; y++) {
            node[1] = y;
            for (x = 1; x <= length[0]; x++) {
                node[0] = x;
                flag = *getFlag(flagField, node, n);

                /* We are interested only in INTERFACE cells now */
                if (flag == INTERFACE) {
                    fraction = *getFraction(fractionField, node, n);

                    if (fraction > 1 + eps) {
                        filledCells[nFilled][0] = node[0];
                        filledCells[nFilled][1] = node[1];
                        filledCells[nFilled][2] = node[2];
                        nFilled++;

                    } else if (fraction < -eps) {
                        emptiedCells[nEmptied][0] = node[0];
                        emptiedCells[nEmptied][1] = node[1];
                        emptiedCells[nEmptied][2] = node[2];
                        nEmptied++;
                    }
                }
            }
        }
    }


    // Update neighbors of filled and emptied cells in order to have closed interface layer
    performFill(collideField, flagField, n, filledCells, nFilled, emptiedCells, &nEmptied);
    performEmpty(collideField, flagField, n, emptiedCells, nEmptied);

    for (i = 0; i < (n[0] * n[1] * n[2]); ++i){
        free(filledCells[i]);
        free(emptiedCells[i]);
    }
    free(filledCells);
    free(emptiedCells);
    

    // TODO: Redistribute mass
}






