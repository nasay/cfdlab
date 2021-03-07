#include "LBDefinitions.h"
#include "helper.h"
#include "flag.h"
#include "utils.h"
#include "computeCellValues.h"

#include <vector>

void makeAvgDistFn(Fields &fields, int * n, int * cell) {
    /*
    A GAS cell that is promoted to INTERFACE needs an initial distribution function, which 
    is calculated via f_eq(rho_avg, v_avg), 
    where v_avg, rho_avg are averaged from the neighboring FLUID and INTERFACE cells.

    collideField  An array of DFs for each cell in the domain, excluding boundary cells
    n             The dimensions of the domain, including boundary cells
    cell          The coordinates of the cell in need of a DF
    */

    int i, nNeighbors, flag;
    float density, density_avg, nNeighbors_inv;
    Velocity velocity, velocity_avg;
    Cell neighbor;

    auto cellDF = getEl(fields.collide, cell, 0, n);

    nNeighbors = 0;
    density_avg = 0;
    velocity_avg.x = 0;
    velocity_avg.y = 0;
    velocity_avg.z = 0;

    // for each i <- lattice direction
    for (i = 0; i < Q; i++) {

        // Retrieve coordinates of neighbor in direction i
        neighbor.x = cell[0] + LATTICEVELOCITIES[i].x;
        neighbor.y = cell[1] + LATTICEVELOCITIES[i].y;
        neighbor.z = cell[2] + LATTICEVELOCITIES[i].z;

        // Do not overstep boundaries
        if (neighbor.x < 1 || neighbor.x > n[0]-2 || neighbor.y == 0 || 
            neighbor.y == n[1] || neighbor.z == 0 || neighbor.z == n[2]) {
            continue;
        }

        flag = *getFlag(fields, neighbor, n);

        if (flag != FLUID && flag != INTERFACE) {
            continue;
        }

        // Retrieve distribution function of that neighbor
        auto neighborDF = getEl(fields.collide, neighbor, 0, n);

        // Extract density, velocity from that neighbor
        computeDensity(neighborDF, &density);
        computeVelocity(neighborDF, &density, velocity);

        nNeighbors++;
        density_avg += density;
        velocity_avg.x += velocity.x;
        velocity_avg.y += velocity.y;
        velocity_avg.z += velocity.z;
    }

    nNeighbors_inv = 1.0 / nNeighbors;

    density_avg *= nNeighbors_inv;
    velocity_avg.x *= nNeighbors_inv;
    velocity_avg.y *= nNeighbors_inv;
    velocity_avg.z *= nNeighbors_inv;

    float feq[Q];
    /* Set cell DF as f_eq with average velocity */ 
    computeFeq(&density_avg, velocity_avg, feq);

    for (i = 0; i < Q; i++)
      *(cellDF + i) = feq[i];
}


/**
   Remove cell from emptied cell list
 */
void removeFromEmptyList(std::vector<Cell> &emptiedCells, int *nEmptied, int *targetCell) {
    int j = 0;
    while (j < *nEmptied && !(
        emptiedCells[j].x == targetCell[0] &&
        emptiedCells[j].y == targetCell[1] &&
        emptiedCells[j].z == targetCell[2])) {
        j++;
    }

    // Mark element as deleted (since coordinates should be non-negative)
    emptiedCells[j].x = -1;
    emptiedCells[j].y = -1;
    emptiedCells[j].z = -1;
}


void performFill(Fields &fields, int * n, std::vector<Cell> &filledCells, int nFilled, std::vector<Cell> &emptiedCells, int * nEmptied, int n_threads) {
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

    int i, k;
    Cell neighbor;

    // for each k <- cell that has been updated
#pragma omp parallel for schedule(dynamic) private(i, neighbor, flag) num_threads(n_threads)
    for (k = 0; k < nFilled; k++) {

        // Update the cell's own flag
        *getFlag(fields, filledCells[k], n) = FLUID;

        // Update each neighbor to ensure mesh integrity
        for (i = 0; i < Q; i++) {            // for each i <- lattice direction

            // Retrieve coordinates of neighbor in direction i
            neighbor.x = filledCells[k].x + LATTICEVELOCITIES[i].x;
            neighbor.y = filledCells[k].y + LATTICEVELOCITIES[i].y;
            neighbor.z = filledCells[k].z + LATTICEVELOCITIES[i].z;


            // Check if neighbor is on the domain boundary, in which case we ignore it
            // TODO: See if this actually makes things faster, if not, delete it. 
            //       Unless someone sets the FLUID flag as a domain boundary condition, 
            //       which would probably cause lots of problems, this won't do anything.
            if (neighbor.x == 0 || neighbor.x == n[0] || neighbor.y == 0 || 
                neighbor.y == n[1] || neighbor.z == 0 || neighbor.z == n[2]) {
                continue;
            }

            // Retrieve the flag corresponding to neighbor
            auto flag = getFlag(fields, neighbor, n);

            // If neighbor needs to be updated
            if (*flag == GAS) {

                *flag = INTERFACE;

                // update distribution function from average of neighbors
                makeAvgDistFn(fields, n, neighbor);

                // Remove this neighbor from 'empty' list
                removeFromEmptyList(emptiedCells, nEmptied, neighbor);
            }
        }
    }
}

void performEmpty(Fields &fields, int * n, std::vector<Cell> &updatedCells, int nUpdated, int n_threads) {
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

    int i, k;
    Cell neighbor;

#pragma omp parallel for schedule(dynamic) private(i, neighbor, flag) num_threads(n_threads)
    // for each k <- cell that has been updated
    for (k = 0; k < nUpdated; k++) {
        if (updatedCells[k].x == -1) continue;

        *getFlag(fields, updatedCells[k], n) = GAS;

        // Heal interface by updating neighbors
        for (i = 0; i < Q; i++) {            // for each i <- lattice direction

            // Retrieve coordinates of neighbor in direction i
            neighbor.x = updatedCells[k].x + LATTICEVELOCITIES[i].x;
            neighbor.y = updatedCells[k].y + LATTICEVELOCITIES[i].y;
            neighbor.z = updatedCells[k].z + LATTICEVELOCITIES[i].z;

            // Check if neighbor is on the domain boundary, in which case we ignore it
            // TODO: See if this actually makes things faster, if not, delete it. 
            //       Unless someone sets the FLUID flag as a domain boundary condition, 
            //       which would probably cause lots of problems, this won't do anything.
            if (neighbor.x == 0 || neighbor.x == n[0] || neighbor.y == 0 || 
                neighbor.y == n[1] || neighbor.z == 0 || neighbor.z == n[2]) {
                continue;
            }

            // Retrieve the flag corresponding to neighbor
            auto flag = getFlag(fields, neighbor, n);

            // Swap out flag
            if (*flag == FLUID) {
                *flag = INTERFACE;
            }
        }
    }
}

void updateFlagField(Fields &fields, std::vector<Cell> &filledCells, std::vector<Cell> &emptiedCells, int * length, int n_threads) {
    int x, y, z, flag, nFilled = 0, nEmptied = 0;
    Cell node;
    float fraction, eps = 1e-3;
    int n[3] = { length[0] + 2, length[1] + 2, length[2] + 2 };

    /*
      Updating flags for INTERFACE cells:
         if fraction > 1 set flag to FLUID;
         if fraction < 0 set flag to GAS.
      Saving all emptied and filled cells to emptiedCells and filledCells arrays.
    */
    for (z = 1; z <= length[2]; z++) {
        node.z = z;
        for (y = 1; y <= length[1]; y++) {
            node.y = y;
            for (x = 1; x <= length[0]; x++) {
                node.x = x;
                flag = *getFlag(fields, node, n);

                /* We are interested only in INTERFACE cells now */
                if (flag == INTERFACE) {
                    fraction = *getFraction(fields, node, n);

                    if (fraction > 1 + eps) {
                      filledCells[nFilled] = std::move(node);
                      nFilled++;

                    } else if (fraction < -eps) {
                      emptiedCells[nEmptied] = std::move(node);
                      nEmptied++;
                    }
                }
            }
        }
    }

    // Update neighbors of filled and emptied cells in order to have closed interface layer
    performFill(fields, n, filledCells, nFilled, emptiedCells, &nEmptied, n_threads);
    performEmpty(fields, n, emptiedCells, nEmptied, n_threads);
}

