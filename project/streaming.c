#include "streaming.h"
#include "helper.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"
#include "utils.h"

/* doStremingCell: performs the streaming operation for one cell, in fact each cell receives all the streaming from the neighbor cells */
void doStremingCell(Fields &fields, int * node, std::vector<float>::iterator el, int * n, int isInterface, int isFluid, float exchange) {
    int i, flag;
    int source_node[3];
    float fi_nb, se;

    for (i = 0; i < Q; i++) {
        /* neighboring cell from which particles are obtained */
        source_node[0] = node[0] - LATTICEVELOCITIES[i][0];
        source_node[1] = node[1] - LATTICEVELOCITIES[i][1];
        source_node[2] = node[2] - LATTICEVELOCITIES[i][2];

        /* Amount of particles that goes to this cell */
        fi_nb = *getEl(fields.collide, source_node, i, n);
        if (isFluid) {
            *(el + i) = fi_nb;
        }
        if (isInterface) {
            /* Obtain the flag of the neighbor cell */
            flag = *getFlag (fields, source_node, n); 
            if (flag == GAS) {
                float velocity[3], feq[Q], rho_ref = 1;
                /* get pointer to the fluid cell */
                auto fluidCell = getEl (fields.collide, node, 0, n);
                /* compute velocity of the fluid cell */
                computeVelocity(fluidCell, &rho_ref, velocity);
                /* compute f-equilibrium of the fluid cell */
                computeFeq(&rho_ref, velocity, feq);
                /* set boundary */
                *(el + i) = feq[Q - i -1] + feq[i] - fluidCell[Q - 1 - i];
            } else {
                *(el + i) = fi_nb;
            }
            /* If the neighbor cell is fluid or interface then update the mass value */
            if (flag == FLUID || flag == INTERFACE) {
                se = fi_nb - *getEl(fields.collide, source_node, Q - 1 - i, n);
                *getMass(fields, node, n) += exchange * se * (*getFraction(fields, node, n) + *getFraction(fields, source_node, n)) * 0.5;
            }
        }
    }
}

void doStreaming(Fields &fields, int * length, int n_threads, float exchange){
    int x, y, z, flag, isFluid, isInterface;
    int node[3];
    int n[3] = { length[0] + 2, length[1] + 2, length[2] + 2 };

    /* Loop for inner cells */
    #pragma omp parallel for schedule(dynamic) private(x, node, isFluid, flag, isInterface, el) num_threads(n_threads) collapse(2)
    for (z = 1; z <= length[2]; z++) {
        for (y = 1; y <= length[1]; y++) {
            node[2] = z;
            node[1] = y;
            for (x = 1; x <= length[0]; x++) {
                /* Obtain the pointer and the flag for each element */
                node[0] = x;
                auto el = getEl (fields.stream, node, 0, n);
                flag = *getFlag (fields, node, n);
                isFluid = flag == FLUID;
                isInterface = flag == INTERFACE;
                /* Make streaming just for fluid and interface cells */
                if (isFluid || isInterface) {
                    doStremingCell (fields, node, el, n, isInterface, isFluid, exchange);
                }
            }
        }
    }
}
