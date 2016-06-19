#include "streaming.h"
#include "helper.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"

void doStremingCell(double * collideField, double * streamField, int * flagField, double * massField, double * fractionField, int * node, int * n, int isFluid) {
    int i;
    int source_node[3];
    double fi_nb, se;

    for (i = 0; i < Q; i++) {
        /* neighboring cell from which particles are obtained */
        source_node[0] = node[0] - LATTICEVELOCITIES[i][0];
        source_node[1] = node[1] - LATTICEVELOCITIES[i][1];
        source_node[2] = node[2] - LATTICEVELOCITIES[i][2];

        /* no interaction between GAS cells */
        if (*getFlag(flagField, source_node, n) == GAS) {
            continue;
        }

        if (isFluid) {
            /* If it is fluid then mass is equal to density */
            double density;
            computeDensity(streamField, &density);
            *getMass(massField, node, n) = density;
        } else {
        /* Amount of particles that goes to this cell */
        fi_nb = *getEl(collideField, source_node, i, n);
        *getEl(streamField, node, i, n) = fi_nb;

        /*
          Update mass field according formula (4.3) (PhD thesis):
          dm_i(x, t + dt) = se * (E(x + dt*e_i, t) + E(x, t)) / 2
          se = f_i_inv(x + dt*e_i, t) - f_i(x,t)

          m(x, t + dt) = m(x, t) + sum(dm_i(x, t + dt));
          i in formula corresponds to Q - 1 - i in code and vice versa
        */
        se = fi_nb - *getEl(collideField, source_node, Q - 1 - i, n);
        *getMass(massField, node, n) += se * (*getFraction(fractionField, node, n) + *getFraction(fractionField, source_node, n)) / 2;
        }
    }

}

void doStreaming(double * collideField, double * streamField, int * flagField, double * massField, double * fractionField, int * length){
    int x, y, z, flag, isFluid, isInterface;
    int node[3];

    int n[3] = { length[0] + 2, length[1] + 2, length[2] + 2 };

    /* Loop for inner cells */
    for (z = 1; z <= length[2]; z++) {
        node[2] = z;
        for (y = 1; y <= length[1]; y++) {
            node[1] = y;
            for (x = 1; x <= length[0]; x++) {
                node[0] = x;
                flag = *getFlag(flagField, node, n);
                isFluid = flag == FLUID;
                isInterface = flag == INTERFACE;

                if (isFluid || isInterface) {
                    doStremingCell(collideField, streamField, flagField, massField, fractionField, node, n, isFluid);
                }
            }
        }
    }
}
