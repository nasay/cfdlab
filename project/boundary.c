#include "boundary.h"
#include "LBDefinitions.h"
#include "helper.h"
#include "computeCellValues.h"

/**
 * Set NOSLIP condition
 */
void setNoSlip(float * collideField, int * flagField, int * node, int * n) {
    int i, coord_dest[3], flag;
    float * cell_ptr;

    /* for each lattice */
    for (i = 0; i < Q; i++) {
        /* compute a cell where lattice is pointing*/
        coord_dest[0] = node[0] + LATTICEVELOCITIES[i][0];
        coord_dest[1] = node[1] + LATTICEVELOCITIES[i][1];
        coord_dest[2] = node[2] + LATTICEVELOCITIES[i][2];

        /* does the pointed cell lay in our domain? */
        if (coord_dest[0] < n[0] && coord_dest[1] < n[1] && coord_dest[2] < n[2] &&
            coord_dest[0] >= 0 && coord_dest[1] >= 0 && coord_dest[2] >= 0) {
            flag = *getFlag(flagField, coord_dest, n);
            /* if pointed cell is FLUID or INTERFACE */
            if (flag == FLUID || flag == INTERFACE) {
                /* get pointer to the i-th lattice of boundary cell */
                cell_ptr = getEl(collideField, node, i, n);
    
                /* NOSLIP */
                /* set i-th lattice to inverse lattice of the computed inner cell */
                *cell_ptr= *getEl(collideField, coord_dest, Q-1-i, n);
            }
        }
    }
}

/**
 * Set MOVING_WALL condition
 */
void setMovingWall(float * collideField, int * flagField,  const float * const wallVelocity, int * node, int * n) {
    int i, coord_dest[3], flag;
    float * cell_ptr;
    float dotProd;
    float density;
    
/* for each lattice */
    for (i = 0; i < Q; i++) {
        /* compute a cell where lattice is pointing*/
        coord_dest[0] = node[0] + LATTICEVELOCITIES[i][0];
        coord_dest[1] = node[1] + LATTICEVELOCITIES[i][1];
        coord_dest[2] = node[2] + LATTICEVELOCITIES[i][2];

        /* does the pointed cell lay in our domain? */
        if (coord_dest[0] < n[0] && coord_dest[1] < n[1] && coord_dest[2] < n[2] &&
            coord_dest[0] >= 0 && coord_dest[1] >= 0 && coord_dest[2] >= 0 ) {

            flag = *getFlag(flagField, coord_dest, n);
            /* if pointed cell is FLUID */
            if (flag == FLUID || flag == INTERFACE) {
                /* get pointer to the i-th lattice of boundary cell */
                cell_ptr = getEl(collideField, node, i, n);

                /* NOSLIP */
                /* set i-th lattice to inverse lattice of the computed inner cell */
                *cell_ptr= *getEl(collideField, coord_dest, 18-i, n);

                dotProd = 0;

                /* compute inner product of wall velocity and i-th lattice velocity */
                for (int j = 0; j < 3; ++j)	{
                    dotProd += LATTICEVELOCITIES[i][j] * wallVelocity[j];
                }

                computeDensity(getEl(collideField, coord_dest, 0, n), &density);

                /* Set boundary i-th lattice with respect to the formula for MOVING_WALL */
                *cell_ptr += 2.0 * LATTICEWEIGHTS[i] * dotProd * density / (C_S*C_S);
            }
        }
    }
}

/**
 * Set OUTFLOW condition
 */
void setOutflow(float * collideField,
                int * flagField,
                const float * const ro_ref,
                int * node,
                int * n) {
    int i, coord_dest[3], flag;
    float * cell_ptr;
    float feq[Q];
    float velocity[D];
    float * fluidCell;

    /* for each lattice */
    for (i = 0; i < Q; i++) {
        /* compute a cell where lattice is pointing*/
        coord_dest[0] = node[0] + LATTICEVELOCITIES[i][0];
        coord_dest[1] = node[1] + LATTICEVELOCITIES[i][1];
        coord_dest[2] = node[2] + LATTICEVELOCITIES[i][2];

        /* does the pointed cell lay in our domain? */
        if (coord_dest[0] < n[0] && coord_dest[1] < n[1] && coord_dest[2] < n[2] &&
            coord_dest[0] >= 0 && coord_dest[1] >= 0 && coord_dest[2] >= 0 ) {
            flag = *getFlag(flagField, coord_dest, n);
       
            if (flag == FLUID || flag == INTERFACE) {
                /* get pointer to the fluid cell */
                fluidCell = getEl(collideField, coord_dest, 0, n);

                /* compute velocity of the fluid cell */
                computeVelocity(fluidCell, ro_ref, velocity);

                /* compute f-equilibrium of the fluid cell */
                computeFeq(ro_ref, velocity, feq);

                /* pointer to the i-th lattice of the boundary cell */
                cell_ptr = getEl(collideField, node, i, n);

                /* set boundary */
                *cell_ptr = feq[Q - i -1] + feq[i] - fluidCell[Q - 1 - i];
            }
        }
    }
}

/**
 * Set INFLOW condition
 */
void setInflow(float * collideField,
               int * flagField,
               const char * const scenario,
               const float * const Re,
               const float * const ro_ref,
               const float * const ro_in,
               const float * const inVelocity,
               int * node,
               int * n) {
    int i, coord_dest[3], flag;
    float * cell_ptr;
    float feq[Q];
    float velocity[3];

    /* If scenario is parabolic */
    if (strcmp(scenario, PARABOLIC_SCENARIO) == 0) {
        velocity[0] = 0;
        velocity[1] = 0;
        velocity[2] = - 0.5 * (*Re) * (*ro_in - *ro_ref) / n[0] * node[0] * (node[0] - n[2]);
    } else {
        velocity[0] = inVelocity[0];
        velocity[1] = inVelocity[1];
        velocity[2] = inVelocity[2];
    }

    /* for each lattice */
    for (i = 0; i < Q; i++) {
        /* compute a cell where lattice is pointing*/
        coord_dest[0] = node[0] + LATTICEVELOCITIES[i][0];
        coord_dest[1] = node[1] + LATTICEVELOCITIES[i][1];
        coord_dest[2] = node[2] + LATTICEVELOCITIES[i][2];

        /* does the pointed cell lay in our domain? */
        if (coord_dest[0] < n[0] && coord_dest[1] < n[1] && coord_dest[2] < n[2] &&
            coord_dest[0] >= 0 && coord_dest[1] >= 0 && coord_dest[2] >= 0) {
            flag = *getFlag(flagField, coord_dest, n); 
       
            if (flag == FLUID || flag == INTERFACE) {
            
                computeFeq(ro_ref, velocity, feq);

                cell_ptr = getEl(collideField, node, i, n);

                *cell_ptr = feq[i];
            }
        }
    }
}

void setFreeSlip(float * collideField, int * flagField, int * node, int * n) {
    int i, j, k, coord_dest[3], non_fluid_cell[3], flag;
    float * cell_ptr, sum, lv0, lv1, lv2;

    for (i = 0; i < Q; i++) {
        /* Initialize the cell with a flag, that will later make possible to know if some lattice was modifided*/
        *getEl(collideField, node, i, n) = 0;    
    }

    for (i = 0; i < Q; i++) {
        lv0 = LATTICEVELOCITIES[i][0];
        lv1 = LATTICEVELOCITIES[i][1];
        lv2 = LATTICEVELOCITIES[i][2];
        
        sum = lv0 * lv0 + lv1 * lv1 + lv2 * lv2;

        /* In this part we are interested only in the face of the cell, thus the lattice has just one component */
        if (sum == 1.0){
            coord_dest[0] = node[0] + LATTICEVELOCITIES[i][0];
            coord_dest[1] = node[1] + LATTICEVELOCITIES[i][1];
            coord_dest[2] = node[2] + LATTICEVELOCITIES[i][2];
            /* If the pointed cell does not fall out of bounds */
            if (coord_dest[0] < n[0] && coord_dest[1] < n[1] && coord_dest[2] < n[2] &&
                coord_dest[0] >= 0 && coord_dest[1] >= 0 && coord_dest[2] >= 0) {            
                flag = *getFlag(flagField, coord_dest, n);
                /* if pointed cell is FLUID */
                if (flag == FLUID || flag == INTERFACE) {    
                    for (j = 0; j < Q; j++) {
                        /* looking for a direction with one of the components inverse to the direction of the face */
                        if(LATTICEVELOCITIES[i][0]*LATTICEVELOCITIES[j][0] == -1.0 ||
                           LATTICEVELOCITIES[i][1]*LATTICEVELOCITIES[j][1] == -1.0 ||
                           LATTICEVELOCITIES[i][2]*LATTICEVELOCITIES[j][2] == -1.0) {
                            
                            /* If the selected direction of the fluid cell falls on another fluid cell, they will interact in the streaming step */
                            non_fluid_cell[0] = coord_dest[0] + LATTICEVELOCITIES[j][0];
                            non_fluid_cell[1] = coord_dest[1] + LATTICEVELOCITIES[j][1];
                            non_fluid_cell[2] = coord_dest[2] + LATTICEVELOCITIES[j][0];

                            flag = *getFlag(flagField, non_fluid_cell, n);
                            if (flag != FLUID && flag != INTERFACE) {
                                for (k = 0; k < Q; k++) {
                                    /* Search for a (unique) direcrion in the boundary cell which is the reflection of the fluid cell */
                                    if(( LATTICEVELOCITIES[k][0]*LATTICEVELOCITIES[i][0] == 1.0 || 
                                         LATTICEVELOCITIES[k][1]*LATTICEVELOCITIES[i][1] == 1.0 || 
                                         LATTICEVELOCITIES[k][2]*LATTICEVELOCITIES[i][2] == 1.0 ) &&
                                       ( LATTICEVELOCITIES[k][0]*LATTICEVELOCITIES[j][0] == 1.0 || 
                                         LATTICEVELOCITIES[k][1]*LATTICEVELOCITIES[j][1] == 1.0 || 
                                         LATTICEVELOCITIES[k][2]*LATTICEVELOCITIES[j][2] == 1.0 )) {
                                        cell_ptr = getEl(collideField, node, k, n);
                                        *cell_ptr = *getEl(collideField, coord_dest, j, n);
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    /* for each lattice */
    for (i = 0; i < Q; i++) {
        /*  If the lattice was not modifided in the previous process, this happens for the inverse direction of the latice going out of the face
        *   is also possible that the boundary does not share a fece with the fluid but it shares an edge or vertex.
        *   In those cases the boundary behabes as non slip, bouncing back everything*/
        if (*getEl(collideField, node, i, n) == 0){
            /* compute a cell where lattice is pointing*/
            coord_dest[0] = node[0] + LATTICEVELOCITIES[i][0];
            coord_dest[1] = node[1] + LATTICEVELOCITIES[i][1];
            coord_dest[2] = node[2] + LATTICEVELOCITIES[i][2];
    
            if (coord_dest[0] < n[0] && coord_dest[1] < n[1] && coord_dest[2] < n[2] &&
                coord_dest[0] >= 0 && coord_dest[1] >= 0 && coord_dest[2] >= 0) {            
                flag = *getFlag(flagField, coord_dest, n); 
                /* if pointed cell is FLUID */
                if (flag == FLUID || flag == INTERFACE) {
                    /* get pointer to the i-th lattice of boundary cell */
                    cell_ptr = getEl(collideField, node, i, n);
        
                    /* NOSLIP */
                    /* set i-th lattice to inverse lattice of the computed inner cell */
                    *cell_ptr= *getEl(collideField, coord_dest, Q-1-i, n);
                }
            }
        }
    }
}

void boundaryCell(float * collideField,
                  int * flagField,
                  const char * const scenario,
                  const float * const Re,
                  const float * const ro_ref,
                  const float * const ro_in,
                  const float * const velocity,
                  int * node,
                  int flag,
                  int * n) {

    if (flag == GAS) {
        return;
    } else if (flag == NOSLIP) {
        setNoSlip(collideField, flagField, node, n);
    } else if (flag == MOVING_WALL) {
        setMovingWall(collideField, flagField, velocity, node, n);
    } else if (flag == INFLOW) {
        setInflow(collideField, flagField, scenario, Re, ro_ref, ro_in, velocity, node, n);
    } else if (flag == OUTFLOW) {
        setOutflow(collideField, flagField, ro_ref, node, n);
    } else if (flag == FREESLIP) {
        setFreeSlip(collideField, flagField, node, n);
    } else if (flag == PRESSURE_IN) {
        setOutflow(collideField, flagField, ro_in, node, n);
    }
}

void treatBoundary(float *collideField,
                   int *flagField,
                   const char * const scenario,
                   const float * const Re,
                   const float * const ro_ref,
                   const float * const ro_in,
                   const float * const velocity,
                   int * length, 
                   int n_threads) { 
    int x, y, z, flag, node[3];
    int n[3] = { length[0] + 2, length[1] + 2, length[2] + 2 };

#pragma omp parallel for schedule(dynamic) collapse(2) private(node, x, flag) num_threads(n_threads)
    for (z = 0; z < n[2]; z++) {
        for (y = 0; y < n[1]; y++) {
            node[2] = z;
            node[1] = y;
            for (x = 0; x < n[0]; x++) {
                node[0] = x;
                flag = *getFlag(flagField, node, n);
                if (flag != FLUID && flag != OBSTACLE && flag != INTERFACE) {
                    boundaryCell(collideField, flagField, scenario, Re, ro_ref, ro_in, velocity, node, flag, n);
                }
            }
        }
    }
}
