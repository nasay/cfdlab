#include "helper.h"
#include "computeCellValues.h"
#include "LBDefinitions.h"
#include "checks.h"


void check_in_rank(Fields &fields, int * length, int t){
    int node[D], i, x, y, z;
    int n[D] = { length[0] + 2, length[1] + 2, length[2] + 2 };
    float density, norm_v=0;
    Velocity velocity;

    /* Pragma usted to parallelize the loop */
    #pragma omp parallel for private(node, y, x, density, norm_v, currentCell, velocity, i)
    /* For each cell in the cavity */
    for (z = 1; z < length[2]; ++z){
        node[2] = z;
        for (y = 1; y < length[1]; ++y){
            node[1] = y;
            for (x = 1; x < length[0]; ++x){
                node[0] = x;
                /* Just make the check if the cells are fluid */
                if (*getFlag (fields, node, n) == FLUID){
                    auto currentCell = getEl(fields.collide, node, 0, n);
                    computeDensity(currentCell, &density);
                    computeVelocity(currentCell, &density, velocity);
                    norm_v = 0;
                    /* Compute the norm of the velocity */
                        norm_v += velocity.x*velocity.x;
                        norm_v += velocity.y*velocity.y;
                        norm_v += velocity.z*velocity.z;
                    /* If the density is outside the normal values, trow a warning */
                    if(density>1.1 || density<0.9)
                        printf("Warning: In timestep %d position %d %d %d is an anormal density of %f \n", t, node[0],node[1],node[2],density);
                    /* If the velocity is grather than the speed of sound, trow a warning */
                    if(norm_v>(3 * C_S * C_S))
                        printf("Warning: In timestep %d position %d %d %d is an anormal velocity of %f \n", t, node[0],node[1],node[2],norm_v);
                }
            }
        }
    }
}

void check_flags(Fields &fields, int* length, int flag1, int flag2, int t){
    int node[D], node2[D], i, j, x, y, z;
    int n[D] = { length[0] + 2, length[1] + 2, length[2] + 2 };
    /* For each cell in the cavity */
    for (z = 1; z < length[2]; ++z){
        node[2] = z;
        for (y = 1; y < length[1]; ++y){
            node[1] = y;
            for (x = 1; x < length[0]; ++x){
                node[0] = x;
                /* Run the checks only if the flag coincides with flag1*/
                if (*getFlag (fields, node, n) == flag1){
                    /*For each neighbour cell */
                    for (i = 0; i < Q; ++i){
                      node2[0]=node[0]+LATTICEVELOCITIES[i].x;
                      node2[1]=node[1]+LATTICEVELOCITIES[i].y;
                      node2[2]=node[2]+LATTICEVELOCITIES[i].z;
                        /* If the cell flag coincides with flag2, trow a warning*/
                        if (*getFlag (fields, node, n) == flag2)
                            printf("Warning: In timestep %d position %d %d %d a cell %d is adjacent to a forbiden cell %d in %d %d %d \n", t, node[0],node[1],node[2],flag1, node2[0],node2[1],node2[2], flag2);
                    }
                }
            }
        }
    }
}

void check_mass(Fields &fields, int* length, int t){
    float tot_mass = 0;
    int node[D], x, y, z, flag;
    int n[D] = { length[0] + 2, length[1] + 2, length[2] + 2 };
    /* For each cell in the cavity */
    for (z = 1; z < length[2]; ++z){
        node[2] = z;
        for (y = 1; y < length[1]; ++y){
            node[1] = y;
            for (x = 1; x < length[0]; ++x){
                node[0] = x;
                flag = *getFlag(fields, node, n);
                /* If the cell is fluid or interface, add is mass to a global counter of the cavity */ 
                if (flag == FLUID || flag == INTERFACE) {
                    tot_mass += *getMass(fields, node, n);
                }
            }
        }
    }
    printf("On timestep %d : Total mass: %f \n",t,tot_mass);
}

void run_checks(Fields &fields, int * length, int t ){
    check_in_rank(fields, length, t);
    check_flags(fields, length, FLUID, GAS,   t);
    check_mass( fields, length, t);
}
