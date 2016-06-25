#include "LBDefinitions.h"
#include "computeCellValues.h"
#include <stdio.h>

void computeDensity(const double *const currentCell, double *density){
    /*
     * Computes the macroscopic density within currentCell
     * rho(x,t) = sum(f_i, i=0:Q-1)
     */

    int i;

    *density = 0;
    for (i=0; i<Q; i++) {
        *density += currentCell[i];
    }
}

void computeVelocity(const double * const currentCell, const double * const density, double *velocity){
    /*
     * Computes the velocity within currentCell 
     * u(x,t) = sum(f[i] * c[i] for i in [0:Q-1]) / rho
     */

    int i,j;

    for (j=0; j<3; j++) {
        velocity[j] = 0;
        for (i=0; i<Q; i++) {
            velocity[j] += LATTICEVELOCITIES[i][j] * currentCell[i];
        }
        velocity[j] /= *density;
    }
}

void computeFeq(const double * const density, const double * const velocity, double *feq){
    /*
     * feq[i] = w[i] * rho * (1 + c_i*u / (c_s^2) + ... 
     */

    int i;
    double ci_dot_u_cs2, u_dot_u_cs2;
    double ci0, ci1, ci2, u0, u1, u2;

    u0 = velocity[0];
    u1 = velocity[1];
    u2 = velocity[2];
    u_dot_u_cs2 = (u0*u0 + u1*u1 + u2*u2) / (2 * C_S2);
 
    for (i=0; i<Q; i++) {
        
        ci0 = LATTICEVELOCITIES[i][0];
        ci1 = LATTICEVELOCITIES[i][1];
        ci2 = LATTICEVELOCITIES[i][2];
       
        ci_dot_u_cs2 = ci0*u0 + ci1*u1 + ci2*u2 / C_S2;

        feq[i] = (1 + ci_dot_u_cs2 + (ci_dot_u_cs2 * ci_dot_u_cs2) / 2 - u_dot_u_cs2) * (*density) * LATTICEWEIGHTS[i];
    }
}

