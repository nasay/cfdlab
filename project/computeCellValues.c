#include "LBDefinitions.h"
#include "utils.h"
#include "computeCellValues.h"
#include <stdio.h>

void computeDensity(std::vector<float>::iterator currentCell, float *density){
    /*
     * Computes the macroscopic density within currentCell
     * rho(x,t) = sum(f_i, i=0:Q-1)
     */

    int i;

    *density = 0;
    for (i=0; i<Q; i++) {
      *density += (*(currentCell + i));
    }
}

void computeVelocity(std::vector<float>::iterator currentCell, const float * const density, float *velocity){
    /*
     * Computes the velocity within currentCell 
     * u(x,t) = sum(f[i] * c[i] for i in [0:Q-1]) / rho
     */

    int i;
    float density_inv = 1 / (*density);

    velocity[0] = 0;
    velocity[1] = 0;
    velocity[2] = 0;

    for (i=0; i<Q; i++) {
      velocity[0] += LATTICEVELOCITIES[i].x * (*(currentCell + i));
      velocity[1] += LATTICEVELOCITIES[i].y * (*(currentCell + i));
      velocity[2] += LATTICEVELOCITIES[i].z * (*(currentCell + i));
    }

    velocity[0] *= density_inv;
    velocity[1] *= density_inv;
    velocity[2] *= density_inv;
}

void computeFeq(const float * const density, const float * const velocity, float * const feq){
    /*
     * feq[i] = w[i] * rho * (1 + c_i*u / (c_s^2) + ...
     */

    int i;
    float ci_dot_u_cs2[Q];
    float ci0[Q], ci1[Q], ci2[Q], a1[Q], a2[Q], a3[Q];

    const float u0 = velocity[0];
    const float u1 = velocity[1];
    const float u2 = velocity[2];
    const float u_dot_u_cs2 = (u0*u0 + u1*u1 + u2*u2) * 0.5 * C_S2_inv;

#pragma GCC ivdep
    for (i=0; i<Q/2 + 1; i++) {
        ci0[i] = LATTICEVELOCITIES[i].x;
        ci1[i] = LATTICEVELOCITIES[i].y;
        ci2[i] = LATTICEVELOCITIES[i].z;
    }
#pragma GCC ivdep
    for (i=0; i<Q/2 + 1; i++) {
        a1[i] = ci0[i]*u0;
        a2[i] = ci1[i]*u1;
        a3[i] = ci2[i]*u2;
        ci_dot_u_cs2[i] = (a1[i] + a2[i] + a3[i]) * C_S2_inv;
        feq[i] = (1 + ci_dot_u_cs2[i] + (ci_dot_u_cs2[i] * ci_dot_u_cs2[i]) * 0.5 - u_dot_u_cs2) * (*density)* LATTICEWEIGHTS[i];
    }

#pragma GCC ivdep
    for (i = 0; i < Q / 2; i++) {
        feq[Q - i - 1] = feq[i] - 2 * ci_dot_u_cs2[i] * (*density) * LATTICEWEIGHTS[i];
    }

}

