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

void computeVelocity(std::vector<float>::iterator currentCell, const float * const density, Velocity &velocity){
    /*
     * Computes the velocity within currentCell 
     * u(x,t) = sum(f[i] * c[i] for i in [0:Q-1]) / rho
     */

    int i;
    float density_inv = 1 / (*density);

    velocity.x = 0;
    velocity.y = 0;
    velocity.z = 0;

    for (i=0; i<Q; i++) {
      velocity.x += LATTICEVELOCITIES[i].x * (*(currentCell + i));
      velocity.y += LATTICEVELOCITIES[i].y * (*(currentCell + i));
      velocity.z += LATTICEVELOCITIES[i].z * (*(currentCell + i));
    }

    velocity.x *= density_inv;
    velocity.y *= density_inv;
    velocity.z *= density_inv;
}

void computeFeq(const float * const density, const Velocity &velocity, float * const feq){
    /*
     * feq[i] = w[i] * rho * (1 + c_i*u / (c_s^2) + ...
     */

    int i;
    float ci_dot_u_cs2[Q];
    float ci0[Q], ci1[Q], ci2[Q], a1[Q], a2[Q], a3[Q];

    const float u0 = velocity.x;
    const float u1 = velocity.y;
    const float u2 = velocity.z;
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

