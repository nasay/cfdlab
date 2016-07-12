#include "helper.h"
#include "collision.h"
#include <omp.h>

/* Dotprod function, calcuales the dot product of two vectors */
float dotProd (const float * array1, float * array2, int length){
    int i;
    float product =0;
    
    for (i = 0; i < length; ++i){
        product += array1[i] * array2[i];
    }
    return product;
}

/* computeExternal : Calculates the influence of external forces in the lattices */
float computeExternal (int i, float density, float * extForces){
    return dotProd(LATTICEVELOCITIES[i], extForces, D) * density * LATTICEWEIGHTS[i];
}

/** computes the post-collision distribution functions according to the BGK update rule and
 *  stores the results again at the same position.
 */
void computePostCollisionDistributions(int *node, float * currentCell, int* flagField, float* fractionField, const float * const tau_inv, const float *const feq, const float *const feqAtm, float density, float * extForces, int * n){
    int i;
    float fi;

    for (i=0; i<Q; i++) {
        fi = currentCell[i];
        currentCell[i] = fi - (fi - feq[i]) * (*tau_inv) + computeExternal(i, density, extForces);
    }
}

void doCollision(float *collideField, int *flagField, float * massField, float * fractionField, const float * const tau, int * length, float * extForces, int n_threads){
    float density, densityAtm =1;
    float velocity[D];
    float feq[Q], feqAtm[Q];
    float * currentCell;
    int x, y, z, node[3], flag, isFluid;
    int n[3] = { length[0] + 2, length[1] + 2, length[2] + 2 };
    const float tau_inv = 1 / *tau;

    #pragma omp parallel for collapse(2) schedule(dynamic) private(x, node, density, feq, velocity, densityAtm, feqAtm, isFluid, flag, currentCell) num_threads(n_threads)
    // Loop over inner cells
    for (z = 1; z <= length[2]; z++) {
        for (y = 1; y <= length[1]; y++) {
            node[2] = z;
            node[1] = y;
            for (x = 1; x <= length[0]; x++) {
                node[0] = x;
                flag = *getFlag(flagField, node, n);
                isFluid = flag == FLUID;
                /* Perform the colision step just for fluid or interface cells */
                if (isFluid || flag == INTERFACE) {
                    /* Calculate the density velocity, and equilibrium distributions 
                    *  for the density calculated an the atmosferic density 
                    *  with all that data calculate the postcollision distributions for the cell*/
                    currentCell = getEl(collideField, node, 0, n);
                    computeDensity(currentCell, &density);
                    computeVelocity(currentCell, &density, velocity);
                    computeFeq(&density, velocity, feq);
                    computeFeq(&densityAtm, velocity, feqAtm);
                    computePostCollisionDistributions(node, currentCell, flagField, fractionField, &tau_inv, feq, feqAtm, density, extForces, n);
                    /* Update fluid fraction */
                    if (isFluid) {
                        *getMass(massField, node, n) = density;
                        *getFraction(fractionField, node, n) = 1; 
                    } else {
                        *getFraction(fractionField, node, n) = *getMass(massField, node, n) / density;
                    }
                }
            }
        }
    }
}
