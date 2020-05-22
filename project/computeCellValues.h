#ifndef _COMPUTECELLVALUES_H_
#define _COMPUTECELLVALUES_H_
#include <vector>

/** computes the density from the particle distribution functions stored at currentCell.
 *  currentCell thus denotes the address of the first particle distribution function of the
 *  respective cell. The result is stored in density.
 */
void computeDensity(std::vector<float>::iterator currentCell, float *density);

/** computes the velocity within currentCell and stores the result in velocity */
void computeVelocity(std::vector<float>::iterator currentCell, const float * const density,float *velocity);

/** computes the equilibrium distributions for all particle distribution functions of one
 *  cell from density and velocity and stores the results in feq.
 */
void computeFeq(const float * const density, const float * const velocity, float *feq);

#endif

