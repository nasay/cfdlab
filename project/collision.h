#ifndef _COLLISION_H_
#define _COLLISION_H_

#include "computeCellValues.h"
#include "utils.h"

/** carries out the whole local collision process. Computes density and velocity and
 *  equilibrium distributions. Carries out BGK update.
 */
void doCollision(Fields &fields, const float * const tau,int * length, float * extForces, int n_threads);
#endif
