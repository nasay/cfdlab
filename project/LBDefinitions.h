#include <math.h>

#ifndef _LBDEFINITIONS_H_
#define _LBDEFINITIONS_H_

#define Q 19
#define D 3

/* Math.sqrt(3.0) = 1.732050807568877 */
static const float C_S = 1.0 / 1.732050807568877;
static const float C_S2_inv = 3.0; 

static const int FLUID = 0;
static const int INTERFACE = 1;
static const int MOVING_WALL = 2;
static const int INFLOW = 3;
static const int OUTFLOW = 4;
static const int PRESSURE_IN = 5;
static const int OBSTACLE = 6;
static const int FREESLIP = 7;
static const int NOSLIP = 8;
static const int GAS = 9;

static const char * const PARABOLIC_SCENARIO = "parabolic";
static const char * const CONSTANT_SCENARIO = "constant";

#endif

