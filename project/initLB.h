#ifndef _INITLB_H_
#define _INITLB_H_
#include "helper.h"
#include "utils.h"


/* reads the parameters for the lid driven cavity scenario from a config file */
int readParameters(Config &config,
                   int argc,                           /* number of arguments. Should equal 2 (program + name of config file */
                   char *argv[]);                      /* argv[1] shall contain the path to the config file */


/* initialises the particle distribution functions and the flagfield */
void initialiseFields(float *collideField, float *streamField,int *flagField, float * massField, float * fractionField, int *length, int * boundaries, int r, char *argv[], double * num_fluid_cell);

/* Initializes oen cell in the fields, this is called by initialiseFields*/
void initialiseCell(float *collideField, float *streamField, int *flagField, int *n, int * node, int flag);
#endif

