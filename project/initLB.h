#ifndef _INITLB_H_
#define _INITLB_H_
#include "helper.h"
#include "utils.h"


/* reads the parameters for the lid driven cavity scenario from a config file */
int readParameters(Config &config,
                   int argc,                           /* number of arguments. Should equal 2 (program + name of config file */
                   char *argv[]);                      /* argv[1] shall contain the path to the config file */


/* initialises the particle distribution functions and the flagfield */
void initialiseFields(Fields &fields, int *length, int *boundaries, int r, char *argv[], double * num_fluid_cell);

/* Initializes oen cell in the fields, this is called by initialiseFields*/
void initialiseCell(Fields &fields, int *n, int * node, int flag);
#endif

