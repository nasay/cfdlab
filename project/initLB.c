#include "helper.h"
#include "initLB.h"
#include "LBDefinitions.h"


/**
 * Reads the parameters for the lid driven cavity scenario from a config file.
 * Throws an error if number of program arguments does not equal 2.
 **/
int readParameters(
    int *length,                        /* reads domain size. Parameter name: "xlength" */
    double *tau,                        /* relaxation parameter tau. Parameter name: "tau" */
    double *velocity,               /* velocity of the lid. Parameter name: "characteristicvelocity" */
    double *extForces,                  /* External force, like gravity or electromagnetic */
    int *timesteps,                     /* number of timesteps. Parameter name: "timesteps" */
    int *timestepsPerPlotting,          /* timesteps between subsequent VTK plots. Parameter name: "vtkoutput" */
    int argc,                           /* number of arguments. Should equal 2 (program + name of config file */
    char *argv[],                       /* argv[1] shall contain the path to the config file */
    char *problem,                      /* specifies the considered problem scenario (parabolic or constant inflow )*/
    double *ro_ref,                     /* reference density nomally set to 1 */
    double *ro_in,                      /* density of inflow/outflow */
    int *boundaries,                     /* definition of the type of boundaries on each one of the walls, for definitions see LBDefinitios.h*/
    int * r
    ){

    if (argc != 2) {
        ERROR("number of arguments is incorrect");
    }
    char filename[20] ;
    strcpy(filename,argv[1]);
    const char *szFileName = strcat( filename,".dat"); /* Concatenate .dat so the imput is independent of extesion*/

    read_int(szFileName, "zlength", &length[2]);
    read_int(szFileName, "ylength", &length[1]);
    read_int(szFileName, "xlength", &length[0]);
    read_int(szFileName, "radius", r);
    read_double(szFileName,"tau", tau);

    read_double(szFileName, "velocity_x", &velocity[0]);
    read_double(szFileName, "velocity_y", &velocity[1]);
    read_double(szFileName, "velocity_z", &velocity[2]);
    
    read_double(szFileName, "forces_x", &extForces[0]);
    read_double(szFileName, "forces_y", &extForces[1]);
    read_double(szFileName, "forces_z", &extForces[2]);
    
    read_int(szFileName, "timesteps", timesteps);
    read_int(szFileName, "vtkoutput", timestepsPerPlotting);
    read_string(szFileName, "problem", problem);
    read_double(szFileName,"ro_ref", ro_ref);
    read_double(szFileName,"ro_in", ro_in);

    read_int(szFileName, "wall_x0", &boundaries[0]);
    read_int(szFileName, "wall_xmax", &boundaries[1]);
    read_int(szFileName, "wall_y0", &boundaries[2]);
    read_int(szFileName, "wall_ymax", &boundaries[3]);
    read_int(szFileName, "wall_z0", &boundaries[4]);
    read_int(szFileName, "wall_zmax", &boundaries[5]);

    if (strcmp(problem, PARABOLIC_SCENARIO) != 0 && strcmp(problem, CONSTANT_SCENARIO) != 0) {
        ERROR("Unrecognized scenario");
    }

    return 0;
}

void initialiseCell(double *collideField, double *streamField, int *flagField, int *n, int * node, int flag) {	
    int i;

    *getFlag(flagField, node, n) = flag;    /*asigns the flag to the specified cell */

    for (i = 0; i < Q; ++i){
        *getEl(streamField, node, i, n) = LATTICEWEIGHTS[i];    /*Insert on each cell the initial value */
        *getEl(collideField, node, i, n) = LATTICEWEIGHTS[i];
    }
}

/*
  For INTERFACE cells sets initial mass to sum of the distributions from the neighboring FLUID cells.
  Set fraction to be equal to the mass, since our initial density is equal to 1.
 */
void initializeMass(double * collideField, int * flagField, double *massField, double * fractionField, int * node, int * n) {
    int neighbor_node[3], i;
    double * mass;

    mass = getMass(massField, node, n);
    *mass = 0.0;

    for (i = 0; i < Q; i++) {
        neighbor_node[0] = node[0] + LATTICEVELOCITIES[i][0];
        neighbor_node[1] = node[1] + LATTICEVELOCITIES[i][1];
        neighbor_node[2] = node[2] + LATTICEVELOCITIES[i][2];

        if (neighbor_node[0] < n[0] && neighbor_node[1] < n[1] && neighbor_node[2] < n[2] &&
            neighbor_node[0] >= 0 && neighbor_node[1] >= 0 && neighbor_node[2] >= 0) {            
           
            if (*getFlag(flagField, neighbor_node, n) == FLUID) {
                *mass += *getEl(collideField, neighbor_node, Q - 1- i, n);
            }
        }
    }

    *getFraction(fractionField, node, n) = *mass;
}

/**
 * Checks that in input image there was no too thin boundaries
 */
void checkForbiddenPatterns(int ** image, int * length) {
    int x, z;
    int c1, c2, c3, c4, result;

    for (x = 0; x <= length[0]; x++) {
        for (z = 0; z <= length[2]; z++) {
            c1 = image[x][z] == FLUID;
            c2 = image[x + 1][z] == FLUID;
            c3 = image[x][z + 1] == FLUID;
            c4 = image[x + 1][z + 1] == FLUID;
            result = c1 << 3 | c2 << 2 | c3 << 1 | c4;

            if (result == 6 || result == 9) {  // result == 0b0110 or 0b1001
                ERROR("forbidden boundary");
            }
        }
    }
}

void initDropletFlags(int * flagField, int * n, int r) {
    int x0 = n[0] / 2; 
    int y0 = n[1] / 2; 
    int z0 = n[2] - r - 2;

    int node[3], neighbor_node[3];

    int x, y, z, i;

    for (z = z0 - r; z <= z0 + r; z++) {
        node[2] = z;
        for (y = y0 - r; y <= y0 + r; y++) {
            node[1] = y;
            for (x = x0 - r; x <= x0 + r; x++) {
                node[0] = x;
                /* Initialize sphere of FLUID with radius r */
                if ((x-x0)*(x-x0) + (y-y0)*(y-y0) + (z-z0)*(z-z0) < r*r) {
                    *getFlag(flagField, node, n) = FLUID;
                    /* Check whether its neighbors are in the sphere. If not mark them as INTERFACE */
                    for (i = 0; i < Q; i++) {
                        neighbor_node[0] = node[0] + LATTICEVELOCITIES[i][0];
                        neighbor_node[1] = node[1] + LATTICEVELOCITIES[i][1];
                        neighbor_node[2] = node[2] + LATTICEVELOCITIES[i][2];

                        if ((neighbor_node[0]-x0)*(neighbor_node[0]-x0) +
                            (neighbor_node[1]-y0)*(neighbor_node[1]-y0) +
                            (neighbor_node[2]-z0)*(neighbor_node[2]-z0) >= r*r) {
                            *getFlag(flagField, node, n) = INTERFACE;
                        }
                    }
                }
            }
        }
    }
}

void initialiseFields(double *collideField, double *streamField, int *flagField, double * massField, double * fractionField, int * length, int * boundaries, int r, char *argv[]){
    int x, y, z, node[3];
    int ** image;
    char filename[20];
    int n[3] = { length[0] + 2, length[1] + 2, length[2] + 2 };

    strcpy(filename,argv[1]);                   /*Copy the value from argv to filename to not modify the original one*/ 
    strcat( filename,".pgm"); /* Concatenate .pgm so the imput is independent of extension*/
    image = read_pgm(filename);
    checkForbiddenPatterns(image, length);

    /* 
    *   Definition of the fields 
    *   The structure of the looping makes possible to define the boudaries
    */

    /* z = 0; */
    node[2] = 0;
    for (y = 0; y <= length[1] + 1; ++y){
        node[1] = y;
        for (x = 0; x <= length[0] + 1; ++x){			
            node[0] = x;
            initialiseCell(collideField, streamField, flagField, n, node, boundaries[4]);   /* Loop for the down wall of the cavity*/
        }
    }
    
    for (z = 1; z <= length[2]; ++z){
        node[2] = z;
        /* y = 0; */
        node[1] = 0;
        for (x = 0; x <= length[0] + 1; ++x){
            node[0] = x;
            initialiseCell(collideField, streamField, flagField, n, node, boundaries[2]);   /* Loop for the front wall of the cavity*/
        }

        for (y = 1; y <= length[1]; ++y){
            node[1] = y;
            /* x = 0; */
            node[0] = 0;
            initialiseCell(collideField, streamField, flagField, n, node, boundaries[0]);   /* Loop for the left wall of the cavity*/
            for (x = 1; x <= length[0]; ++x){				
                node[0] = x;
                initialiseCell(collideField, streamField, flagField, n, node, image[x][z]);           /* Loop for the interior points*/
            }
            /* x = length[0]+1; */
            node[0] = length[0] + 1;
            initialiseCell(collideField, streamField, flagField, n, node, boundaries[1]);   /* Loop for the right wall of the cavity*/
        }

        /* y = length[1]+1; */
        node[1] = length[1] + 1;
        for (x = 0; x <= length[0] + 1; ++x){
            node[0] = x;
            initialiseCell(collideField, streamField, flagField, n, node, boundaries[3]);   /* Loop for the back wall of the cavity*/
        }
    }

    /* z = length[2] + 1; */
    node[2] = length[2] + 1;
    for (y = 0; y <= length[1] + 1; ++y){
        node[1] = y;
        for (x = 0; x < length[0] + 2; ++x){
            node[0] = x;
            initialiseCell(collideField, streamField, flagField, n, node, boundaries[5]);   /* Loop for the up wall of the cavity*/
        }
    }
 
    if (strcmp(argv[1], "droplet") == 0) {
        initDropletFlags(flagField, n, r);
    }
   
    for (z = 1; z <= length[2]; z++) {
        node[2] = z;
        for (y = 1; y <= length[1]; y++) {
            node[1] = y;
            for (x = 1; x <= length[0]; x++) {
                node[0] = x;

                if (*getFlag(flagField, node, n) == INTERFACE) {
                    initializeMass(collideField, flagField, massField, fractionField, node, n);
                }
            }
        }
        
    }
    free_imatrix(image, 0, length[2] + 2, 0, length[0] + 2);
}

