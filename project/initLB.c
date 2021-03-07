#include "helper.h"
#include "initLB.h"
#include "LBDefinitions.h"

#include <regex.h>

/**
 * Reads the parameters for the lid driven cavity scenario from a config file.
 * Throws an error if number of program arguments does not equal 2.
 **/
int readParameters(Config &config, int argc, char *argv[])
{

    if (argc != 3) {
        printf("Usage: %s filename #threads\n",argv[0] );
        ERROR("number of arguments is incorrect");
    }
    config.n_threads = atoi(argv[2]);
    char path[80] = "examples/";

    strcat(path,argv[1]);
    const char *szFileName = strcat(path,".dat"); /* Concatenate .dat so the imput is independent of extesion*/

    read_int(szFileName, "zlength", &(config.length[2]));
    read_int(szFileName, "ylength", &(config.length[1]));
    read_int(szFileName, "xlength", &(config.length[0]));
    read_int(szFileName, "radius", &config.r);
    read_float(szFileName,"tau", &config.tau);
    read_float(szFileName,"exchange_factor", &config.exchange);

    read_float(szFileName, "velocity_x", &(config.velocity.x));
    read_float(szFileName, "velocity_y", &(config.velocity.y));
    read_float(szFileName, "velocity_z", &(config.velocity.z));

    read_float(szFileName, "forces_x", &(config.extForces[0]));
    read_float(szFileName, "forces_y", &(config.extForces[1]));
    read_float(szFileName, "forces_z", &(config.extForces[2]));

    read_int(szFileName, "timesteps", &config.timesteps);
    read_int(szFileName, "vtkoutput", &config.timestepsPerPlotting);
    read_string(szFileName, "problem", config.problem);
    read_float(szFileName,"ro_ref", &config.ro_ref);
    read_float(szFileName,"ro_in", &config.ro_in);

    read_int(szFileName, "wall_x0", &(config.boundaries[0]));
    read_int(szFileName, "wall_xmax", &(config.boundaries[1]));
    read_int(szFileName, "wall_y0", &(config.boundaries[2]));
    read_int(szFileName, "wall_ymax", &(config.boundaries[3]));
    read_int(szFileName, "wall_z0", &(config.boundaries[4]));
    read_int(szFileName, "wall_zmax", &(config.boundaries[5]));

    if (strcmp (config.problem, PARABOLIC_SCENARIO) != 0
        && strcmp (config.problem, CONSTANT_SCENARIO) != 0) {
        ERROR("Unrecognized scenario");
    }

    return 0;
}

/**
   Initialize cell flag and fluid distributions
 */
void initialiseCell(Fields &fields, int *n, int * node, int flag) {	
    int i;

    *getFlag (fields, node, n) = flag;    /*asigns the flag to the specified cell */

    for (i = 0; i < Q; ++i){
        *getEl(fields.stream, node, i, n) = LATTICEWEIGHTS[i];    /*Insert on each cell the initial value */
        *getEl(fields.collide, node, i, n) = LATTICEWEIGHTS[i];
    }
}

/*
  For FLUID cells set fluid fraction to 1.
  For INTERFACE cells sets initial mass to sum of the distributions from the neighboring FLUID cells.
  Set fraction to be equal to the mass, since our initial density is equal to 1.
 */
void initializeCellMass(Fields &fields, int * node, int * n) {
    int neighbor_node[3], i;

    auto flag = *getFlag (fields, node, n);
    if (flag == FLUID) {
        *getFraction (fields, node, n) = 1.0;
        return;
    }

    if (flag != INTERFACE) {
        return;
    }

    auto mass = getMass (fields, node, n);
    *mass = 0.0;

    for (i = 0; i < Q; i++) {
        neighbor_node[0] = node[0] + LATTICEVELOCITIES[i].x;
        neighbor_node[1] = node[1] + LATTICEVELOCITIES[i].y;
        neighbor_node[2] = node[2] + LATTICEVELOCITIES[i].z;

        if (neighbor_node[0] < n[0] && neighbor_node[1] < n[1] && neighbor_node[2] < n[2] &&
            neighbor_node[0] >= 0 && neighbor_node[1] >= 0 && neighbor_node[2] >= 0) {

            if (*getFlag (fields, neighbor_node, n) == FLUID) {
                *mass += *getEl(fields.collide, neighbor_node, Q - 1- i, n);
            }
        }
    }

    *getFraction (fields, node, n) = *mass;
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

/**
   Add droplet to the flag field.
   Droplet center will be in (n0/2, n1/2, n2 - r - 2), where r  is droplet radius.
 */
void initDropletFlags(Fields &fields, int * n, int r) {
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
                    *getFlag(fields, node, n) = FLUID;
                    /* Check whether its neighbors are in the sphere. If not mark them as INTERFACE */
                    for (i = 0; i < Q; i++) {
                        neighbor_node[0] = node[0] + LATTICEVELOCITIES[i].x;
                        neighbor_node[1] = node[1] + LATTICEVELOCITIES[i].y;
                        neighbor_node[2] = node[2] + LATTICEVELOCITIES[i].z;

                        if ((neighbor_node[0]-x0)*(neighbor_node[0]-x0) +
                            (neighbor_node[1]-y0)*(neighbor_node[1]-y0) +
                            (neighbor_node[2]-z0)*(neighbor_node[2]-z0) >= r*r) {
                            *getFlag (fields, node, n) = INTERFACE;
                        }
                    }
                }
            }
        }
    }
}

void initMartiniFlags(Fields &fields, int *n, int r) {
    int node[3], neighbor_node[3], i;

    int x0 = n[0] / 2;
    int y0 = n[1] / 2;
    int z0 = n[2] / 2;
    int x, y, z, sum, eps = 1;

    for (z = z0-r; z <= z0; z++) {
        node[2] = z;
        for (y = y0 - r; y <= y0 + r; y++) {
            node[1] = y;
            for (x = x0 - r; x <= x0 + r; x++) {
                node[0] = x;
                sum = (x-x0)*(x-x0) + (y-y0)*(y-y0) + (z-z0)*(z-z0);
                /* Initialize sphere of FLUID with radius r */
                if (sum <= (r+eps)*(r+eps) && sum >= (r-eps)*(r-eps)) {
                    *getFlag(fields, node, n) = OBSTACLE;
                }
            }
        }
    }

    for (z = z0-r; z <= z0; z++) {
        node[2] = z;
        for (y = y0 - r; y <= y0 + r; y++) {
            node[1] = y;
            for (x = x0 - r; x <= x0 + r; x++) {
                node[0] = x;

                if (*getFlag (fields, node, n) != OBSTACLE) {
                    continue;
                }

                /* Check whether its neighbors are in the sphere. If not mark them as NOSLIP */
                for (i = 0; i < Q; i++) {
                    neighbor_node[0] = node[0] + LATTICEVELOCITIES[i].x;
                    neighbor_node[1] = node[1] + LATTICEVELOCITIES[i].y;
                    neighbor_node[2] = node[2] + LATTICEVELOCITIES[i].z;
                    auto flag = getFlag(fields, neighbor_node, n);

                    if (*flag != GAS) {
                        continue;
                    }

                    *flag = NOSLIP;
                }
            }
        }
    }
}

/*
  Initialize flags from the image and set initial distributions for stream and collide field.
 */
void initialiseFlagsAndDF(Fields &fields, int ** image, int * length, int * n, int * walls, double * num_fluid_cells) {
    int x, y, z, node[3];

   /* 
    *   Definition of the fields 
    *   The structure of the looping makes possible to define the boundaries
    */

    /* z = 0; */
    node[2] = 0;
    for (y = 0; y <= length[1] + 1; ++y){
        node[1] = y;
        for (x = 0; x <= length[0] + 1; ++x) {
            node[0] = x;
            initialiseCell(fields, n, node, walls[4]);   /* Loop for the down wall of the cavity*/
        }
    }

    for (z = 1; z <= length[2]; ++z){
        node[2] = z;
        /* y = 0; */
        node[1] = 0;
        for (x = 0; x <= length[0] + 1; ++x){
            node[0] = x;
            initialiseCell(fields, n, node, walls[2]);   /* Loop for the front wall of the cavity*/
        }

        for (y = 1; y <= length[1]; ++y){
            node[1] = y;
            /* x = 0; */
            node[0] = 0;
            initialiseCell(fields, n, node, walls[0]);   /* Loop for the left wall of the cavity*/
            for (x = 1; x <= length[0]; ++x){
                node[0] = x;
                initialiseCell(fields, n, node, image[x][z]);           /* Loop for the interior points*/
                if (image[x][z] != GAS)
                    (*num_fluid_cells) ++;
            }
            /* x = length[0]+1; */
            node[0] = length[0] + 1;
            initialiseCell(fields, n, node, walls[1]);   /* Loop for the right wall of the cavity*/
        }

        /* y = length[1]+1; */
        node[1] = length[1] + 1;
        for (x = 0; x <= length[0] + 1; ++x){
            node[0] = x;
            initialiseCell(fields, n, node, walls[3]);   /* Loop for the back wall of the cavity*/
        }
    }

    /* z = length[2] + 1; */
    node[2] = length[2] + 1;
    for (y = 0; y <= length[1] + 1; ++y){
        node[1] = y;
        for (x = 0; x < length[0] + 2; ++x){
            node[0] = x;
            initialiseCell(fields, n, node, walls[5]);   /* Loop for the up wall of the cavity*/
        }
    }
}

/*
  Initializes collide, stream, flag, mass and fraction fields.
  Update number of fluid cells.
 */
void initialiseFields (Fields &fields, int * length, int * boundaries, int r, char *argv[], double * num_fluid_cells){
    int x, y, z, node[3], reti;
    int ** image;
    char path[80] = "examples/";
    int n[3] = { length[0] + 2, length[1] + 2, length[2] + 2 };
    regex_t regex;

    /*Copy the value from argv to filename to not modify the original one*/ 
    strcat(path, argv[1]);

    /* Concatenate .pgm so the imput is independent of extension*/
    strcat(path,".pgm");

    image = read_pgm(path);
    checkForbiddenPatterns(image, length);

    initialiseFlagsAndDF(fields, image, length, n, boundaries, num_fluid_cells);

    reti = regcomp(&regex, "droplet", 0);
    if (reti) {
        ERROR("could not compile regular expression");
    }

    /* Check whether it is droplet example */
    reti = regexec(&regex, argv[1], 0, NULL, 0);
    if ( ! reti) {
        initDropletFlags(fields, n, r);
    }

    initMartiniFlags(fields, n, r);

    /* Set initial mass and fraction for INTERFACE cells */
    for (z = 1; z <= length[2]; z++) {
        node[2] = z;
        for (y = 1; y <= length[1]; y++) {
            node[1] = y;
            for (x = 1; x <= length[0]; x++) {
                node[0] = x;

                initializeCellMass (fields, node, n);
            }
        }

    }

    free_imatrix(image, 0, length[2] + 2, 0, length[0] + 2);
}
