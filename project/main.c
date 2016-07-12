#ifndef _MAIN_C_
#define _MAIN_C_

#include <omp.h>

#include "collision.h"
#include "streaming.h"
#include "initLB.h"
#include "visualLB.h"
#include "boundary.h"
#include "flag.h"
#include "checks.h"

int main(int argc, char *argv[]){
    int length[3], timesteps, timestepsPerPlotting, boundaries[6], r, t, n_threads, i;

    float tau, extForces[3], exchange;
    float velocity[3], ro_in, ro_ref;
    float *swap=NULL;

    double start_time, total_time = 0;
    double elapsed_time, mlups;

    char problem[10];

    #ifdef DEBUG
    double streamingTime = 0.0;
    double collisionTime = 0.0;
    double boundaryTime = 0.0;
    double flagTime = 0.0;
    #endif


    /* Read the file of parameters */
    readParameters(length, &tau, velocity, extForces, &timesteps, &timestepsPerPlotting,
                   argc, argv, problem, &ro_ref, &ro_in, boundaries, &r, &n_threads, &exchange);

    /* Domain size */
    int n[3] = { length[0] + 2, length[1] + 2, length[2] + 2 };

    double num_fluid_cells = 2 * (n[0] + n[1] + n[2]);

    /* Collide field */
    float *collideField = (float *) malloc((size_t)(Q * n[0] * n[1] * n[2]) * sizeof(float));

    /* Streaming field */
    float *streamField = (float *) malloc((size_t)(Q * n[0] * n[1] * n[2]) * sizeof(float));

    /* Flag field */
    int *flagField = (int *) malloc((size_t)(n[0] * n[1] * n[2]) * sizeof(int));

    /* Mass field */
    float * massField = (float *) malloc((size_t)(n[0] * n[1] * n[2]) * sizeof(float));

    /* fluid fraction field */
    float * fractionField = (float *) malloc((size_t)(n[0] * n[1] * n[2]) * sizeof(float));

    /* Cells that were filled during current timestep */
    int **filledCells = (int **) malloc((size_t)(n[0] * n[1] * n[2] * sizeof( int * )));

    /* Cells that were filled during current timestep */
    int **emptiedCells = (int **) malloc((size_t)(n[0] * n[1] * n[2] * sizeof( int * )));

    for (i = 0; i < (n[0] * n[1] * n[2]); ++i){
        filledCells[i] = (int *) malloc((size_t)( 3 * sizeof( int )));
        emptiedCells[i] = (int *) malloc((size_t)( 3 * sizeof( int )));
    }

    float viscosity = C_S * C_S * (tau - 0.5);
    float Re = 1 / viscosity;

    if (collideField == NULL || streamField == NULL || flagField == NULL ||
        massField == NULL || fractionField == NULL ||
        filledCells == NULL || emptiedCells == NULL) {
        ERROR("Unable to allocate matrices.");
    }

    initialiseFields(collideField, streamField, flagField,
                     massField, fractionField,
                     length, boundaries, r, argv, & num_fluid_cells);

    #ifdef DEBUG
    boundaryTime -= omp_get_wtime();
    #endif
    treatBoundary(collideField, flagField, problem, &Re, &ro_ref, &ro_in, velocity, length, n_threads);
    #ifdef DEBUG
    boundaryTime += omp_get_wtime();
    #endif

    start_time = omp_get_wtime();  // Start the timer for the lattice updates

    for (t = 0; t < timesteps; t++) {

        /* Streaming step */
        #ifdef DEBUG
        streamingTime -= omp_get_wtime();
        #endif
        doStreaming(collideField, streamField, flagField, massField, fractionField, length, n_threads, exchange);
        #ifdef DEBUG
        streamingTime += omp_get_wtime();
        #endif

        swap = collideField;
        collideField = streamField;
        streamField = swap;

        /* Collision step */
        #ifdef DEBUG
        collisionTime -= omp_get_wtime();
        #endif
        doCollision(collideField, flagField, massField, fractionField, &tau, length, extForces, n_threads);
        #ifdef DEBUG
        collisionTime += omp_get_wtime();
        #endif

        /* Updating flags */
        #ifdef DEBUG
        flagTime -= omp_get_wtime();
        #endif
        updateFlagField(collideField, flagField, fractionField, filledCells, emptiedCells, length, n_threads);
        #ifdef DEBUG
        flagTime += omp_get_wtime();
        #endif

        /* Updating boundaries */
        #ifdef DEBUG
        boundaryTime -= omp_get_wtime();
        #endif
        treatBoundary(collideField, flagField, problem, &Re, &ro_ref, &ro_in, velocity, length, n_threads);
        #ifdef DEBUG
        boundaryTime += omp_get_wtime();
        #endif

        if (t % timestepsPerPlotting == 0) {
            total_time += omp_get_wtime() - start_time ; // Add elapsed ticks to total_time
            #ifdef DEBUG
                run_checks(collideField, massField, flagField, length, t );
            #endif
            writeVtkOutput(collideField, flagField, argv[1], t, length);
            printf("Time step %i finished, vtk file was created\n", t);
            start_time = omp_get_wtime();  // Start the timer for the lattice updates
        }
    }

    /* Compute average mega-lattice-updates-per-second in order to judge performance */
    elapsed_time = total_time; 
    mlups = num_fluid_cells * timesteps / (elapsed_time*1000000);
    printf("Average MLUPS = %f\nElapsed time (excluding vtk writes) = %10.2f\n", mlups, elapsed_time);

    #ifdef DEBUG
    printf("Streaming time: %10.2f\n", streamingTime);
    printf("Collide time:   %10.2f\n", collisionTime);
    printf("Flag time:      %10.2f\n", flagTime);
    printf("Boundary time:  %10.2f\n", boundaryTime);
    #endif

    /* Free allocated memory */

    for (i = 0; i < (n[0] * n[1] * n[2]); ++i){
        free(filledCells[i]);
        free(emptiedCells[i]);
    }

    free(filledCells);
    free(emptiedCells);

    free(collideField);
    free(streamField);
    free(flagField);
    free(massField);
    free(fractionField);

    return 0;
}

#endif

