#ifndef _MAIN_C_
#define _MAIN_C_

#include "collision.h"
#include "streaming.h"
#include "initLB.h"
#include "visualLB.h"
#include "boundary.h"
#include "flag.h"
#include "utils.h"
#include "checks.h"

#include <chrono>

int main(int argc, char *argv[]){
    float *swap = nullptr;
    Config config;

    std::chrono::duration<double> total_time {0.0};
    std::chrono::duration<double> elapsed_time;
    double mlups;

    char problem[10];

    #ifdef DEBUG
    std::chrono::duration<double> streamingTime {0.0};
    std::chrono::duration<double> collisionTime {0.0};
    std::chrono::duration<double> boundaryTime {0.0};
    std::chrono::duration<double> flagTime {0.0};
    #endif

    /* Read the file of parameters */
    readParameters (config, argc, argv);

    /* Domain size */
    int n[3] = { config.length[0] + 2, config.length[1] + 2, config.length[2] + 2 };

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

    for (int i = 0; i < (n[0] * n[1] * n[2]); ++i){
        filledCells[i] = (int *) malloc((size_t)( 3 * sizeof( int )));
        emptiedCells[i] = (int *) malloc((size_t)( 3 * sizeof( int )));
    }

    float viscosity = C_S * C_S * (config.tau - 0.5);
    float Re = 1 / viscosity;

    if (collideField == NULL || streamField == NULL || flagField == NULL ||
        massField == NULL || fractionField == NULL ||
        filledCells == NULL || emptiedCells == NULL) {
        ERROR("Unable to allocate matrices.");
    }

    initialiseFields(collideField, streamField, flagField,
                     massField, fractionField,
                     config.length, config.boundaries, config.r, argv, & num_fluid_cells);

    #ifdef DEBUG
    boundaryTime -= std::chrono::steady_clock::now ();
    #endif
    treatBoundary(collideField, flagField, problem, &Re, &config.ro_ref, &config.ro_in, config.velocity, config.length, config.n_threads);
    #ifdef DEBUG
    boundaryTime += std::chrono::steady_clock::now ();
    #endif

    // Start the timer for the lattice updates.
    auto start_time = std::chrono::steady_clock::now();

    for (int t = 0; t < config.timesteps; t++) {

        /* Streaming step */
        #ifdef DEBUG
        streamingTime -= std::chrono::steady_clock::now ();
        #endif
        doStreaming(collideField, streamField, flagField, massField, fractionField, config.length, config.n_threads, config.exchange);
        #ifdef DEBUG
        streamingTime += std::chrono::steady_clock::now ();
        #endif

        swap = collideField;
        collideField = streamField;
        streamField = swap;

        /* Collision step */
        #ifdef DEBUG
        collisionTime -= std::chrono::steady_clock::now ();
        #endif
        doCollision(collideField, flagField, massField, fractionField, &config.tau, config.length, config.extForces, config.n_threads);
        #ifdef DEBUG
        collisionTime += std::chrono::steady_clock::now ();
        #endif

        /* Updating flags */
        #ifdef DEBUG
        flagTime -= std::chrono::steady_clock::now ();
        #endif
        updateFlagField(collideField, flagField, fractionField, filledCells, emptiedCells, config.length, config.n_threads);
        #ifdef DEBUG
        flagTime += std::chrono::steady_clock::now ();
        #endif

        /* Updating boundaries */
        #ifdef DEBUG
        boundaryTime -= std::chrono::steady_clock::now ();
        #endif
        treatBoundary(collideField, flagField, problem, &Re, &config.ro_ref, &config.ro_in, config.velocity, config.length, config.n_threads);
        #ifdef DEBUG
        boundaryTime += std::chrono::steady_clock::now ();
        #endif

        if (t % config.timestepsPerPlotting == 0) {
            total_time += std::chrono::steady_clock::now () - start_time ; // Add elapsed ticks to total_time
            #ifdef DEBUG
                run_checks(collideField, massField, flagField, config.length, t );
            #endif
            writeVtkOutput(collideField, flagField, argv[1], t, config.length);
            printf("Time step %i finished, vtk file was created\n", t);
            start_time = std::chrono::steady_clock::now ();  // Start the timer for the lattice updates
        }
    }

    /* Compute average mega-lattice-updates-per-second in order to judge performance */
    elapsed_time = total_time;
    mlups = num_fluid_cells * config.timesteps / (elapsed_time.count () * 1000000);
    printf("Average MLUPS = %f\nElapsed time (excluding vtk writes) = %10.2f\n", mlups, elapsed_time);

    #ifdef DEBUG
    printf("Streaming time: %10.2f\n", streamingTime);
    printf("Collide time:   %10.2f\n", collisionTime);
    printf("Flag time:      %10.2f\n", flagTime);
    printf("Boundary time:  %10.2f\n", boundaryTime);
    #endif

    /* Free allocated memory */

    for (int i = 0; i < (n[0] * n[1] * n[2]); ++i){
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

