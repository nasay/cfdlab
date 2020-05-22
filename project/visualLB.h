#ifndef _VISUALLB_H_
#define _VISUALLB_H_
#include <stdio.h>
#include <stdlib.h>
#include "helper.h"

/**
 * Method for writing header information in vtk format. 
 * 
 * The name of the file consists of the problem name (szProblem) 
 * and of the current time step. It gets the suffix .vtk. 
 * 
 * @param szProblem      File pointer for writing info.  
 * @param timeStepNumber Number of the current time step to be printed.  
 * @param xlength Length in x-direction
 * @param ylength Length in y-direction
 * @param imax    Maximum number of entries (?) in x-direction
 * @param jmax    Maximum number of entries (?) in y-direction
 * @param dx      Mesh size in x-direction
 * @param dy      Mesh size in x-direction
 * @param U       Velocities in x-direction
 * @param V       Velocities in y-direction
 * @param P       Pressure data
 * 
 * @author Tobias Neckel
 */
void write_vtkFile(const char *szProblem,
                   int    t,
                   int * length,
                   Fields &fields);

/**
 * Method for writing header information in vtk format. 
 * 
 * @param fp      File pointer for writing info.  
 * @param imax    Maximum number of entries (minus 2) in x-direction
 * @param jmax    Maximum number of entries (minus 2) in y-direction
 * @param dx      mesh size dx
 * @param dy      mesh size dy
 * 
 * @author Tobias Neckel
 */
void write_vtkHeader(FILE *fp, int * length);

/**
 * Method for writing grid coordinate information in vtk format. 
 * 
 * @param fp      File pointer for writing info.  
 * @param imax    Maximum number of entries (minus 2) in x-direction
 * @param jmax    Maximum number of entries (minus 2) in y-direction
 * @param dx      mesh size dx
 * @param dy      mesh size dy
 * 
 * @author Tobias Neckel
 */
void write_vtkPointCoordinates(FILE *fp, int * length); 

/** writes the density and velocity field (derived from the distributions in collideField)
 *  to a file determined by 'filename' and timestep 't'. You can re-use parts of the code
 *  from visual.c (VTK output for Navier-Stokes solver) and modify it for 3D datasets.
 */
void writeVtkOutput(Fields &fields,
                    const char * filename,
                    unsigned int t,
                    int * length);

#endif

