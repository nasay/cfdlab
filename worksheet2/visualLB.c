#include "visualLB.h"
#include "computeCellValues.h"
#include "helper.h"

void write_vtkFile(const char *szProblem,
                   int    t,
                   double xlength,
                   double *collideField) {
  
    int x, y, z;
    char szFileName[80];
    FILE *fp=NULL;
    int n = xlength + 2;
    double velocity[3];
    double density;
    double * el = NULL;

    sprintf( szFileName, "%s.%i.vtk", szProblem, t );
    fp = fopen( szFileName, "w");
    if( fp == NULL )		       
    {
        char szBuff[80];
        sprintf( szBuff, "Failed to open %s", szFileName );
        ERROR( szBuff );
        return;
    }

    write_vtkHeader(fp, xlength);
    write_vtkPointCoordinates(fp, xlength);

    fprintf(fp,"POINT_DATA %i \n", n*n*n);//(n-2)*(n-2)*(n-2));
  
    fprintf(fp,"\n");
    fprintf(fp, "VECTORS velocity float\n");
    for(z = 0; z < n; z++) {
        for(y = 0; y < n; y++) {
            for(x = 0; x < n; x++) {
                /* TODO compute velocity */
                el = getEl(collideField, x, y, z, 0, n);
                computeDensity(el, &density);
                /* TODO */
                computeVelocity(el, &density, velocity);

                fprintf(fp, "%f %f %f\n", velocity[0], velocity[1], velocity[2]);
            }
        }
    }

    fprintf(fp,"\n");

    fprintf(fp, "SCALARS density float 1 \n"); 
    fprintf(fp, "LOOKUP_TABLE default \n");

    for(z = 0; z < n; z++) {
        for(y = 0; y < n; y++) {
            for(x = 0; x < n; x++) {
                /* TODO compute density */
                computeDensity(getEl(collideField, x, y, z, 0, n), &density);
                fprintf(fp, "%f\n", density);
            }
        }
    }

    if( fclose(fp) )
    {
        char szBuff[80];
        sprintf( szBuff, "Failed to close %s", szFileName );
        ERROR( szBuff );
    }
}


void write_vtkHeader( FILE *fp, int xlength) {
    int n = xlength + 2;

    if( fp == NULL ) {
        char szBuff[80];
        sprintf( szBuff, "Null pointer in write_vtkHeader" );
    
        ERROR( szBuff );
    
        return;
    }

    fprintf(fp,"# vtk DataFile Version 2.0\n");
    fprintf(fp,"generated by CFD-lab course output (written by Tobias Neckel) \n");
    fprintf(fp,"ASCII\n");
    fprintf(fp,"\n");	
    fprintf(fp,"DATASET STRUCTURED_GRID\n");
    fprintf(fp,"DIMENSIONS  %i %i %i \n", n, n, n);
    fprintf(fp,"POINTS %i float\n", n*n*n);//xlength*xlength*xlength);
    fprintf(fp,"\n");
}


void write_vtkPointCoordinates(FILE *fp, int xlength) {
    int originX = 0;
    int originY = 0;
    int originZ = 0;

    int x = 0;
    int y = 0;
    int z = 0;
    int n = xlength + 2;

    for(z = 0; z < n; z++) {
        for(y = 0; y < n; y++) {
            for(x = 0; x < n; x++) {
                /* dx = dy = dz = 1 */
                fprintf(fp, "%d %d %d\n", originX + x, originY + y, originZ + z);
            }
        }
    }
}

void writeVtkOutput(double * collideField, const int * const flagField, const char * filename, unsigned int t, int xlength) {
    write_vtkFile("ws2", t, xlength, collideField);
}

