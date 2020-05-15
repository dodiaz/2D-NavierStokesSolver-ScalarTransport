#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <string.h>


/* ================================================================================
= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
======================================= MAIN ======================================
= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
================================================================================ */


int main() {
    /* --------------------- Initializing Variables  ------------------------------ */
    int i = 0, j = 0, k = 0; /* Counters for all "for loops" */
    int iter = 0;

    double Re = 20;  /* Problem parameters */
    double dt = 0.0001;
    double H = 1;
    double U_inlet = 1;
    double Diff = pow(10,-4);

    double epsilon = pow(10, -3);

    int Nx = 500;
    int Ny = 100;
    double nx = Nx;
    double ny = Ny;
    double dx = H / 20;
    double dy = H / 20;

    


    /* ---------------------- Initializing Arrays -------------------------------- */

    double** u = (double**)calloc(ny, sizeof(double*));
    double** v = (double**)calloc(ny, sizeof(double*));
    double** u_star = (double**)calloc(ny, sizeof(double*));
    double** v_star = (double**)calloc(ny, sizeof(double*));

    double** H_old = (double**)calloc(ny, sizeof(double*));

    double** p = (double**)calloc(ny, sizeof(double*));             //cell centered
    double** phi = (double**)calloc(ny, sizeof(double*));



    for (i = 0; i < nx; i++) {
        u[i] = (double*)calloc(nx, sizeof(double));
        v[i] = (double*)calloc(nx, sizeof(double));
        u_star[i] = (double*)calloc(nx, sizeof(double));
        v_star[i] = (double*)calloc(nx, sizeof(double));

        H_old[i] = (double*)calloc(nx, sizeof(double));

        p[i] = (double*)calloc(nx, sizeof(double));
        phi[i] = (double*)calloc(nx, sizeof(double));
        
        printf("%d \n", i);

    }

    

    
    //sample for loop
    for(i = 0; i < Nx; i++) {
        for (j = 0; j < Ny; j++) {
            u[j][i] = 0;
            printf("Hello world");
        }
    }
    



    /* ----------------------- Freeing Variables -------------------------------- */
    for (i = 0; i < Nx; i++) {
        free(u[i]);
        free(v[i]);
        free(u_star[i]);
        free(v_star[i]);

        free(H_old[i]);

        free(p[i]);
        free(phi[i]);
    }

    free(u);
    free(v);
    free(u_star);
    free(v_star);

    free(H_old);

    free(p);
    free(phi);

}
