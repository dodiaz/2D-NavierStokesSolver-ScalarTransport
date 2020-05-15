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
    /* -------------------------------------------------------------------------
    -------------------------- Initializing Variables --------------------------
    ------------------------------------------------------------------------- */
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



    /* Initializing variables for this step 1 */
    double Hx;
    double Hy;
    double u_cc; /* cc = Cell cented velocity */
    double v_cc;
    double u_cc_im1; /* im1 = "i minus 1"; cell centered velocity at 1 position back in the i direction */
    double v_cc_jm1; /* jm1 = "j minus 1"; cell centered velocity at 1 position back in the j direction */
    double u_s;  /* s = staggered velocity (i-1/2, j-1/2) from the cc values */
    double v_s;
    double u_s_ip1; /* ip1 = "i plus 1"; staggered velocity at 1 position forward in the i direction */
    double u_s_jp1; /* jp1 = "j plus 1"; staggered velocity at 1 position forward in the j direction */
    double v_s_ip1;
    double v_s_jp1;

    


    /* -------------------------------------------------------------------------
    -------------------------- Initializing Arrays -----------------------------
    ------------------------------------------------------------------------- */

    double** u = (double**)calloc(ny, sizeof(double*));
    double** v = (double**)calloc(ny, sizeof(double*));
    double** u_star = (double**)calloc(ny, sizeof(double*));
    double** v_star = (double**)calloc(ny, sizeof(double*));

    double** Hx_old = (double**)calloc(ny, sizeof(double*));
    double** Hy_old = (double**)calloc(ny, sizeof(double*));

    double** p = (double**)calloc(ny, sizeof(double*));             //cell centered
    double** phi = (double**)calloc(ny, sizeof(double*));

    double** step1_mat_x = (double**)calloc(ny, sizeof(double*));
    double** step1_mat_y = (double**)calloc(ny, sizeof(double*));


    double* du_s = (double*)calloc(nx*ny, sizeof(double));
    double* dv_s = (double*)calloc(nx*ny, sizeof(double));
    double* du_ss = (double*)calloc(nx*ny, sizeof(double));
    double* dv_ss = (double*)calloc(nx*ny, sizeof(double));

    double* step1_mat_x_vec = (double*)calloc(nx*ny, sizeof(double));
    double* step1_mat_y_vec = (double*)calloc(nx*ny, sizeof(double));


    for (i = 0; i < ny; i++) {
        u[i] = (double*)calloc(nx, sizeof(double));
        v[i] = (double*)calloc(nx, sizeof(double));
        u_star[i] = (double*)calloc(nx, sizeof(double));
        v_star[i] = (double*)calloc(nx, sizeof(double));

        Hx_old[i] = (double*)calloc(nx, sizeof(double));
        Hy_old[i] = (double*)calloc(nx, sizeof(double));

        p[i] = (double*)calloc(nx, sizeof(double));
        phi[i] = (double*)calloc(nx, sizeof(double));

        step1_mat_x[i] = (double*)calloc(nx, sizeof(double));
        step1_mat_y[i] = (double*)calloc(nx, sizeof(double));
        

    }

    

    
    //sample for loop
    for(i = 0; i < Nx; i++) {
        for (j = 0; j < Ny; j++) {
            u[j][i] = 0;
            printf("(%d, %d)", j, i);
        }
        printf("\n");
    }



    /* -------------------------------------------------------------------------
    ----------------------------------- Step 1 ---------------------------------
    ------------------------------------------------------------------------- */

    for (j = 0; j < Ny; j++) {
        for (i = 0; i < Nx; i++) {
            /* Calculating step1_mat_x */
            if (i > 0 && j > 0 && i < (Nx - 1) && j < (Ny - 1)) {  /* At interior points except at the edges  */
                u_cc = (u[j][i] + u[j][i + 1]) / 2;
                u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

                u_s = (u[j][i] + u[j - 1][i]) / 2;
                v_s = (v[j][i] + v[j][i - 1]) / 2;
                u_s_jp1 = (u[j + 1][i] + u[j][i]) / 2;
                v_s_jp1 = (v[j + 1][i] + v[j + 1][i - 1]) / 2;

                Hx = ((pow(u_cc, 2) - pow(u_cc_im1, 2)) / dx) + (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

                step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * ((u[j][i + 1] - 2 * u[j][i] + u[j][i - 1]) / pow(dx, 2) + (u[j + 1][i] - 2 * u[j][i] + u[j - 1][i]) / pow(dy, 2)));
            }
            else if (i == 0) {      /* At the left wall  */
                u_cc = (u[j][i] + u[j][i + 1]) / 2;
                u_cc_im1 = U_inlet;

                u_s = U_inlet;
                v_s = 0;
                u_s_jp1 = (u[j + 1][i] + u[j][i]) / 2;
                v_s_jp1 = 0;

                Hx = ((pow(u_cc, 2) - pow(u_cc_im1, 2)) / dx) + (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

                step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * ((u[j][i + 1] - 2 * u[j][i] + u[j][i]) / pow(dx, 2) + (u[j + 1][i] - 2 * u[j][i] + u[j - 1][i]) / pow(dy, 2)));
            }
            else if (i == Nx - 1) {     /* At the right wall  */
                u_cc = u[j][i];
                u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

                u_s = (u[j][i] + u[j - 1][i]) / 2;
                v_s = (v[j][i] + v[j][i - 1]) / 2;
                u_s_jp1 = (u[j + 1][i] + u[j][i]) / 2;
                v_s_jp1 = (v[j + 1][i] + v[j + 1][i - 1]) / 2;

                Hx = ((pow(u_cc, 2) - pow(u_cc_im1, 2)) / dx) + (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

                step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * ((u[j][i] - 2 * u[j][i] + u[j][i - 1]) / pow(dx, 2) + (u[j + 1][i] - 2 * u[j][i] + u[j - 1][i]) / pow(dy, 2)));
            }
            else if (j == 0) {      /* At the bottom wall  */
                u_cc = (u[j][i] + u[j][i + 1]) / 2;
                u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

                u_s = u[j][i];     //this could also be set to zero
                v_s = (v[j][i] + v[j][i - 1]) / 2;
                u_s_jp1 = (u[j + 1][i] + u[j][i]) / 2;
                v_s_jp1 = (v[j + 1][i] + v[j + 1][i - 1]) / 2;

                Hx = ((pow(u_cc, 2) - pow(u_cc_im1, 2)) / dx) + (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

                step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * ((u[j][i + 1] - 2 * u[j][i] + u[j][i - 1]) / pow(dx, 2) + (u[j + 1][i] - 2 * u[j][i] + u[j][i]) / pow(dy, 2)));
            }
            else if (j == Ny - 1) {    /* At the top lid  */
                u_cc = (u[j][i] + u[j][i + 1]) / 2;
                u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

                u_s = (u[j][i] + u[j - 1][i]) / 2;
                v_s = (v[j][i] + v[j][i - 1]) / 2;
                u_s_jp1 = 0;
                v_s_jp1 = 0;

                Hx = ((pow(u_cc, 2) - pow(u_cc_im1, 2)) / dx) + (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

                step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * ((u[j][i + 1] - 2 * u[j][i] + u[j][i - 1]) / pow(dx, 2) + (u[j][i] - 2 * u[j][i] + u[j - 1][i]) / pow(dy, 2)));
            }






            /* Calculating step1_mat_y */
            if (i > 0 && j > 0 && i < (Nx - 1) && j < (Ny - 1)) {  /* At interior points except at the edges  */
                v_cc = (v[j][i] + v[j + 1][i]) / 2;
                v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

                u_s = (u[j][i] + u[j - 1][i]) / 2;
                v_s = (v[j][i] + v[j][i - 1]) / 2;
                u_s_ip1 = (u[j][i + 1] + u[j - 1][i + 1]) / 2;
                v_s_ip1 = (v[j][i + 1] + v[j][i]) / 2;

                Hy = ((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / dy) + (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

                step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) * ((v[j][i + 1] - 2 * v[j][i] + v[j][i - 1]) / pow(dx, 2) + (v[j + 1][i] - 2 * v[j][i] + v[j - 1][i]) / pow(dy, 2)));
            }
            else if (j == 0) {    /* At the bottom wall  */
                v_cc = (v[j][i] + v[j + 1][i]) / 2;
                v_cc_jm1 = 0;

                u_s = 0;
                v_s = (v[j][i] + v[j][i - 1]) / 2;
                u_s_ip1 = 0;
                v_s_ip1 = 0;

                Hy = ((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / dy) + (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

                step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) * ((v[j][i + 1] - 2 * v[j][i] + v[j][i - 1]) / pow(dx, 2) + (v[j + 1][i] - 2 * v[j][i] + v[j][i]) / pow(dy, 2)));
            }
            else if (j == Ny - 1) {     /* At the top wall  */
                v_cc = (v[j][i] + 0) / 2;
                v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

                u_s = (u[j][i] + u[j - 1][i]) / 2;
                v_s = (v[j][i] + v[j][i - 1]) / 2;
                u_s_ip1 = (u[j][i + 1] + u[j - 1][i + 1]) / 2;
                v_s_ip1 = (v[j][i + 1] + v[j][i]) / 2;

                Hy = ((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / dy) + (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

                step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) * ((v[j][i + 1] - 2 * v[j][i] + v[j][i - 1]) / pow(dx, 2) + (v[j][i] - 2 * v[j][i] + v[j - 1][i]) / pow(dy, 2)));
            }
            else if (i == 0) {      /* At left wall  */
                v_cc = (v[j][i] + v[j + 1][i]) / 2;
                v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

                u_s = U_inlet;
                v_s = 0;
                u_s_ip1 = (u[j][i + 1] + u[j - 1][i + 1]) / 2;
                v_s_ip1 = (v[j][i + 1] + v[j][i]) / 2;

                Hy = ((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / dy) + (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

                step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) * ((v[j][i + 1] - 2 * v[j][i] + v[j][i]) / pow(dx, 2) + (v[j + 1][i] - 2 * v[j][i] + v[j - 1][i]) / pow(dy, 2)));
            }
            else if (i == Nx - 1) {     /* At the right wall  */
                v_cc = (v[j][i] + v[j + 1][i]) / 2;
                v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

                u_s = (u[j][i] + u[j - 1][i]) / 2;
                v_s = (v[j][i] + v[j][i - 1]) / 2;
                u_s_ip1 = u_s;
                v_s_ip1 = v_s;

                Hy = ((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / dy) + (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

                step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) * ((v[j][i] - 2 * v[j][i] + v[j][i - 1]) / pow(dx, 2) + (v[j + 1][i] - 2 * v[j][i] + v[j - 1][i]) / pow(dy, 2)));
            }
        }
    }

    //deal with corner points for step1_mat_x and step1_mat_y


    //vectorize step1_mat_x and step1_mat_y


    //tridiagonal solve for du_ss and dv_ss


    //tridiagonal solve for du_s and dv_s


    //solve for u_star



    /* -------------------------------------------------------------------------
    ---------------------------------- Step 2 ----------------------------------
    ------------------------------------------------------------------------- */









    /* -------------------------------------------------------------------------
    ---------------------------------- Step 3 ----------------------------------
    ------------------------------------------------------------------------- */










    /* -------------------------------------------------------------------------
    ----------------------------- Scalar Transport -----------------------------
    ------------------------------------------------------------------------- */



    /* -------------------------------------------------------------------------
    ---------------------------- Freeing Variables -----------------------------
    ------------------------------------------------------------------------- */
    for (i = 0; i < Ny; i++) {
        free(u[i]);
        free(v[i]);
        free(u_star[i]);
        free(v_star[i]);

        free(Hx_old[i]);
        free(Hy_old[i]);

        free(p[i]);
        free(phi[i]);

        free(step1_mat_x[i]);
        free(step1_mat_y[i]);
    }

    free(u);
    free(v);
    free(u_star);
    free(v_star);

    free(Hx_old);
    free(Hy_old);

    free(p);
    free(phi);

    free(step1_mat_x);
    free(step1_mat_y);


    free(du_s);
    free(dv_s);
    free(du_ss);
    free(dv_s);

    free(step1_mat_x_vec);
    free(step1_mat_y_vec);

}