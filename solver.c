#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <string.h>





double* TriDiag_GaussElim(int size, double dx_or_dy, double dt, double Re, double* d, int right_outlet, int solving_for_du_s) {

    ////input size should be Nx
    ////should pass in d of whole row except for first point (should be size:(size - 1) )
    ////right_outlet = 1 if the right ourlet boundary is to the right of what we are looking at

    int i;

    double* b = (double*)calloc(size-1, sizeof(double));
    double* a = (double*)calloc(size-1, sizeof(double));
    double* c = (double*)calloc(size-1, sizeof(double));


    //initialize tridiagonal matrix
    for (i = 0; i < size - 1; i++) { 
        b[i] = 1 + dt / (Re * pow(dx_or_dy, 2) );
        a[i] = -dt / (2 * Re * pow(dx_or_dy,2) );
        c[i] = -dt / (2 * Re * pow(dx_or_dy,2) );
    }


    //if we have the outlet boundary in our domain
    if (right_outlet == 1)  {
        b[size - 2] = 1;
    }

    //if we are solving for du_s
    if (solving_for_du_s == 1) {
        b[0] = 1 + 3 * dt / (2 * Re * pow(dx_or_dy, 2));
    }



    //forward elimination
    for (i = 1; i < size - 1; i++) {
        b[i] = b[i] - c[i-1] * a[i] / b[i - 1];
        
        d[i] = d[i] - d[i-1] * a[i] / b[i-1]; 
    }


    //back substitution
    d[size - 2] = d[size - 2] / b[size - 2];

    for (i = size - 3; i > -1; i--) {
        d[i] = ( d[i] - c[i] * d[i+1] ) / b[i];
    }



    free(a);
    free(b);
    free(c);


    return d;

}
















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

    printf("\n \n");



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

    double** Hx_old = (double**)calloc(ny, sizeof(double*));
    double** Hy_old = (double**)calloc(ny, sizeof(double*));

    double** p = (double**)calloc(ny, sizeof(double*));             //cell centered
    double** phi = (double**)calloc(ny, sizeof(double*));

    double** step1_mat_x = (double**)calloc(ny, sizeof(double*));
    double** step1_mat_y = (double**)calloc(ny, sizeof(double*));


    double** du_s = (double**)calloc(ny, sizeof(double*));
    double** dv_s = (double**)calloc(ny, sizeof(double*));
    double** du_ss = (double**)calloc(ny, sizeof(double*));
    double** dv_ss = (double**)calloc(ny, sizeof(double*));


    double* temp_vec_x_small = (double*)calloc(99, sizeof(double));
    double* temp_vec_x_medium = (double*)calloc(19 * 20 - 1, sizeof(double));
    double* temp_vec_x_medium2 = (double*)calloc(19 * 20, sizeof(double));
    double* temp_vec_x_long = (double*)calloc(nx - 1, sizeof(double));

    double* temp_vec_y_small = (double*)calloc(40, sizeof(double));
    double* temp_vec_y_small2 = (double*)calloc(39, sizeof(double));
    double* temp_vec_y_long = (double*)calloc(ny, sizeof(double));
    double* temp_vec_y_long2 = (double*)calloc(ny - 1, sizeof(double));



    for (i = 0; i < Ny; i++) {
        u[i] = (double*)calloc(nx, sizeof(double));
        v[i] = (double*)calloc(nx, sizeof(double));

        Hx_old[i] = (double*)calloc(nx, sizeof(double));
        Hy_old[i] = (double*)calloc(nx, sizeof(double));

        p[i] = (double*)calloc(nx, sizeof(double));
        phi[i] = (double*)calloc(nx, sizeof(double));

        du_s[i] = (double*)calloc(nx, sizeof(double));
        dv_s[i] = (double*)calloc(nx, sizeof(double));
        du_ss[i] = (double*)calloc(nx, sizeof(double));
        dv_ss[i] = (double*)calloc(nx, sizeof(double));

        step1_mat_x[i] = (double*)calloc(nx, sizeof(double));
        step1_mat_y[i] = (double*)calloc(nx, sizeof(double));
        

    }



    //these variables must be a little larger for the multigrid to work
    double** u_star = (double**)calloc(ny, sizeof(double*));
    double** v_star = (double**)calloc(ny + 1, sizeof(double*));

    for (i = 0; i < Ny; i++) {
        u_star[i] = (double*)calloc(nx + 1, sizeof(double));
    }

    for (i = 0; i < Ny + 1; i++) {
        v_star[i] = (double*)calloc(nx, sizeof(double));
    }

    

    
    //sample for loop
    for(i = 0; i < Nx; i++) {
        for (j = 0; j < Ny; j++) {
            u[j][i] = 0;
            
        }
        
    }



    /* -------------------------------------------------------------------------
    ----------------------------------- Step 1 ---------------------------------
    ------------------------------------------------------------------------- */


    //--//--//--//--//-- Solve for step1_mat_x and step1_mat_y along boundaries and interior points --\\--\\--\\--\\--\\ 


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
            
                Hx_old[j][i] = Hx;
            }
            else if (i == 0 && j != 0 && j != Ny - 1) {      /* At the left wall  */ //////checked
                u_cc = (u[j][i] + u[j][i + 1]) / 2;
                u_cc_im1 = U_inlet;     ///technically don't use this value below

                u_s = U_inlet;
                v_s = 0;
                u_s_jp1 = (u[j + 1][i] + u[j][i]) / 2;
                v_s_jp1 = 0;

                Hx = ((pow(u_cc, 2) - pow(U_inlet, 2)) * 2 / dx) + (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

                step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * (          ( (u[j][i+1] - u[j][i]) / dx - (u_cc - U_inlet) * 2 / dx )/dx                    + (u[j + 1][i] - 2 * u[j][i] + u[j - 1][i]) / pow(dy, 2)));
            
                Hx_old[j][i] = Hx;
            }
            else if (i == Nx - 1 && j != 0 && j != Ny - 1) {     /* At the right wall  */  //////checked
                u_cc = u[j][i];
                u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

                u_s = (u[j][i] + u[j - 1][i]) / 2;
                v_s = (v[j][i] + v[j][i - 1]) / 2;
                u_s_jp1 = (u[j + 1][i] + u[j][i]) / 2;
                v_s_jp1 = (v[j + 1][i] + v[j + 1][i - 1]) / 2;

                Hx = ((pow(u_cc, 2) - pow(u_cc_im1, 2)) / dx) + (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

                step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * ((u[j][i] - 2 * u[j][i] + u[j][i - 1]) / pow(dx, 2) + (u[j + 1][i] - 2 * u[j][i] + u[j - 1][i]) / pow(dy, 2)));
            
                Hx_old[j][i] = Hx;
            }
            else if (j == 0&& i != 0 && i != Nx - 1) {      /* At the bottom wall  */   ///////checked 
                u_cc = (u[j][i] + u[j][i + 1]) / 2;
                u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

                u_s = 0;     
                v_s = 0;
                u_s_jp1 = (u[j + 1][i] + u[j][i]) / 2;
                v_s_jp1 = (v[j + 1][i] + v[j + 1][i - 1]) / 2;

                Hx = ((pow(u_cc, 2) - pow(u_cc_im1, 2)) / dx) + (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

                step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * ((u[j][i + 1] - 2 * u[j][i] + u[j][i - 1]) / pow(dx, 2) +    ( (u[j+1][i] - u[j][i]) / dy - (u[j][i] - u_s) *(2/dy) ) / dy   ));
            
                Hx_old[j][i] = Hx;
            }
            else if (j == Ny - 1 && i != 0 && i != Nx - 1) {    /* At the top wall  */    ///////checked
                u_cc = (u[j][i] + u[j][i + 1]) / 2;
                u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

                u_s = (u[j][i] + u[j - 1][i]) / 2;
                v_s = (v[j][i] + v[j][i - 1]) / 2;
                u_s_jp1 = 0;
                v_s_jp1 = 0;

                Hx = ((pow(u_cc, 2) - pow(u_cc_im1, 2)) / dx) + (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

                step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * ((u[j][i + 1] - 2 * u[j][i] + u[j][i - 1]) / pow(dx, 2) +         ( (u_s_jp1 - u[j][i]) *(2/dy) - (u[j][i] - u[j-1][i]) / dy) / dy        ));
            
                Hx_old[j][i] = Hx;
            }






            /* Calculating step1_mat_y */
            if (i > 0 && j > 0 && i < (Nx - 1) && j < (Ny - 1)) {  // At interior points except at the edges  
                v_cc = (v[j][i] + v[j + 1][i]) / 2;
                v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

                u_s = (u[j][i] + u[j - 1][i]) / 2;
                v_s = (v[j][i] + v[j][i - 1]) / 2;
                u_s_ip1 = (u[j][i + 1] + u[j - 1][i + 1]) / 2;
                v_s_ip1 = (v[j][i + 1] + v[j][i]) / 2;

                Hy = ((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / dy) + (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

                step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) * ((v[j][i + 1] - 2 * v[j][i] + v[j][i - 1]) / pow(dx, 2) + (v[j + 1][i] - 2 * v[j][i] + v[j - 1][i]) / pow(dy, 2)));
            
                Hy_old[j][i] = Hy;
            }
            else if (j == 0 && i != 0 && i != Nx - 1) {    /* At the bottom wall  */    /////corrected
                v_cc = (v[j][i] + v[j + 1][i]) / 2;
                v_cc_jm1 = 0;

                u_s = 0;
                v_s = 0;
                u_s_ip1 = 0;
                v_s_ip1 = 0;

                Hy = ((pow(v_cc, 2) - pow(v_cc_jm1, 2))*2 / dy) + (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

                step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) * ((v[j][i + 1] - 2 * v[j][i] + v[j][i - 1]) / pow(dx, 2) + ( (v[j+1][i] - v[j][i])/dy - (v_cc - v[j][i])*2/dy ) /dy        ));
            
                Hy_old[j][i] = Hy;
            }
            else if (j == Ny - 1 && i != 0 && i != Nx - 1) {     /* At the top wall  */    ///////corrected
                v_cc = (v[j][i] + 0) / 2;
                v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

                u_s = (u[j][i] + u[j - 1][i]) / 2;
                v_s = (v[j][i] + v[j][i - 1]) / 2;
                u_s_ip1 = (u[j][i + 1] + u[j - 1][i + 1]) / 2;
                v_s_ip1 = (v[j][i + 1] + v[j][i]) / 2;

                Hy = ((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / dy) + (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

                step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) * ((v[j][i + 1] - 2 * v[j][i] + v[j][i - 1]) / pow(dx, 2) +     (0 - 2 * v[j][i] + v[j - 1][i]) / pow(dy, 2)));
            
                Hy_old[j][i] = Hy;
            }
            else if (i == 0 && j != 0 && j != Ny - 1) {      /* At left wall  */  ///////corrected
                v_cc = (v[j][i] + v[j + 1][i]) / 2;
                v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

                u_s = U_inlet;
                v_s = 0;
                u_s_ip1 = (u[j][i + 1] + u[j - 1][i + 1]) / 2;
                v_s_ip1 = (v[j][i + 1] + v[j][i]) / 2;

                Hy = ((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / dy) + (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

                step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) *          (     ( (v[j][i+1] - v[j][i])/dx - (v[j][i] - v_s)*2/dx )  / dx         + (v[j + 1][i] - 2 * v[j][i] + v[j - 1][i]) / pow(dy, 2)));
            
                Hy_old[j][i] = Hy;
            }
            else if (i == Nx - 1 && j != 0 && j != Ny - 1) {     /* At the right wall  */    /////corrected
                v_cc = (v[j][i] + v[j + 1][i]) / 2;
                v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

                u_s = (u[j][i] + u[j - 1][i]) / 2;
                v_s = (v[j][i] + v[j][i - 1]) / 2;
                u_s_ip1 = u_s;
                v_s_ip1 = v_s;

                Hy = ((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / dy) + (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

                step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) * (     (v[j][i] - 2 * v[j][i] + v[j][i - 1]) / pow(dx, 2)        + (v[j + 1][i] - 2 * v[j][i] + v[j - 1][i]) / pow(dy, 2)));
            
                Hy_old[j][i] = Hy;
            }
        }
    }






    //--//--//--//--//-- deal with corner points for step1_mat_x and step1_mat_y --\\--\\--\\--\\--\\ 


    /////////////////top left corner/////////////////////
    i = 0;
    j = Ny - 1;

    //x data
    u_cc = (u[j][i] + u[j][i + 1]) / 2; 
    u_cc_im1 = U_inlet;     

    u_s = U_inlet;
    v_s = 0;
    u_s_jp1 = 0;
    v_s_jp1 = 0;

    Hx = ((pow(u_cc, 2) - pow(U_inlet, 2)) * 2 / dx) + (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

    step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * (          ( (u[j][i+1] - u[j][i]) / dx - (u_cc - U_inlet) * 2 / dx )/dx                    +  ( (u_s_jp1 - u[j][i]) *(2/dy) - (u[j][i] - u[j-1][i]) / dy) / dy    ));

    Hx_old[j][i] = Hx;
            

    //y data
    v_cc = (v[j][i] + 0) / 2;
    v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

    u_s = U_inlet;
    v_s = 0;
    u_s_ip1 = (u[j][i + 1] + u[j - 1][i + 1]) / 2;
    v_s_ip1 = (v[j][i + 1] + v[j][i]) / 2;

    Hy = ((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / dy) + (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

    step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) * (    ( (v[j][i+1] - v[j][i])/dx - (v[j][i] - v_s)*2/dx )  / dx +     (0 - 2 * v[j][i] + v[j - 1][i]) / pow(dy, 2)));

    Hy_old[j][i] = Hy;

    /////////////////top right corner/////////////////
    i = Nx - 1;
    j = Ny - 1;

    //x data
    u_cc = u[j][i];
    u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

    u_s = (u[j][i] + u[j - 1][i]) / 2;
    v_s = (v[j][i] + v[j][i - 1]) / 2;
    u_s_jp1 = 0;
    v_s_jp1 = 0;

    Hx = ((pow(u_cc, 2) - pow(u_cc_im1, 2)) / dx) + (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

    step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * (       (u[j][i] - 2 * u[j][i] + u[j][i - 1]) / pow(dx, 2)   +     ( (u_s_jp1 - u[j][i]) *(2/dy) - (u[j][i] - u[j-1][i]) / dy) / dy      ));

    Hx_old[j][i] = Hx;
            


    //y data
    v_cc = (v[j][i] + 0) / 2;
    v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

    u_s = (u[j][i] + u[j - 1][i]) / 2;
    v_s = (v[j][i] + v[j][i - 1]) / 2;
    u_s_ip1 = u_s;
    v_s_ip1 = v_s;

    Hy = ((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / dy) + (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

    step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) * (      (v[j][i] - 2 * v[j][i] + v[j][i - 1]) / pow(dx, 2)       +     (0 - 2 * v[j][i] + v[j - 1][i]) / pow(dy, 2)));

    Hy_old[j][i] = Hy;



    /////////////////bottom left corner/////////////////
    i = 0;
    j = 0;

    //x data
    u_cc = (u[j][i] + u[j][i + 1]) / 2;
    u_cc_im1 = U_inlet;     

    u_s = 0;
    v_s = 0;
    u_s_jp1 = (u[j + 1][i] + u[j][i]) / 2;
    v_s_jp1 = 0;

    Hx = ((pow(u_cc, 2) - pow(U_inlet, 2)) * 2 / dx) + (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

    step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * (          ( (u[j][i+1] - u[j][i]) / dx - (u_cc - U_inlet) * 2 / dx )/dx          +   ( (u[j+1][i] - u[j][i]) / dy - (u[j][i] - u_s) *(2/dy) ) / dy  ));

    Hx_old[j][i] = Hx;
            

    //y data
    v_cc = (v[j][i] + v[j + 1][i]) / 2;
    v_cc_jm1 = 0;

    u_s = 0;
    v_s = 0;
    u_s_ip1 = 0;
    v_s_ip1 = 0;

    Hy = ((pow(v_cc, 2) - pow(v_cc_jm1, 2))*2 / dy) + (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

    step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) * (    ( (v[j][i+1] - v[j][i])/dx - (v[j][i] - v_s)*2/dx )  / dx      + ( (v[j+1][i] - v[j][i])/dy - (v_cc - v[j][i])*2/dy ) /dy        ));

    Hy_old[j][i] = Hy;


    /////////////////bottom right corner/////////////////
    i = Nx - 1;
    j = 0;

    //x data
    u_cc = u[j][i];
    u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

    u_s = 0;
    v_s = 0;
    u_s_jp1 = (u[j + 1][i] + u[j][i]) / 2;
    v_s_jp1 = (v[j + 1][i] + v[j + 1][i - 1]) / 2;

    Hx = ((pow(u_cc, 2) - pow(u_cc_im1, 2)) / dx) + (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

    step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * (   (u[j][i] - 2 * u[j][i] + u[j][i - 1]) / pow(dx, 2)    +       ( (u[j+1][i] - u[j][i]) / dy - (u[j][i] - u_s) *(2/dy) ) / dy   ));

    Hx_old[j][i] = Hx;
            

    //y data
    v_cc = (v[j][i] + v[j + 1][i]) / 2;
    v_cc_jm1 = 0;

    u_s = 0;
    v_s = 0;
    u_s_ip1 = 0;
    v_s_ip1 = 0;

    Hy = ((pow(v_cc, 2) - pow(v_cc_jm1, 2))*2 / dy) + (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

    step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) * (    (v[j][i] - 2 * v[j][i] + v[j][i - 1]) / pow(dx, 2) + ( (v[j+1][i] - v[j][i])/dy - (v_cc - v[j][i])*2/dy ) /dy        ));

    Hy_old[j][i] = Hy;




    

    
    //--//--//--//--//-- Deal with boundary conditions for square obstacle --\\--\\--\\--\\--\\ 


    /////////////////////Top boundary of square/////////////////////////

    // u barycenter

    j = 60;
    for (i = 100; i < 120; i++) {
        u_cc = (u[j][i] + u[j][i + 1]) / 2;
        u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

        u_s = 0;     
        v_s = 0;
        u_s_jp1 = (u[j + 1][i] + u[j][i]) / 2;
        v_s_jp1 = (v[j + 1][i] + v[j + 1][i - 1]) / 2;

        Hx = ((pow(u_cc, 2) - pow(u_cc_im1, 2)) / dx) + (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

        step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * ((u[j][i + 1] - 2 * u[j][i] + u[j][i - 1]) / pow(dx, 2) +    ( (u[j+1][i] - u[j][i]) / dy - (u[j][i] - u_s) *(2/dy) ) / dy   ));
            
        Hx_old[j][i] = Hx;
    }

    /////////////////////Right boundary of square/////////////////////////

    // v barycenter
    i = 120;
    for (j = 40; j < 60; j++) {
        v_cc = (v[j][i] + v[j + 1][i]) / 2;
        v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

        u_s = 0;
        v_s = 0;
        u_s_ip1 = (u[j][i + 1] + u[j - 1][i + 1]) / 2;
        v_s_ip1 = (v[j][i + 1] + v[j][i]) / 2;

        Hy = ((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / dy) + (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

        step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) *          (     ( (v[j][i+1] - v[j][i])/dx - (v[j][i] - v_s)*2/dx )  / dx         + (v[j + 1][i] - 2 * v[j][i] + v[j - 1][i]) / pow(dy, 2)));
    
        Hy_old[j][i] = Hy;
    }

    /////////////////////Bottom boundary of square/////////////////////////

    //u barycenter
    j = 39;
    for (i = 100; i < 120; i++) {
        u_cc = (u[j][i] + u[j][i + 1]) / 2;
        u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

        u_s = (u[j][i] + u[j - 1][i]) / 2;
        v_s = (v[j][i] + v[j][i - 1]) / 2;
        u_s_jp1 = 0;
        v_s_jp1 = 0;

        Hx = ((pow(u_cc, 2) - pow(u_cc_im1, 2)) / dx) + (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

        step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * ((u[j][i + 1] - 2 * u[j][i] + u[j][i - 1]) / pow(dx, 2) +         ( (u_s_jp1 - u[j][i]) *(2/dy) - (u[j][i] - u[j-1][i]) / dy) / dy        ));
        
        Hx_old[j][i] = Hx;
              
    }

    //v barycenter
    j = 39;
    for (i = 100; i < 120; i++) {
        v_cc = (v[j][i] + 0) / 2;
        v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

        u_s = (u[j][i] + u[j - 1][i]) / 2;
        v_s = (v[j][i] + v[j][i - 1]) / 2;
        u_s_ip1 = (u[j][i + 1] + u[j - 1][i + 1]) / 2;
        v_s_ip1 = (v[j][i + 1] + v[j][i]) / 2;

        Hy = ((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / dy) + (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

        step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) * ((v[j][i + 1] - 2 * v[j][i] + v[j][i - 1]) / pow(dx, 2) +     (0 - 2 * v[j][i] + v[j - 1][i]) / pow(dy, 2)));
    
        Hy_old[j][i] = Hy;
    }



    /////////////////////Left boundary of square/////////////////////////

    //u barycenter
    i = 99;
    for (j = 40; j < 60; j++) {
        u_cc = (u[j][i] + 0 ) / 2;
        u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

        u_s = (u[j][i] + u[j - 1][i]) / 2;
        v_s = (v[j][i] + v[j][i - 1]) / 2;
        u_s_jp1 = (u[j + 1][i] + u[j][i]) / 2;
        v_s_jp1 = (v[j + 1][i] + v[j + 1][i - 1]) / 2;

        Hx = ((pow(u_cc, 2) - pow(u_cc_im1, 2)) / dx) + (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

        step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * ((u[j][i + 1] - 2 * u[j][i] + u[j][i - 1]) / pow(dx, 2) + (u[j + 1][i] - 2 * u[j][i] + u[j - 1][i]) / pow(dy, 2)));
    
        Hx_old[j][i] = Hx;
    }


    //v barycenter
    i = 99;
    for (j = 40; j < 60; j++) {
        v_cc = (v[j][i] + v[j + 1][i]) / 2;
        v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

        u_s = (u[j][i] + u[j - 1][i]) / 2;
        v_s = (v[j][i] + v[j][i - 1]) / 2;
        u_s_ip1 = 0;
        v_s_ip1 = 0;

        Hy = ((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / dy) + (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

        step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) * (     (v[j][i + 1] - 2 * v[j][i] + v[j][i - 1]) / pow(dx, 2)        + (v[j + 1][i] - 2 * v[j][i] + v[j - 1][i]) / pow(dy, 2)));
    
        Hy_old[j][i] = Hy;
    }




    
    


    
    //--//--//--//--//-- Tridiagonal solve for du_ss --\\--\\--\\--\\--\\ 


    //rows below the bottom of the square
    for (j = 0; j < 40; j++) {

        //fill in temp_vec with row of pixels except first pixel
        for(i = 0; i < Nx - 1; i++) {
            temp_vec_x_long[i] = step1_mat_x[j][i+1];
        }
 
        //gaussian elimination
        TriDiag_GaussElim(Nx, dx, dt, Re, temp_vec_x_long, 1, 0);

        //update du_ss
        for(i = 0; i < Nx - 1; i++) {
            du_ss[j][i+ 1] = temp_vec_x_long[i];
        }
        
    }


    //rows to left and right of square
    for (j = 40; j < 60; j++) {

        /////////left
        
        for(i = 0; i < 99; i++) {
            temp_vec_x_small[i] = step1_mat_x[j][i+1];
        }
        
        TriDiag_GaussElim(100, dx, dt, Re, temp_vec_x_small, 0, 0);

        for(i = 0; i < 99; i++) {
            du_ss[j][i+ 1] = temp_vec_x_small[i];
        }


        ////////right

        for(i = 0; i < 19 * 20 - 1; i++) {
            temp_vec_x_medium[i] = step1_mat_x[j][i+121];
        }
 
        TriDiag_GaussElim(19*20, dx, dt, Re, temp_vec_x_medium, 1, 0);

        for(i = 0; i < 19 * 20 - 1; i++) {
            du_ss[j][i+ 121] = temp_vec_x_medium[i];
        }
        

    }

    

    //rows above the top of the square
    for (j = 60; j < Ny; j++) {

        for(i = 0; i < Nx - 1; i++) {
            temp_vec_x_long[i] = step1_mat_x[j][i+1];
        }

        TriDiag_GaussElim(Nx, dx, dt, Re, temp_vec_x_long, 1, 0);

        for(i = 0; i < Nx - 1; i++) {
            du_ss[j][i+ 1] = temp_vec_x_long[i];
        }

    }

    
    

    //--//--//--//--//-- Tridiagonal solve for dv_ss --\\--\\--\\--\\--\\ 

    
    //rows below the bottom of the square
    for (j = 1; j < 40; j++) {

        //fill in temp_vec with row of pixels except first pixel
        for(i = 0; i < Nx - 1; i++) {
            temp_vec_x_long[i] = step1_mat_y[j][i+1];
        }
 
        //gaussian elimination
        TriDiag_GaussElim(Nx, dx, dt, Re, temp_vec_x_long, 1, 0);

        //update dv_ss
        for(i = 0; i < Nx - 1; i++) {
            dv_ss[j][i+ 1] = temp_vec_x_long[i];
        }
        
    }


    //rows to left and right of square
    for (j = 40; j < 60; j++) {

        /////////left
        
        for(i = 0; i < 99; i++) {
            temp_vec_x_small[i] = step1_mat_y[j][i+1];
        }
        
        TriDiag_GaussElim(100, dx, dt, Re, temp_vec_x_small, 0, 0);

        for(i = 0; i < 99; i++) {
            dv_ss[j][i+ 1] = temp_vec_x_small[i];
        }


        ////////right (this one was changed to incorporate the values at the boundary)

        for(i = 0; i < 19 * 20; i++) {
            temp_vec_x_medium2[i] = step1_mat_y[j][i+120];
        }
 
        TriDiag_GaussElim(19*20 + 1, dx, dt, Re, temp_vec_x_medium2, 1, 0);

        for(i = 0; i < 19 * 20; i++) {
            dv_ss[j][i+ 120] = temp_vec_x_medium2[i];
        }
        

    }


    //rows above the top of the square
    for (j = 60; j < Ny; j++) {

        for(i = 0; i < Nx - 1; i++) {
            temp_vec_x_long[i] = step1_mat_y[j][i+1];
        }

        TriDiag_GaussElim(Nx, dx, dt, Re, temp_vec_x_long, 1, 0);

        for(i = 0; i < Nx - 1; i++) {
            dv_ss[j][i+ 1] = temp_vec_x_long[i];
        }

    }






    
    //--//--//--//--//-- Tridiagonal solve for du_s --\\--\\--\\--\\--\\ 


    //columns to the left of the square
    for (i = 1; i < 100; i++) {

        for (j = 0; j < Ny; j++) {
            temp_vec_y_long[j] = du_ss[j][i];
        }

        TriDiag_GaussElim(Ny + 1, dy, dt, Re, temp_vec_y_long, 0, 1);

        for (j = 0; j < Ny; j++) {
            du_s[j][i] = temp_vec_y_long[j];
        }

    }


    

    //columns above and below the square

    for (i = 100; i < 120; i++) {

        //below
        for (j = 0; j < 40; j++) {
            temp_vec_y_small[j] = du_ss[j][i];
        }

        TriDiag_GaussElim(41, dy, dt, Re, temp_vec_y_small, 0, 1);

        for (j = 0; j < 40; j++) {
            du_s[j][i] = temp_vec_y_small[j];
        }

        //above
        for (j = 0; j < 40; j++) {
            temp_vec_y_small[j] = du_ss[j + 60][i];
        }

        TriDiag_GaussElim(41, dy, dt, Re, temp_vec_y_small, 0, 1);

        for (j = 0; j < 40; j++) {
            du_s[j + 60][i] = temp_vec_y_small[j];
        }

    }




    //columns to the right of the square
    for (i = 120; i < Nx; i++) {

        for (j = 0; j < Ny; j++) {
            temp_vec_y_long[j] = du_ss[j][i];
        }

        TriDiag_GaussElim(Ny + 1, dy, dt, Re, temp_vec_y_long, 0, 1);

        for (j = 0; j < Ny; j++) {
            du_s[j][i] = temp_vec_y_long[j];
        }

    }


    




    //--//--//--//--//-- Tridiagonal solve for dv_s --\\--\\--\\--\\--\\ 


    //columns to the left of the square
    for (i = 1; i < 100; i++) {

        for (j = 0; j < Ny - 1; j++) {
            temp_vec_y_long2[j] = dv_ss[j + 1][i];
        }

        TriDiag_GaussElim(Ny, dy, dt, Re, temp_vec_y_long2, 0, 0);

        for (j = 0; j < Ny - 1; j++) {
            dv_s[j + 1][i] = temp_vec_y_long2[j];
        }

    }


    //columns above and below the square

    for (i = 100; i < 120; i++) {

        //below
        for (j = 0; j < 39; j++) {
            temp_vec_y_small2[j] = dv_ss[j + 1][i];
        }

        TriDiag_GaussElim(40, dy, dt, Re, temp_vec_y_small2, 0, 0);

        for (j = 0; j < 39; j++) {
            dv_s[j + 1][i] = temp_vec_y_small2[j];
        }



        //above
        for (j = 0; j < 39; j++) {
            temp_vec_y_small2[j] = dv_ss[j + 61][i];
        }

        TriDiag_GaussElim(40, dy, dt, Re, temp_vec_y_small2, 0, 0);

        for (j = 0; j < 39; j++) {
            dv_s[j + 61][i] = temp_vec_y_small2[j];
        }


    }


    //columns to the right of the square
    for (i = 120; i < Nx; i++) {

        for (j = 0; j < Ny - 1; j++) {
            temp_vec_y_long2[j] = dv_ss[j + 1][i];
        }

        TriDiag_GaussElim(Ny, dy, dt, Re, temp_vec_y_long2, 0, 0);

        for (j = 0; j < Ny - 1; j++) {
            dv_s[j + 1][i] = temp_vec_y_long2[j];
        }

    }

    


    //--//--//--//--//-- Update u_star and v_star --\\--\\--\\--\\--\\ 

    for (j = 0; j < Ny; j++) {
        for (i = 0; i < Nx; i++) {
            u_star[j][i] = u[j][i] + du_s[j][i];
            v_star[j][i] = v[j][i] + dv_s[j][i];
        }
    }

    // update right boundary of u_star
    for (j = 0; j < Ny; j++) {
        u_star[j][Nx] = u_star[j][Nx - 1];
    }




    


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

        free(Hx_old[i]);
        free(Hy_old[i]);

        free(p[i]);
        free(phi[i]);

        free(step1_mat_x[i]);
        free(step1_mat_y[i]);

        free(du_s[i]);
        free(dv_s[i]);
        free(du_ss[i]);
        free(dv_ss[i]);

    }
    
    free(u);
    free(v);

    free(Hx_old);
    free(Hy_old);

    free(p);
    free(phi);

    free(du_s);
    free(dv_s);
    free(du_ss);
    free(dv_ss);

    free(step1_mat_x);
    free(step1_mat_y);

    free(temp_vec_x_small);
    free(temp_vec_x_medium);
    free(temp_vec_x_medium2);
    free(temp_vec_x_long);

    free(temp_vec_y_small);
    free(temp_vec_y_small2);
    free(temp_vec_y_long);
    free(temp_vec_y_long2);


    for (i = 0; i < Ny; i++) {
        free(u_star[i]);
    }

    for (i = 0; i < Ny + 1; i++) {
        free(v_star[i]);
    }

    free(u_star);
    free(v_star);

    printf("Wow the code finished working with no stopping errors, congrats!");
    

}