#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <string.h>


//function used to print the data to a text file
void print_current_data(int step, double** u, double** v, double** phi, int NX, int NY, char method[]) {
    int i, j;

    char filename1[40] = "step_"; // printing u data
    itoa(step, filename1 + 4, 10);
    strcat(filename1, "_");
    strcat(filename1, method);
    strcat(filename1, "_u_data.txt");

    FILE* fpointer1 = fopen(filename1, "w");
    for (j = 0; j < NY; j++) {
        for (i = 0; i < NX; i++) {
            fprintf(fpointer1, "%.20lf ", u[j][i]);
        }
        fprintf(fpointer1, "\n");
    }
    fclose(fpointer1);


    char filename2[25] = "step_"; // printing v data
    itoa(step, filename2 + 4, 10);
    strcat(filename2, "_");
    strcat(filename2, method);
    strcat(filename2, "_v_data.txt");

    FILE* fpointer2 = fopen(filename2, "w");
    for (j = 0; j < NY; j++) {
        for (i = 0; i < NX; i++) {
            fprintf(fpointer2, "%.20lf ", v[j][i]);
        }
        fprintf(fpointer2, "\n");
    }
    fclose(fpointer2);



    char filename3[30] = "step_"; // printing phi data
    itoa(step, filename3 + 4, 10);
    strcat(filename3, "_");
    strcat(filename3, method);
    strcat(filename3, "_phi_data.txt");

    FILE* fpointer3 = fopen(filename3, "w");
    for (j = 0; j < NY; j++) {
        for (i = 0; i < NX; i++) {
            fprintf(fpointer3, "%.20lf ", phi[j][i]);
        }
        fprintf(fpointer3, "\n");
    }
    fclose(fpointer3);

    return;
}








double* TriDiag_GaussElim(int size, double dx_or_dy, double dt, double Re, double* d, int right_outlet, int first_last_half, int first_half) {

    ////input size should be Nx
    ////should pass in d of whole row except for first point (should be size:(size - 1) )
    ////right_outlet = 1 if the right outlet boundary is to the right of what we are looking at
    ////first_last_half == 1 if the first and last datapoints in row/column used and there is a wall a half-step away from both endpoints with u* or u** = 0 at the wall
    ////first_half == 1 if the first datapoint in a row/column is used and there is a wall a half-step away from the point with the u* or u** = 0 at that wall

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
        b[size - 2] = 1 + dt / (2 * Re * pow(dx_or_dy,2));
    }

    //if we are solving for du_s
    if (first_last_half == 1) {
        b[0] = 1 + 3 * dt / (2 * Re * pow(dx_or_dy, 2));
        b[size - 2] = 1 + 3 * dt / (2 * Re * pow(dx_or_dy, 2));
    }

    //if the first datapoint in a row/column is used and there is a wall a half-step away from the point with the u* or u** value being zero at that wall
    if (first_half == 1) {
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




/* ==========================================================================
======================= n-Step Gauss-Seidel Solver ==========================
========================================================================== */

int GS_nstep(double** f, double** phi, int Nx, int Ny, double D_x, double D_y, double epsilon, int nGS, int void_len, int Bj, int Bi) {

    ///////// Initilaizations ////////
    int i, j;
    double nx = Nx;
    double ny = Ny;

    double f_norm;
    double integral;
    double lambda = pow(D_x, -2);  // for now, this only works for square meshes
    double RHS;    
    double laplace_phi_minus_f_norm;
   
    int step = 1;
    int max_num_steps = nGS;

    double** laplace_phi = (double**)calloc(Ny, sizeof(double*));
    for (j = 0; j < Ny; j++) {
        laplace_phi[j] = (double*)calloc(Nx, sizeof(double));
    }


    ///////// Solving ////////

    // compute f_norm 
    f_norm = 0;
    for (j = 0; j < Ny; j++) {
        for (i = 0; i < Nx; i++) {

            if (  j < Bj ||  j > Bj + void_len - 1  || i < Bi || i > Bi + void_len - 1  ) {     // if not inside the void 
                f_norm = f_norm + pow(f[j][i], 2);
            }

        }
    }
   
    f_norm = sqrt(f_norm);
     
    


    do {
        
        step += 1;

        /* ================================================= CALCULATING PHI ================================================= */

        ///// Loop through entire space 
        for (j = Ny - 1; j > -1; j--) {
            for (i = Nx - 1; i > -1; i--) {

                if (i > 0 && i < Nx - 1 && j > 0 && j < Ny - 1) { // points that are not along the outer edges
                    
                    if ( j >= Bj && j <= Bj + void_len - 1 && i >= Bi && i <= Bi + void_len - 1 ) { // inside the void
                        
                        
                    }
                    else if (i == Bi - 1 && j >= Bj && j <= Bj + void_len - 1) { // at the left boundary of the void
                        
                        phi[j][i] = (phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) / 3 - f[j][i] / (3 * lambda);
                        
                    }
                    else if (i == Bi + void_len && j >= Bj && j <= Bj + void_len - 1) { // right boundary of the void
                        
                        phi[j][i] = (phi[j][i + 1] + phi[j - 1][i] + phi[j + 1][i]) / 3 - f[j][i] / (3 * lambda);
                        
                    }
                    else if (j == Bj + void_len && i >= Bi && i <= Bi + void_len - 1) { // top boundary of the void
                        
                        phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j + 1][i]) / 3 - f[j][i] / (3 * lambda);
                        
                    }
                    else if (j == Bj - 1 && i >= Bi && i <= Bi + void_len - 1) { // bottom boundary of the void
                        
                        phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i]) / 3 - f[j][i] / (3 * lambda);
                        
                    }
                    else { // all other interior points not next to an edge
                        
                        phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) / 4 - f[j][i] / (4 * lambda);
                   
                    }

                }
                else if (j == 0 && i != 0 && i != Nx - 1) { // bottom edge
                    
                    phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j + 1][i]) / 3 - f[j][i] / (3 * lambda);
                                           
                }
                else if (j == Ny - 1 && i != 0 && i != Nx - 1) { // top edge
                    
                    phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i]) / 3 - f[j][i] / (3 * lambda);
                       
                }
                else if (i == 0 && j != 0 && j != Ny - 1) { // left edge 
                    
                    phi[j][i] = (phi[j][i + 1] + phi[j - 1][i] + phi[j + 1][i]) / 3 - f[j][i] / (3 * lambda);

                }
                else if (i == Nx - 1 && j != 0 && j != Ny - 1) { // right edge
                    
                    phi[j][i] = (phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) / 3 - f[j][i] / (3 * lambda);

                }
                else if (i == 0 && j == 0) { // bottom left corner 
                    
                    phi[0][0] = (phi[1][0] + phi[0][1]) / 2 - f[0][0] / (2 * lambda);

                }
                else if (i == Nx - 1 && j == 0) { // bottom right corner
                    
                    phi[0][Nx - 1] = (phi[1][Nx - 1] + phi[0][Nx - 2]) / 2 - f[0][Nx - 1] / (2 * lambda);

                }
                else if (i == 0 && j == Ny - 1) { // top left corner
                    
                    phi[Ny - 1][0] = (phi[Ny - 1][1] + phi[Ny - 2][0]) / 2 - f[Ny - 1][0] / (2 * lambda);

                }
                else if (i == Nx - 1 && j ==  Ny - 1) { //top right corner
                    
                    phi[Ny - 1][Nx - 1] = (phi[Ny - 1][Nx - 2] + phi[Ny - 2][Nx - 1]) / 2 - f[Ny - 1][Nx - 1] / (2 * lambda);

                }

                

            }

            
        }






        
        
        
    
        if (step % 5 == 2) {

            
            
            
            /* ================================================= CALCULATING LAPLACE PHI ================================================= */

            /* Loop through entire space */
            for (j = Ny - 1; j > -1; j--) {
                for (i = Nx - 1; i > -1; i--) {

                    if (i > 0 && i < Nx - 1 && j > 0 && j < Ny - 1) { // points that are not along the outer edges
                        
                        if ( j >= Bj && j <= Bj + void_len - 1 && i >= Bi && i <= Bi + void_len - 1 ) { // inside the void
                            
                        }
                        else if (i == Bi - 1 && j >= Bj && j <= Bj + void_len - 1) { // left boundary of void
                            
                            laplace_phi[j][i] = (phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - 3 * lambda * phi[j][i];
                        }
                        else if (i == Bi + void_len && j >= Bj && j <= Bj + void_len - 1) { // right boundary of void
                            
                            laplace_phi[j][i] = (phi[j][i + 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - 3 * lambda * phi[j][i];
                        }
                        else if (j == Bj + void_len && i >= Bi && i <= Bi + void_len - 1) { // top boundary of void
                            
                            laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j + 1][i]) * lambda - 3 * lambda * phi[j][i];
                            
                        }
                        else if (j == Bj - 1 && i >= Bi && i <= Bi + void_len - 1) { // bottom boundary of void
                            
                            laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i]) * lambda - 3 * lambda * phi[j][i];
                        }
                        else { // general interior points not next to a boundary
                            
                            laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - 4 * lambda * phi[j][i];
                        
                        }

                    }
                    else if (j == 0 && i != 0 && i != Nx - 1) { // bottom edge
                        
                        laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j + 1][i]) * lambda - 3 * lambda * phi[j][i];

                    }
                    else if (j == Ny - 1 && i != 0 && i != Nx - 1) { // top edge
                        
                        laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i]) * lambda - 3 * lambda * phi[j][i];

                    }
                    else if (i == 0 && j != 0 && j != Ny - 1) { // left edge
                        
                        laplace_phi[j][i] = (phi[j][i + 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - 3 * lambda * phi[j][i];

                    }
                    else if (i == Nx - 1 && j != 0 && j != Ny - 1) { // right edge
                       
                        laplace_phi[j][i] = (phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - 3 * lambda * phi[j][i];

                    }
                    else if (i == 0 && j == 0) { // bottom left corner
                        
                        laplace_phi[j][i] = (phi[j][i + 1] + phi[j + 1][i]) * lambda - 2 * lambda * phi[j][i];

                    }
                    else if (i == Nx - 1 && j == 0) { // bottom right corner
                        
                        laplace_phi[j][i] = (phi[j][i - 1] + phi[j + 1][i]) * lambda - 2 * lambda * phi[j][i];

                    }
                    else if (i == 0 && j == Ny - 1) { // top left corner
                        
                        laplace_phi[j][i] = (phi[j][i + 1] + phi[j - 1][i]) * lambda - 2 * lambda * phi[j][i];

                    }
                    else if (i == Nx - 1 && j ==  Ny - 1) { // top right corner
                        
                        laplace_phi[j][i] = (phi[j][i - 1] + phi[j - 1][i]) * lambda - 2 * lambda * phi[j][i];

                    }

                }
            }  


            /* compute the norm */
            laplace_phi_minus_f_norm = 0;

            for (j = 0; j < Ny; j++) {
                for (i = 0; i < Nx; i++) {

                    if ( j < Bj ||  j > Bj + void_len - 1 || i < Bi || i > Bi + void_len - 1  ) {     /* if not inside the void... */
                        laplace_phi_minus_f_norm += pow((laplace_phi[j][i] - f[j][i]), 2);
                    }

                }
            }


            laplace_phi_minus_f_norm = sqrt(laplace_phi_minus_f_norm);



        }

        



        if (f_norm == 0) {
            RHS = epsilon;
        }
        else {
            RHS = epsilon * f_norm;
        }

        

        // break out of the loop if the max number of steps has been reached 
        if (step == max_num_steps) {

            break;
        }

        
        

    } while (laplace_phi_minus_f_norm > RHS);



    // Impose condition that the integral over the domain is equal to zero 
    integral = 0;
    for (j = 0; j < Ny; j++) {
        for (i = 0; i < Nx; i++) {
           
            if ( j < Bj ||  j > Bj + void_len - 1 || i < Bi || i > Bi + void_len - 1  ) {     // if not inside the void
                integral += phi[j][i] * D_x * D_y;
            }

        }
    }

    for (j = 0; j < Ny; j++) {
        for (i = 0; i < Nx; i++) {

            if ( j < Bj ||  j > Bj + void_len - 1 || i < Bi || i > Bi + void_len - 1  ) {     // if not inside the void
                phi[j][i] -= integral / (Nx * Ny - pow(void_len, 2));
            }
       
        }
    }

   

    // Freeing arrays made in this function 
    for (i = 0; i < Ny; i++) {
        free(laplace_phi[i]);
    }

    free(laplace_phi);

    return 0;

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

    printf("\n");

    // Changeable parameters (parameters that change the problem we are solving) 
    double Re = 400;  
    double H = 1;
    double U_inlet = 1;
    double Diff = pow(10,-4);
    
    int Nx = 500;
    int Ny = 100;
    double nx = Nx;
    double ny = Ny;
    
    int H_mesh = 20;
    double dx = H / H_mesh;
    double dy = H / H_mesh;

    int Bj = 40; // y location of the bottom corner of the void
    int Bi = 100; // x location of the bottom corner of the void
    int void_len = 20; // represents how many cells are in the void (void_len * void_len)

    

    // parameters that define the methods we use to solve the NS equations
    char convective_method[] = "quick";  // options are "upwind", "centraldiff", or "quick"
    int max_time_steps = 100;
    double dt = 0.001;
    double epsilon = pow(10, -3);
    double nGS = 500; // Number of Gauss-Seidel steps to take 

    

    // counter for all "for loops"
    int i = 0, j = 0, k = 0; 
    int iter = 0;



    // Initializing variables for this step 1 
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



    // For scalar transport 
    double convec;
    double diffu;



    /* -------------------------------------------------------------------------
    -------------------------- Initializing Arrays -----------------------------
    ------------------------------------------------------------------------- */

    // barycentered arrays
    double** u = (double**)calloc(Ny, sizeof(double*));
    double** v = (double**)calloc(Ny, sizeof(double*));

    double** Hx_old = (double**)calloc(Ny, sizeof(double*));
    double** Hy_old = (double**)calloc(Ny, sizeof(double*));

    double** step1_mat_x = (double**)calloc(Ny, sizeof(double*));
    double** step1_mat_y = (double**)calloc(Ny, sizeof(double*));

    double** du_s = (double**)calloc(Ny, sizeof(double*));
    double** dv_s = (double**)calloc(Ny, sizeof(double*));
    double** du_ss = (double**)calloc(Ny, sizeof(double*));
    double** dv_ss = (double**)calloc(Ny, sizeof(double*));

    double* temp_vec_x_small = (double*)calloc(Bi - 1, sizeof(double));                    
    double* temp_vec_x_small2 = (double*)calloc(Bi, sizeof(double));
    double* temp_vec_x_medium = (double*)calloc(Nx - Bi - void_len - 1, sizeof(double));    
    double* temp_vec_x_medium2 = (double*)calloc(Nx - Bi - void_len, sizeof(double));       
    double* temp_vec_x_long = (double*)calloc(Nx - 1, sizeof(double));     
    double* temp_vec_x_long2 = (double*)calloc(Nx, sizeof(double));           

    double* temp_vec_y_small = (double*)calloc(Bj, sizeof(double));
    double* temp_vec_y_small2 = (double*)calloc(Bj - 1, sizeof(double));
    double* temp_vec_y_small3 = (double*)calloc(Ny - Bj - void_len, sizeof(double));
    double* temp_vec_y_small4 = (double*)calloc(Ny - Bj - void_len - 1, sizeof(double));
    double* temp_vec_y_long = (double*)calloc(Ny, sizeof(double));
    double* temp_vec_y_long2 = (double*)calloc(Ny - 1, sizeof(double));


    // cell centered arrays
    double** p = (double**)calloc(Ny, sizeof(double*));             
    double** phi = (double**)calloc(Ny, sizeof(double*));
    double** phi_new = (double**)calloc(Ny, sizeof(double*));
    double** grad_u_star_over_dt = (double**)calloc(Ny, sizeof(double*));



    for (i = 0; i < Ny; i++) {
        u[i] = (double*)calloc(Nx, sizeof(double));
        v[i] = (double*)calloc(Nx, sizeof(double));

        Hx_old[i] = (double*)calloc(Nx, sizeof(double));
        Hy_old[i] = (double*)calloc(Nx, sizeof(double));

        

        step1_mat_x[i] = (double*)calloc(Nx, sizeof(double));
        step1_mat_y[i] = (double*)calloc(Nx, sizeof(double));

        du_s[i] = (double*)calloc(Nx, sizeof(double));
        dv_s[i] = (double*)calloc(Nx, sizeof(double));
        du_ss[i] = (double*)calloc(Nx, sizeof(double));
        dv_ss[i] = (double*)calloc(Nx, sizeof(double));


        p[i] = (double*)calloc(Nx, sizeof(double));
        phi[i] = (double*)calloc(Nx, sizeof(double));
        phi_new[i] = (double*)calloc(Nx, sizeof(double));
        grad_u_star_over_dt[i] = (double*)calloc(Nx, sizeof(double));
    }


    //these variables must be a little larger for the multigrid to work
    double** u_star = (double**)calloc(Ny, sizeof(double*));
    double** v_star = (double**)calloc(Ny + 1, sizeof(double*));


    for (i = 0; i < Ny; i++) {
        u_star[i] = (double*)calloc(Nx + 1, sizeof(double));
    }

    for (i = 0; i < Ny + 1; i++) {
        v_star[i] = (double*)calloc(Nx, sizeof(double));
    }


    //Initialize velocity at the left boundary
    for (j = 0; j < Ny; j++) {
        u[j][0] = U_inlet;
        u_star[j][0] = U_inlet;    
    }



    



    for (iter = 0; iter < max_time_steps; iter++) {

        /* -------------------------------------------------------------------------
        ----------------------------------- Step 1 ---------------------------------
        ------------------------------------------------------------------------- */


        /* ------------- Calculate step1_mat_x and step1_mat_y ------------ */
        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {

        
                ////////// Calculating step1_mat_x (at all but corner points)  /////////
                if (i >= Bi - 1 && i <= Bi + void_len + 1 && j >= Bj - 1 && j <= Bj + void_len) { //void and boundary points
                    
                    if (i >= Bi && i <= Bi + void_len && j >= Bj && j < Bj + void_len) { //points that are inside or at a wall, do nothing
                        
                    }
                    else if (i == Bi - 1 && j >= Bj - 1 && j <= Bj + void_len) { //left boundary of void (and top left and bottom left corner points)
                        
                        if (j == Bj - 1 || j == Bj + void_len) { //top and bottom left corner points
                            u_cc = (u[j][i] + u[j][i + 1]) / 2;
                        }
                        else {
                            u_cc = (u[j][i] + 0 ) / 2;
                        }
                        
                        u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

                        u_s = (u[j][i] + u[j - 1][i]) / 2;
                        v_s = (v[j][i] + v[j][i - 1]) / 2;
                        u_s_jp1 = (u[j + 1][i] + u[j][i]) / 2;
                        v_s_jp1 = (v[j + 1][i] + v[j + 1][i - 1]) / 2;

                        Hx = -((pow(u_cc, 2) - pow(u_cc_im1, 2)) / dx) - (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

                        step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * ((u[j][i + 1] - 2 * u[j][i] + u[j][i - 1]) / pow(dx, 2) + (u[j + 1][i] - 2 * u[j][i] + u[j - 1][i]) / pow(dy, 2)));
                
                        Hx_old[j][i] = Hx;
                        
                    }
                    else if (i == Bi + void_len + 1 && j >= Bj - 1 && j <= Bj + void_len) { //right boundary of void (plus extra corner points that are introduced)
                        u_cc = (u[j][i] + u[j][i + 1]) / 2;
                        u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

                        u_s = (u[j][i] + u[j - 1][i]) / 2;
                        v_s = (v[j][i] + v[j][i - 1]) / 2;
                        u_s_jp1 = (u[j + 1][i] + u[j][i]) / 2;
                        v_s_jp1 = (v[j + 1][i] + v[j + 1][i - 1]) / 2;

                        Hx = -((pow(u_cc, 2) - pow(u_cc_im1, 2)) / dx) - (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

                        step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * ((u[j][i + 1] - 2 * u[j][i] + u[j][i - 1]) / pow(dx, 2) + (u[j + 1][i] - 2 * u[j][i] + u[j - 1][i]) / pow(dy, 2)));
                    
                        Hx_old[j][i] = Hx;
                        
                    }
                    else if (j == Bj + void_len && i >= Bi && i <= Bi + void_len) { //top boundary of void (and top right corner point)
                        u_cc = (u[j][i] + u[j][i + 1]) / 2;
                        u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

                        u_s = 0;    
                        v_s = 0;
                        u_s_jp1 = (u[j + 1][i] + u[j][i]) / 2;
                        v_s_jp1 = (v[j + 1][i] + v[j + 1][i - 1]) / 2;

                        Hx = -((pow(u_cc, 2) - pow(u_cc_im1, 2)) / dx) - (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

                        step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * ((u[j][i + 1] - 2 * u[j][i] + u[j][i - 1]) / pow(dx, 2) +    ( (u[j+1][i] - u[j][i]) / dy - (u[j][i] - u_s) *(2/dy) ) / dy   ));
                        
                        Hx_old[j][i] = Hx;
                        
                    }
                    else if (j == Bj - 1 && i >= Bi && i <= Bi + void_len) { //bottom boundary of void (and bottom right corner point)
                        u_cc = (u[j][i] + u[j][i + 1]) / 2;
                        u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

                        u_s = (u[j][i] + u[j - 1][i]) / 2;
                        v_s = (v[j][i] + v[j][i - 1]) / 2;
                        u_s_jp1 = 0;
                        v_s_jp1 = 0;

                        Hx = -((pow(u_cc, 2) - pow(u_cc_im1, 2)) / dx) - (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

                        step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * ((u[j][i + 1] - 2 * u[j][i] + u[j][i - 1]) / pow(dx, 2) +         ( (u_s_jp1 - u[j][i]) *(2/dy) - (u[j][i] - u[j-1][i]) / dy) / dy        ));
                    
                        Hx_old[j][i] = Hx;
                        
                    }      
                    
                }
                else if (i > 0 && j > 0 && i < (Nx - 1) && j < (Ny - 1)) { // General interior points  
                    u_cc = (u[j][i] + u[j][i + 1]) / 2;
                    u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

                    u_s = (u[j][i] + u[j - 1][i]) / 2;
                    v_s = (v[j][i] + v[j][i - 1]) / 2;
                    u_s_jp1 = (u[j + 1][i] + u[j][i]) / 2;
                    v_s_jp1 = (v[j + 1][i] + v[j + 1][i - 1]) / 2;

                    Hx = -((pow(u_cc, 2) - pow(u_cc_im1, 2)) / dx) - (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

                    step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * ((u[j][i + 1] - 2 * u[j][i] + u[j][i - 1]) / pow(dx, 2) + (u[j + 1][i] - 2 * u[j][i] + u[j - 1][i]) / pow(dy, 2)));
                
                    Hx_old[j][i] = Hx;

                }
                else if (i == 0) { //left wall (and left corners)
                    
                    step1_mat_x[j][i] = 0;

                }
                else if (i == Nx - 1 && j != 0 && j != Ny - 1) { //right wall
                    u_cc = u[j][i];
                    u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

                    u_s = (u[j][i] + u[j - 1][i]) / 2;
                    v_s = (v[j][i] + v[j][i - 1]) / 2;
                    u_s_jp1 = (u[j + 1][i] + u[j][i]) / 2;
                    v_s_jp1 = (v[j + 1][i] + v[j + 1][i - 1]) / 2;

                    Hx = -((pow(u_cc, 2) - pow(u_cc_im1, 2)) / dx) - (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

                    step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * ((u[j][i] - 2 * u[j][i] + u[j][i - 1]) / pow(dx, 2) + (u[j + 1][i] - 2 * u[j][i] + u[j - 1][i]) / pow(dy, 2)));
                    
                    Hx_old[j][i] = Hx;
                    
                    
                }
                else if (j == Ny - 1 && i != 0 && i != Nx - 1) { //top wall
                    u_cc = (u[j][i] + u[j][i + 1]) / 2;
                    u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

                    u_s = (u[j][i] + u[j - 1][i]) / 2;
                    v_s = (v[j][i] + v[j][i - 1]) / 2;
                    u_s_jp1 = 0;
                    v_s_jp1 = 0;

                    Hx = -((pow(u_cc, 2) - pow(u_cc_im1, 2)) / dx) - (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

                    step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * ((u[j][i + 1] - 2 * u[j][i] + u[j][i - 1]) / pow(dx, 2) +         ( (u_s_jp1 - u[j][i]) *(2/dy) - (u[j][i] - u[j - 1][i]) / dy) / dy        ));
                
                    Hx_old[j][i] = Hx;
                    
                }
                else if (j == 0 && i != 0 && i != Nx - 1) { //bottom wall
                    u_cc = (u[j][i] + u[j][i + 1]) / 2;
                    u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

                    u_s = 0;    
                    v_s = 0;
                    u_s_jp1 = (u[j + 1][i] + u[j][i]) / 2;
                    v_s_jp1 = (v[j + 1][i] + v[j + 1][i - 1]) / 2;

                    Hx = -((pow(u_cc, 2) - pow(u_cc_im1, 2)) / dx) - (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

                    step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * ((u[j][i + 1] - 2 * u[j][i] + u[j][i - 1]) / pow(dx, 2) +    ( (u[j+1][i] - u[j][i]) / dy - (u[j][i] - u_s) *(2/dy) ) / dy   ));
                
                    Hx_old[j][i] = Hx;
                    
                }

                


                

                
                


                ////////// Calculating step1_mat_y (at all but corner points)  /////////
                if (i >= Bi - 1 && i <= Bi + void_len && j >= Bj - 1 && j <= Bj + void_len + 1) { //void and boundary points

                    if (i >= Bi && i < Bi + void_len && j >= Bj && j <= Bj + void_len) { //points that are inside or at a wall, do nothing
                        
                    }
                    else if (i == Bi - 1 && j > Bj - 1 && j <= Bj + void_len) {  //left boundary of void (and top left corner point)
                        v_cc = (v[j][i] + v[j + 1][i]) / 2;
                        v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

                        u_s = (u[j][i] + u[j - 1][i]) / 2;
                        v_s = (v[j][i] + v[j][i - 1]) / 2;
                        u_s_ip1 = 0;
                        v_s_ip1 = 0;

                        Hy = -((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / dy) - (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

                        step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) * (   ( (0 - v[j][i])*2/dx - (v[j][i] - v[j][i-1])/dx ) / dx        + (v[j + 1][i] - 2 * v[j][i] + v[j - 1][i]) / pow(dy, 2)));
                
                        Hy_old[j][i] = Hy;
                        
                    }
                    else if (i == Bi + void_len && j >= Bj && j <= Bj + void_len) { //right boundary (and top right corner point)
                        v_cc = (v[j][i] + v[j + 1][i]) / 2;
                        v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

                        u_s = 0;
                        v_s = 0;
                        u_s_ip1 = (u[j][i + 1] + u[j - 1][i + 1]) / 2;
                        v_s_ip1 = (v[j][i + 1] + v[j][i]) / 2;

                        Hy = -((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / dy) - (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

                        step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) *          (     ( (v[j][i+1] - v[j][i])/dx - (v[j][i] - v_s)*2/dx )  / dx         + (v[j + 1][i] - 2 * v[j][i] + v[j - 1][i]) / pow(dy, 2)));
                
                        Hy_old[j][i] = Hy;
                        
                    }
                    else if (j == Bj + void_len + 1 && i >= Bi - 1 && i <= Bi + void_len) { //top boundary (plus extra corner points)
                        v_cc = (v[j][i] + v[j + 1][i]) / 2;
                        v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

                        u_s = (u[j][i] + u[j - 1][i]) / 2;
                        v_s = (v[j][i] + v[j][i - 1]) / 2;
                        u_s_ip1 = (u[j][i + 1] + u[j - 1][i + 1]) / 2;
                        v_s_ip1 = (v[j][i + 1] + v[j][i]) / 2;

                        Hy = -((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / dy) - (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

                        step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) * ((v[j][i + 1] - 2 * v[j][i] + v[j][i - 1]) / pow(dx, 2) + (v[j + 1][i] - 2 * v[j][i] + v[j - 1][i]) / pow(dy, 2)));
                    
                        Hy_old[j][i] = Hy;
                        
                    }
                    else if (j == Bj - 1 && i >= Bi - 1 && i <= Bi + void_len) { //bottom boundary (and bottom left and bottom right points)
                        v_cc = (v[j][i] + v[j + 1][i]) / 2;
                        v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

                        u_s = (u[j][i] + u[j - 1][i]) / 2;
                        v_s = (v[j][i] + v[j][i - 1]) / 2;
                        u_s_ip1 = (u[j][i + 1] + u[j - 1][i + 1]) / 2;
                        v_s_ip1 = (v[j][i + 1] + v[j][i]) / 2;

                        Hy = -((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / dy) - (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

                        step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) * ((v[j][i + 1] - 2 * v[j][i] + v[j][i - 1]) / pow(dx, 2) +     (v[j + 1][i] - 2 * v[j][i] + v[j - 1][i]) / pow(dy, 2)));
                
                        Hy_old[j][i] = Hy;
                        
                    } 

                    
                }
                else if (i > 0 && j > 0 && i < (Nx - 1) && j < (Ny - 1)) { // General interior points  
                    v_cc = (v[j][i] + v[j + 1][i]) / 2;
                    v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

                    u_s = (u[j][i] + u[j - 1][i]) / 2;
                    v_s = (v[j][i] + v[j][i - 1]) / 2;
                    u_s_ip1 = (u[j][i + 1] + u[j - 1][i + 1]) / 2;
                    v_s_ip1 = (v[j][i + 1] + v[j][i]) / 2;

                    Hy = -((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / dy) - (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

                    step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) * ((v[j][i + 1] - 2 * v[j][i] + v[j][i - 1]) / pow(dx, 2) + (v[j + 1][i] - 2 * v[j][i] + v[j - 1][i]) / pow(dy, 2)));
                
                    Hy_old[j][i] = Hy;

                }

                else if (i == 0 && j != 0 && j != Ny - 1) { //left wall
                    
                    v_cc = (v[j][i] + v[j + 1][i]) / 2;
                    v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

                    u_s = U_inlet;
                    v_s = 0;
                    u_s_ip1 = (u[j][i + 1] + u[j - 1][i + 1]) / 2;
                    v_s_ip1 = (v[j][i + 1] + v[j][i]) / 2;

                    Hy = -((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / dy) - (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

                    step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) *          (     ( (v[j][i+1] - v[j][i])/dx - (v[j][i] - v_s)*2/dx )  / dx         + (v[j + 1][i] - 2 * v[j][i] + v[j - 1][i]) / pow(dy, 2)));
                
                    Hy_old[j][i] = Hy;
                    
                }

                else if (i == Nx - 1 && j != 0 && j != Ny - 1) { //right wall
                    v_cc = (v[j][i] + v[j + 1][i]) / 2;
                    v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

                    u_s = (u[j][i] + u[j - 1][i]) / 2;
                    v_s = (v[j][i] + v[j][i - 1]) / 2;
                    u_s_ip1 =  u_s;     
                    v_s_ip1 = v[j][i];          // would have been more natural to use v_s but that led to spurious oscillations at this edge of the domain

                    Hy = -((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / dy) - (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

                    step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) * (     (v[j][i] - 2 * v[j][i] + v[j][i - 1]) / pow(dx, 2)        + (v[j + 1][i] - 2 * v[j][i] + v[j - 1][i]) / pow(dy, 2)));
                
                    Hy_old[j][i] = Hy;
                    
                }

                else if (j == Ny - 1 && i != 0 && i != Nx - 1) { //top wall
                    v_cc = (v[j][i] + 0) / 2;
                    v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

                    u_s = (u[j][i] + u[j - 1][i]) / 2;
                    v_s = (v[j][i] + v[j][i - 1]) / 2;
                    u_s_ip1 = (u[j][i + 1] + u[j - 1][i + 1]) / 2;
                    v_s_ip1 = (v[j][i + 1] + v[j][i]) / 2;

                    Hy = -((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / dy) - (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

                    step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) * ((v[j][i + 1] - 2 * v[j][i] + v[j][i - 1]) / pow(dx, 2) +     (0 - 2 * v[j][i] + v[j - 1][i]) / pow(dy, 2)));
                
                    Hy_old[j][i] = Hy;
                    
                }

                else if (j == 0) { //bottom wall (and bottom corners)

                    step1_mat_y[j][i] = 0;
                    
                    
                }



                


            }
        }


        //////// Calculating step1_mat_x and step1_mat_y at non-fixed corner points /////////


        //// top right corner 
        i = Nx - 1;
        j = Ny - 1;

        // x data
        u_cc = u[j][i];
        u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

        u_s = (u[j][i] + u[j - 1][i]) / 2;
        v_s = (v[j][i] + v[j][i - 1]) / 2;
        u_s_jp1 = 0;
        v_s_jp1 = 0;

        Hx = -((pow(u_cc, 2) - pow(u_cc_im1, 2)) / dx) - (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

        step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * (       (u[j][i] - 2 * u[j][i] + u[j][i - 1]) / pow(dx, 2)   +     ( (u_s_jp1 - u[j][i]) *(2/dy) - (u[j][i] - u[j-1][i]) / dy) / dy      ));

        Hx_old[j][i] = Hx;

        // y data
        v_cc = (v[j][i] + 0) / 2;
        v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

        u_s = (u[j][i] + u[j - 1][i]) / 2;
        v_s = (v[j][i] + v[j][i - 1]) / 2;
        u_s_ip1 = u_s;
        v_s_ip1 = v_s;

        Hy = -((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / dy) - (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

        step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) * (      (v[j][i] - 2 * v[j][i] + v[j][i - 1]) / pow(dx, 2)       +     (0 - 2 * v[j][i] + v[j - 1][i]) / pow(dy, 2)));

        Hy_old[j][i] = Hy;




        //// bottom right corner
        i = Nx - 1;
        j = 0;

        //x data
        u_cc = u[j][i];
        u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

        u_s = 0;
        v_s = 0;
        u_s_jp1 = (u[j + 1][i] + u[j][i]) / 2;
        v_s_jp1 = (v[j + 1][i] + v[j + 1][i - 1]) / 2;

        Hx = -((pow(u_cc, 2) - pow(u_cc_im1, 2)) / dx) - (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

        step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * (   (u[j][i] - 2 * u[j][i] + u[j][i - 1]) / pow(dx, 2)    +       ( (u[j + 1][i] - u[j][i]) / dy - (u[j][i] - u_s) *(2/dy) ) / dy   ));

        Hx_old[j][i] = Hx;




        //// top left corner
        i = 0;
        j = Ny - 1;

        // y data
        v_cc = (v[j][i] + 0) / 2;
        v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

        u_s = U_inlet;
        v_s = 0;
        u_s_ip1 = (u[j][i + 1] + u[j - 1][i + 1]) / 2;
        v_s_ip1 = (v[j][i + 1] + v[j][i]) / 2;

        Hy = -((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / dy) - (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

        step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) * (    ( (v[j][i+1] - v[j][i])/dx - (v[j][i] - v_s)*2/dx )  / dx +     (0 - 2 * v[j][i] + v[j - 1][i]) / pow(dy, 2)));

        Hy_old[j][i] = Hy;





        /* ------------- Tridiagonal solve for du_ss, dv_ss, du_s and dv_s ------------ */


        ///////// Tridiagonal solve for du_ss

        // rows below the bottom of the void
        for (j = 0; j < Bj; j++) {

            //fill in temp_vec with row of pixels except first pixel
            for(i = 0; i < Nx - 1; i++) {
                temp_vec_x_long[i] = step1_mat_x[j][i + 1];
            }

            //gaussian elimination
            TriDiag_GaussElim(Nx, dx, dt, Re, temp_vec_x_long, 1, 0, 0);

            //update du_ss
            for(i = 0; i < Nx - 1; i++) {
                du_ss[j][i + 1] = temp_vec_x_long[i];
            }

        }

        // rows next to the void
        for (j = Bj; j < Bj + void_len; j++) {

            // left of the void
            for(i = 0; i < Bi - 1; i++) {
                temp_vec_x_small[i] = step1_mat_x[j][i + 1];
            }
            
            TriDiag_GaussElim(Bi, dx, dt, Re, temp_vec_x_small, 0, 0, 0);

            for(i = 0; i < Bi - 1; i++) {
                du_ss[j][i + 1] = temp_vec_x_small[i];
            }



            // right of the void
            for(i = 0; i < Nx - Bi - void_len - 1; i++) {
                temp_vec_x_medium[i] = step1_mat_x[j][i + Bi + void_len + 1];
            }

            TriDiag_GaussElim(Nx - Bi - void_len, dx, dt, Re, temp_vec_x_medium, 1, 0, 0);

            for(i = 0; i < Nx - Bi - void_len - 1; i++) {
                du_ss[j][i + Bi + void_len + 1] = temp_vec_x_medium[i];
            }


        }



        // rows above the top of the void
        for (j = Bj + void_len; j < Ny; j++) {

            for(i = 0; i < Nx - 1; i++) {
                temp_vec_x_long[i] = step1_mat_x[j][i + 1];
            }

            TriDiag_GaussElim(Nx, dx, dt, Re, temp_vec_x_long, 1, 0, 0);

            for(i = 0; i < Nx - 1; i++) {
                du_ss[j][i + 1] = temp_vec_x_long[i];
            }

        }




        ///////// Tridiagonal solve for dv_ss

        // rows below the void
        for (j = 1; j < Bj; j++) {

            // fill in temp_vec with row of pixels
            for(i = 0; i < Nx; i++) {
                temp_vec_x_long2[i] = step1_mat_y[j][i];
            }

            // gaussian elimination
            TriDiag_GaussElim(Nx + 1, dx, dt, Re, temp_vec_x_long2, 1, 0, 1);

            // update dv_ss
            for(i = 0; i < Nx; i++) {
                dv_ss[j][i] = temp_vec_x_long2[i];
            }

        }


        // rows next to the void
        for (j = Bj; j <= Bj + void_len; j++) {

            // left
            for(i = 0; i < Bi; i++) {
                temp_vec_x_small2[i] = step1_mat_y[j][i];
            }
            
            TriDiag_GaussElim(Bi + 1, dx, dt, Re, temp_vec_x_small2, 0, 1, 0);

            for(i = 0; i < Bi; i++) {
                dv_ss[j][i] = temp_vec_x_small2[i];
            }


            // right
            for(i = 0; i < Nx - Bi - void_len; i++) {
                temp_vec_x_medium2[i] = step1_mat_y[j][i + Bi + void_len];
            }

            TriDiag_GaussElim(Nx - Bi - void_len + 1, dx, dt, Re, temp_vec_x_medium2, 1, 0, 1);

            for(i = 0; i < Nx - Bi - void_len; i++) {
                dv_ss[j][i + Bi + void_len] = temp_vec_x_medium2[i];
            }

        }


        // rows above the void 
        for (j = Bj + void_len + 1; j < Ny; j++) {

            for(i = 0; i < Nx; i++) {
                temp_vec_x_long2[i] = step1_mat_y[j][i];
            }

            TriDiag_GaussElim(Nx + 1, dx, dt, Re, temp_vec_x_long2, 1, 0, 1);

            for(i = 0; i < Nx; i++) {
                dv_ss[j][i] = temp_vec_x_long2[i];
            }

        }


        ///////// Tridiagonal solve for du_s

        // columns to the left of the square
        for (i = 1; i < Bi; i++) {

            for (j = 0; j < Ny; j++) {
                temp_vec_y_long[j] = du_ss[j][i];
            }

            TriDiag_GaussElim(Ny + 1, dy, dt, Re, temp_vec_y_long, 0, 1, 0);

            for (j = 0; j < Ny; j++) {
                du_s[j][i] = temp_vec_y_long[j];
            }

        }



        // columns above and below the void
        for (i = Bi; i < Bi + void_len + 1; i++) {

            // below
            for (j = 0; j < Bj; j++) {
                temp_vec_y_small[j] = du_ss[j][i];
            }

            TriDiag_GaussElim(Bj + 1, dy, dt, Re, temp_vec_y_small, 0, 1, 0);

            for (j = 0; j < Bj; j++) {
                du_s[j][i] = temp_vec_y_small[j];
            }


            //above
            for (j = 0; j < Ny - Bj - void_len; j++) {
                temp_vec_y_small3[j] = du_ss[j + Bj + void_len][i];
            }

            TriDiag_GaussElim(Ny - Bj - void_len + 1, dy, dt, Re, temp_vec_y_small3, 0, 1, 0);

            for (j = 0; j < Ny - Bj - void_len; j++) {
                du_s[j + Bj + void_len][i] = temp_vec_y_small3[j];
            }


        }




        // columns to the right of the void
        for (i = Bi + void_len + 1; i < Nx; i++) {

            for (j = 0; j < Ny; j++) {
                temp_vec_y_long[j] = du_ss[j][i];
            }

            TriDiag_GaussElim(Ny + 1, dy, dt, Re, temp_vec_y_long, 0, 1, 0);

            for (j = 0; j < Ny; j++) {
                du_s[j][i] = temp_vec_y_long[j];
            }

        }





        ///////// Tridiagonal solve for dv_s

        // columns to the left of the void
        for (i = 0; i < Bi; i++) {

            for (j = 0; j < Ny - 1; j++) {
                temp_vec_y_long2[j] = dv_ss[j + 1][i];
            }

            TriDiag_GaussElim(Ny, dy, dt, Re, temp_vec_y_long2, 0, 0, 0);

            for (j = 0; j < Ny - 1; j++) {
                dv_s[j + 1][i] = temp_vec_y_long2[j];
            }

        }


        // columns above and below the void
        for (i = Bi; i < Bi + void_len; i++) {

            // below
            for (j = 0; j < Bj - 1; j++) {
                temp_vec_y_small2[j] = dv_ss[j + 1][i];
            }

            TriDiag_GaussElim(Bj, dy, dt, Re, temp_vec_y_small2, 0, 0, 0);

            for (j = 0; j < Bj - 1; j++) {
                dv_s[j + 1][i] = temp_vec_y_small2[j];
            }



            // above
            for (j = 0; j < Ny - Bj - void_len - 1; j++) {
                temp_vec_y_small4[j] = dv_ss[j + Bj + void_len + 1][i];
            }

            TriDiag_GaussElim(Ny - Bj - void_len, dy, dt, Re, temp_vec_y_small4, 0, 0, 0);

            for (j = 0; j < Ny - Bj - void_len - 1; j++) {
                dv_s[j + Bj + void_len + 1][i] = temp_vec_y_small4[j];
            }

        }



        // columns to the right of the void
        for (i = Bi + void_len; i < Nx; i++) {

            for (j = 0; j < Ny - 1; j++) {
                temp_vec_y_long2[j] = dv_ss[j + 1][i];
            }

            TriDiag_GaussElim(Ny, dy, dt, Re, temp_vec_y_long2, 0, 0, 0);

            for (j = 0; j < Ny - 1; j++) {
                dv_s[j + 1][i] = temp_vec_y_long2[j];
            }

        }




        
        


        

                
        /* ------------- Update the velocities u_star and v_star ------------ */

        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
                u_star[j][i] = u[j][i] + du_s[j][i];
                v_star[j][i] = v[j][i] + dv_s[j][i];
            }
        }
    

        // Update right boundary of u_star (to be used in step 3)
        for (j = 0; j < Ny; j++) {
            u_star[j][Nx] = u_star[j][Nx - 1];
        }


        




        /* -------------------------------------------------------------------------
        ----------------------------------- Step 2 ---------------------------------
        ------------------------------------------------------------------------- */

        // Calculating the cell-centered divergence of velocity_star divided by dt
        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
            
                grad_u_star_over_dt[j][i] = ((u_star[j][i + 1] - u_star[j][i]) / dx + (v_star[j + 1][i] - v_star[j][i]) / dy) / dt;


            }
        }


        // Solving for cell-centered pressure using multigrid acceleration method 
        
        GS_nstep(grad_u_star_over_dt, p, Nx, Ny, dx, dy, epsilon, nGS, void_len, Bj, Bi);



        





        /* -------------------------------------------------------------------------
        ----------------------------------- Step 3 ---------------------------------
        ------------------------------------------------------------------------- */




        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {

                if ( j >= Bj && j <= Bj + void_len - 1 && i >= Bi && i <= Bi + void_len - 1) { // inside void
                    
                    u[j][i] = 0;
                    v[j][i] = 0;
                }
                else if (i == 0 && j != 0) { // left edge of domain (and top left corner)
    
                    u[j][i] = U_inlet;
                    v[j][i] = v_star[j][i] - dt * (p[j][i] - p[j - 1][i]) / dy;
                }
                else if (i == Nx - 1 && j !=0) { // right edge of domain (and top right corner)
                    
                    u[j][i] = u_star[j][i] - dt * (p[j][i] - p[j][i - 1]) / dx;
                    v[j][i] = v_star[j][i] - dt * (p[j][i] - p[j - 1][i]) / dy;
                }
                else if (j == 0 && i !=0) { // bottom edge of domain (and bottom right corner)
                    
                    u[j][i] = u_star[j][i] - dt * (p[j][i] - p[j][i - 1]) / dx;
                    v[j][i] = 0;
                }
                else if (j == Ny - 1 && i !=0 && i != Nx - 1) { // top edge of domain
                
                    u[j][i] = u_star[j][i] - dt * (p[j][i] - p[j][i - 1]) / dx;
                    v[j][i] = v_star[j][i] - dt * (p[j][i] - p[j - 1][i]) / dy;
                }
                else if (i == 0 && j == 0) { // bottom left corner
                    u[j][i] = U_inlet;
                    v[j][i] = 0;
                }
                else if (i == Bi - 1 && j >= Bj && j <= Bj + void_len - 1) { // left edge of void
                    
                    u[j][i] = u_star[j][i] - dt * (p[j][i] - p[j][i - 1]) / dx;
                    v[j][i] = v_star[j][i] - dt * (p[j][i] - p[j - 1][i]) / dy;
                }
                else if (i == Bi + void_len && j >= Bj && j <= Bj + void_len - 1) { // right edge of void
                   
                    u[j][i] = 0;
                    v[j][i] = v_star[j][i] - dt * (p[j][i] - p[j - 1][i]) / dy;
                }
                else if (j == Bj + void_len && i >= Bi && i <= Bi + void_len - 1) { // top edge of void
                    
                    u[j][i] = u_star[j][i] - dt * (p[j][i] - p[j][i - 1]) / dx;
                    v[j][i] = 0;
                }
                else if (j == Bj - 1 && i >= Bi && i <= Bi + void_len - 1) { // bottom edge of void
                    
                    u[j][i] = u_star[j][i] - dt * (p[j][i] - p[j][i - 1]) / dx;
                    v[j][i] = v_star[j][i] - dt * (p[j][i] - p[j - 1][i]) / dy;
                }
                else { // normal space not next to a bounary
                    
                    u[j][i] = u_star[j][i] - dt * (p[j][i] - p[j][i - 1]) / dx;
                    v[j][i] = v_star[j][i] - dt * (p[j][i] - p[j - 1][i]) / dy;

                }

            }
        }











        

        /* -------------------------------------------------------------------------
        ----------------------------- Scalar Transport -----------------------------
        ------------------------------------------------------------------------- */

        double u_phi_left;
        double u_phi_right;
        double v_phi_top;
        double v_phi_bottom;


        if ( strcmp(convective_method, "upwind") == 0 ) {
            if (iter % 10 == 1) { 
                printf("Step %d with upwind for scalar transport \n", iter);
            }

            //interior points

            for (j = 1; j < Ny - 1; j++) {
                for (i = 1; i < Nx - 1; i++) {

                    diffu = Diff * (   (phi[j][i+1] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - 2*phi[j][i] + phi[j-1][i])/(pow(dy,2))    );

                    if (u[j][i] > 0) { u_phi_left = u[j][i] * phi[j][i - 1]; }
                    else {             u_phi_left = u[j][i] * phi[j][i];     }

                    if (v[j][i] > 0) { v_phi_bottom = v[j][i] * phi[j - 1][i]; }
                    else {             v_phi_bottom = v[j][i] * phi[j][i];     }

                    if (u[j][i + 1] > 0) { u_phi_right = u[j][i + 1] * phi[j][i];     }
                    else {                 u_phi_right = u[j][i + 1] * phi[j][i + 1]; }

                    if (v[j + 1][i] > 0) { v_phi_top = v[j + 1][i] * phi[j][i]; }
                    else {                 v_phi_top = v[j + 1][i] * phi[j + 1][i]; }

                    convec = (u_phi_right - u_phi_left) / dx + (v_phi_top - v_phi_bottom) / dy;

                    phi_new[j][i] = phi[j][i] + dt * (diffu - convec);

                }
            }


            //top boundary points except corners
            j = Ny - 1;

            for (i = 1; i < Nx - 1; i++) {

                diffu = Diff * (   (phi[j][i+1] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   ( - phi[j][i] + phi[j-1][i])/(pow(dy,2))    );

                if (u[j][i] > 0) { u_phi_left = u[j][i] * phi[j][i - 1]; }
                else {             u_phi_left = u[j][i] * phi[j][i];     }

                if (v[j][i] > 0) { v_phi_bottom = v[j][i] * phi[j - 1][i]; }
                else {             v_phi_bottom = v[j][i] * phi[j][i];     }

                if (u[j][i + 1] > 0) { u_phi_right = u[j][i + 1] * phi[j][i];     }
                else {                 u_phi_right = u[j][i + 1] * phi[j][i + 1]; }

                v_phi_top = 0;

                convec = (u_phi_right - u_phi_left) / dx + (v_phi_top - v_phi_bottom) / dy;

                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);

            }


            //right boundary points except corners
            i = Nx - 1;

            for (j = 1; j < Ny - 1; j++) {

                diffu = Diff * (   (phi[j][i] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - 2*phi[j][i] + phi[j-1][i])/(pow(dy,2))    );

                if (u[j][i] > 0) { u_phi_left = u[j][i] * phi[j][i - 1]; }
                else {             u_phi_left = u[j][i] * phi[j][i];     }

                if (v[j][i] > 0) { v_phi_bottom = v[j][i] * phi[j - 1][i]; }
                else {             v_phi_bottom = v[j][i] * phi[j][i];     }

                u_phi_right = u[j][i] * phi[j][i]; // because of neumann BC
                
                if (v[j + 1][i] > 0) { v_phi_top = v[j + 1][i] * phi[j][i]; }
                else {                 v_phi_top = v[j + 1][i] * phi[j + 1][i]; }

                convec = (u_phi_right - u_phi_left) / dx + (v_phi_top - v_phi_bottom) / dy;

                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }


            //left boundary points
            i = 0;

            for (j = 0; j < Ny; j++) {
                phi_new[j][i] = exp( -1 * pow(j * dy - 2.5 * H + 0.5 * dy ,2) );
            }


            //bottom boundary points (except at corners)
            j = 0;

            for (i = 1; i < Nx - 1; i++) {

                diffu = Diff * (   (phi[j][i+1] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - phi[j][i] )/(pow(dy,2))    );

                if (u[j][i] > 0) { u_phi_left = u[j][i] * phi[j][i - 1]; }
                else {             u_phi_left = u[j][i] * phi[j][i];     }

                v_phi_bottom = 0;
                
                if (u[j][i + 1] > 0) { u_phi_right = u[j][i + 1] * phi[j][i];     }
                else {                 u_phi_right = u[j][i + 1] * phi[j][i + 1]; }

                if (v[j + 1][i] > 0) { v_phi_top = v[j + 1][i] * phi[j][i]; }
                else {                 v_phi_top = v[j + 1][i] * phi[j + 1][i]; }

                convec = (u_phi_right - u_phi_left) / dx + (v_phi_top - v_phi_bottom) / dy;

                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }
   

            //top right corner

            j = Ny - 1;
            i = Nx - 1;

            diffu = Diff * (   (- phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (- phi[j][i] + phi[j-1][i])/(pow(dy,2))    );
            
            if (u[j][i] > 0) { u_phi_left = u[j][i] * phi[j][i - 1]; }
            else {             u_phi_left = u[j][i] * phi[j][i];     }

            if (v[j][i] > 0) { v_phi_bottom = v[j][i] * phi[j - 1][i]; }
            else {             v_phi_bottom = v[j][i] * phi[j][i];     }

            u_phi_right = u[j][i] * phi[j][i]; // because of neumann BC

            v_phi_top = 0;

            convec = (u_phi_right - u_phi_left) / dx + (v_phi_top - v_phi_bottom) / dy;

            phi_new[j][i] = phi[j][i] + dt * (diffu - convec);


            //bottom right corner

            j = 0;
            i = Nx - 1;

            diffu = Diff * (   ( - phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - phi[j][i] )/(pow(dy,2))    );
            
            if (u[j][i] > 0) { u_phi_left = u[j][i] * phi[j][i - 1]; }
            else {             u_phi_left = u[j][i] * phi[j][i];     }

            v_phi_bottom = 0;

            u_phi_right = u[j][i] * phi[j][i]; // because of neumann BC

            if (v[j + 1][i] > 0) { v_phi_top = v[j + 1][i] * phi[j][i]; }
            else {                 v_phi_top = v[j + 1][i] * phi[j + 1][i]; }

            convec = (u_phi_right - u_phi_left) / dx + (v_phi_top - v_phi_bottom) / dy;

            phi_new[j][i] = phi[j][i] + dt * (diffu - convec);



            

            //bottom boundary of void
            j = Bj - 1;
            for (i = Bi; i < Bi + void_len; i++) {
                diffu = Diff * (   (phi[j][i+1] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   ( - phi[j][i] + phi[j-1][i])/(pow(dy,2))    );

                if (u[j][i] > 0) { u_phi_left = u[j][i] * phi[j][i - 1]; }
                else {             u_phi_left = u[j][i] * phi[j][i];     }

                if (v[j][i] > 0) { v_phi_bottom = v[j][i] * phi[j - 1][i]; }
                else {             v_phi_bottom = v[j][i] * phi[j][i];     }

                if (u[j][i + 1] > 0) { u_phi_right = u[j][i + 1] * phi[j][i];     }
                else {                 u_phi_right = u[j][i + 1] * phi[j][i + 1]; }

                v_phi_top = 0;
                convec = (u_phi_right - u_phi_left) / dx + (v_phi_top - v_phi_bottom) / dy;
            
                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }

            //top boundary of void
            j = Bj + void_len;
            for (i = Bi; i < Bi + void_len; i++) {
                diffu = Diff * (   (phi[j][i+1] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - phi[j][i] )/(pow(dy,2))    );
                
                if (u[j][i] > 0) { u_phi_left = u[j][i] * phi[j][i - 1]; }
                else {             u_phi_left = u[j][i] * phi[j][i];     }

                v_phi_bottom = 0;

                if (u[j][i + 1] > 0) { u_phi_right = u[j][i + 1] * phi[j][i];     }
                else {                 u_phi_right = u[j][i + 1] * phi[j][i + 1]; }

                if (v[j + 1][i] > 0) { v_phi_top = v[j + 1][i] * phi[j][i]; }
                else {                 v_phi_top = v[j + 1][i] * phi[j + 1][i]; }
                convec = (u_phi_right - u_phi_left) / dx + (v_phi_top - v_phi_bottom) / dy;
            
                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }

           

            //left boundary of void
            i = Bi - 1;
            for (j = Bj; j < Bj + void_len; j++) {
                diffu = Diff * (   ( - phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - 2*phi[j][i] + phi[j-1][i])/(pow(dy,2))    );
                
                if (u[j][i] > 0) { u_phi_left = u[j][i] * phi[j][i - 1]; }
                else {             u_phi_left = u[j][i] * phi[j][i];     }

                if (v[j][i] > 0) { v_phi_bottom = v[j][i] * phi[j - 1][i]; }
                else {             v_phi_bottom = v[j][i] * phi[j][i];     }

                u_phi_right = 0;

                if (v[j + 1][i] > 0) { v_phi_top = v[j + 1][i] * phi[j][i]; }
                else {                 v_phi_top = v[j + 1][i] * phi[j + 1][i]; }
                convec = (u_phi_right - u_phi_left) / dx + (v_phi_top - v_phi_bottom) / dy;
            
                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }

           

            //right boundary of void
            i = Bi + void_len;
            for (j = Bj; j < Bj + void_len; j++) {
                diffu = Diff * (   (phi[j][i+1] - phi[j][i] ) / (pow(dx,2))   +   (phi[j + 1][i] - 2*phi[j][i] + phi[j-1][i])/(pow(dy,2))    );
                
                u_phi_left = 0;

                if (v[j][i] > 0) { v_phi_bottom = v[j][i] * phi[j - 1][i]; }
                else {             v_phi_bottom = v[j][i] * phi[j][i];     }

                if (u[j][i + 1] > 0) { u_phi_right = u[j][i + 1] * phi[j][i];     }
                else {                 u_phi_right = u[j][i + 1] * phi[j][i + 1]; }

                if (v[j + 1][i] > 0) { v_phi_top = v[j + 1][i] * phi[j][i]; }
                else {                 v_phi_top = v[j + 1][i] * phi[j + 1][i]; }
                convec = (u_phi_right - u_phi_left) / dx + (v_phi_top - v_phi_bottom) / dy;
            
                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }

           
            //inside square, set phi_new = 0

            for (j = Bj; j < Bj + void_len; j++) {                        
                for (i = Bi; i < Bi + void_len; i++) {
                    phi_new[j][i] = 0;
                }
            }


            // update phi everywhere

            for (j = 0; j < Ny; j++) {
                for (i = 0; i < Nx; i++) {
                    phi[j][i] = phi_new[j][i];
                }
            }



        }





        if ( strcmp(convective_method, "centraldiff") == 0 ) {
            if (iter % 10 == 1) { 
                printf("Step %d with central diff for scalar transport \n", iter);
            }
            
            //interior points
            for (j = 1; j < Ny - 1; j++) {
                for (i = 1; i < Nx - 1; i++) {

                    diffu = Diff * (   (phi[j][i+1] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - 2*phi[j][i] + phi[j-1][i])/(pow(dy,2))    );

                    convec =  ( u[j][i+1] * (phi[j][i] + phi[j][i+1]) / 2 - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (v[j + 1][i] * (phi[j][i] + phi[j + 1][i])/2 - v[j][i] * (phi[j - 1][i] + phi[j][i])/2 ) / dy;

                    phi_new[j][i] = phi[j][i] + dt * (diffu - convec);

                }
            }

           
            //top boundary points except corners
            j = Ny - 1;

            for (i = 1; i < Nx - 1; i++) {

                diffu = Diff * (   (phi[j][i+1] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   ( - phi[j][i] + phi[j-1][i])/(pow(dy,2))    );

                convec =  ( u[j][i+1] * (phi[j][i] + phi[j][i+1]) / 2 - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (0 - v[j][i] * (phi[j - 1][i] + phi[j][i])/2 ) / dy;

                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }


            //right boundary points except corners
            i = Nx - 1;

            for (j = 1; j < Ny - 1; j++) {

                diffu = Diff * (   (phi[j][i] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - 2*phi[j][i] + phi[j-1][i])/(pow(dy,2))    );

                convec =  ( u[j][i] * (phi[j][i] + phi[j][i]) / 2 - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (v[j + 1][i] * (phi[j][i] + phi[j + 1][i])/2 - v[j][i] * (phi[j - 1][i] + phi[j][i])/2 ) / dy;

                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }


            //left boundary points
            i = 0;

            for (j = 0; j < Ny; j++) {
                phi_new[j][i] = exp( -1 * pow(j * dy - 2.5 * H + 0.5 * dy ,2) );
            }


            //bottom boundary points except corners
            j = 0;

            for (i = 1; i < Nx - 1; i++) {

                diffu = Diff * (   (phi[j][i+1] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - phi[j][i] )/(pow(dy,2))    );

                convec =  ( u[j][i+1] * (phi[j][i] + phi[j][i+1]) / 2 - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (v[j + 1][i] * (phi[j][i] + phi[j + 1][i])/2 - 0 ) / dy;

                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }
   

            //top right corner

            j = Ny - 1;
            i = Nx - 1;

            diffu = Diff * (   (- phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (- phi[j][i] + phi[j-1][i])/(pow(dy,2))    );
            convec =  ( u[j][i] * (phi[j][i] + phi[j][i]) / 2 - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (  0  - v[j][i] * (phi[j - 1][i] + phi[j][i])/2 ) / dy;
            phi_new[j][i] = phi[j][i] + dt * (diffu - convec);


            //bottom right corner

            j = 0;
            i = Nx - 1;

            diffu = Diff * (   ( - phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - phi[j][i] )/(pow(dy,2))    );
            convec =  ( u[j][i] * (phi[j][i] + phi[j][i]) / 2 - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (v[j + 1][i] * (phi[j][i] + phi[j + 1][i])/2 - 0 ) / dy;
            phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
           


            

            //bottom boundary of void
            j = Bj - 1;
            for (i = Bi; i < Bi + void_len; i++) {
                diffu = Diff * (   (phi[j][i+1] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   ( - phi[j][i] + phi[j-1][i])/(pow(dy,2))    );
                convec =  ( u[j][i+1] * (phi[j][i] + phi[j][i+1]) / 2 - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (0 - v[j][i] * (phi[j - 1][i] + phi[j][i])/2 ) / dy;
                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }

            //top boundary of void
            j = Bj + void_len;
            for (i = Bi; i < Bi + void_len; i++) {
                diffu = Diff * (   (phi[j][i+1] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - phi[j][i] )/(pow(dy,2))    );
                convec =  ( u[j][i+1] * (phi[j][i] + phi[j][i+1]) / 2 - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (v[j + 1][i] * (phi[j][i] + phi[j + 1][i])/2 - 0 ) / dy;
                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }
           
            //left boundary of void
            i = Bi - 1;
            for (j = Bj; j < Bj + void_len; j++) {
                diffu = Diff * (   (0 - phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - 2*phi[j][i] + phi[j-1][i])/(pow(dy,2))    );
                convec =  (   0  - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (v[j + 1][i] * (phi[j][i] + phi[j + 1][i])/2 - v[j][i] * (phi[j - 1][i] + phi[j][i])/2 ) / dy;
                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }

            //right boundary of void
            i = Bi + void_len;
            for (j = Bj; j < Bj + void_len; j++) {
                diffu = Diff * (   (phi[j][i+1] - phi[j][i] ) / (pow(dx,2))   +   (phi[j + 1][i] - 2*phi[j][i] + phi[j-1][i])/(pow(dy,2))    );
                convec =  ( u[j][i+1] * (phi[j][i] + phi[j][i+1]) / 2 - 0 ) / dx   +   (v[j + 1][i] * (phi[j][i] + phi[j + 1][i])/2 - v[j][i] * (phi[j - 1][i] + phi[j][i])/2 ) / dy;
                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }


            //inside square, set phi_new = 0

            for (j = Bj; j < Bj + void_len; j++) {                        
                for (i = Bi; i < Bi + void_len; i++) {
                    phi_new[j][i] = 0;
                }
            }



            // update phi everywhere

            for (j = 0; j < Ny; j++) {
                for (i = 0; i < Nx; i++) {
                    phi[j][i] = phi_new[j][i];
                }
            }



        }





        if ( strcmp(convective_method, "quick") == 0 ) {
            printf("Step %d with quick for scalar transport \n", iter);
            
            /*
            if (iter % 10 == 1) { 
                printf("Step %d with quick for scalar transport \n", iter);
            }
            */

            //interior points

            for (j = 2; j < Ny - 2; j++) {
                for (i = 2; i < Nx - 2; i++) {

                    diffu = Diff * (   (phi[j][i+1] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - 2*phi[j][i] + phi[j-1][i])/(pow(dy,2))    );

                    if (u[j][i] > 0) { u_phi_left = u[j][i] * ( phi[j][i - 2]*(-1/6) + phi[j][i - 1]*5/6 + phi[j][i]*2/6 ); }
                    else {             u_phi_left = u[j][i] * ( phi[j][i - 1]*2/6 + phi[j][i]*5/6 + phi[j][i + 1]*(-1/6) ); }

                    if (v[j][i] > 0) { v_phi_bottom = v[j][i] * ( phi[j - 2][i]*(-1/6) + phi[j - 1][i]*5/6 + phi[j][i]*2/6 ); }
                    else {             v_phi_bottom = v[j][i] * ( phi[j - 1][i]*2/6 + phi[j][i]*5/6 + phi[j + 1][i]*(-1/6) ); }

                    if (u[j][i + 1] > 0) { u_phi_right = u[j][i+1] * ( phi[j][i - 1]*(-1/6) + phi[j][i]*5/6 + phi[j][i + 1]*2/6 ); }
                    else {                 u_phi_right = u[j][i+1] * ( phi[j][i]*2/6 + phi[j][i + 1]*5/6 + phi[j][i + 2]*(-1/6) ); }

                    if (v[j + 1][i] > 0) { v_phi_top = v[j+1][i] * ( phi[j - 1][i]*(-1/6) + phi[j][i]*5/6 + phi[j + 1][i]*2/6 ); }
                    else {                 v_phi_top = v[j+1][i] * ( phi[j][i]*2/6 + phi[j + 1][i]*5/6 + phi[j + 2][i]*(-1/6) ); }

                    convec = (u_phi_right - u_phi_left) / dx + (v_phi_top - v_phi_bottom) / dy;

                    phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
                }
            }


            //////////////////just replaced with second order central differencing for boundary points

           
            //top boundary points except corners
            j = Ny - 1;

            for (i = 1; i < Nx - 1; i++) {

                diffu = Diff * (   (phi[j][i+1] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   ( - phi[j][i] + phi[j-1][i])/(pow(dy,2))    );

                convec =  ( u[j][i+1] * (phi[j][i] + phi[j][i+1]) / 2 - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (0 - v[j][i] * (phi[j - 1][i] + phi[j][i])/2 ) / dy;

                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }


            //right boundary points except corners
            i = Nx - 1;

            for (j = 1; j < Ny - 1; j++) {

                diffu = Diff * (   (phi[j][i] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - 2*phi[j][i] + phi[j-1][i])/(pow(dy,2))    );

                convec =  ( u[j][i] * (phi[j][i] + phi[j][i]) / 2 - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (v[j + 1][i] * (phi[j][i] + phi[j + 1][i])/2 - v[j][i] * (phi[j - 1][i] + phi[j][i])/2 ) / dy;

                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }


            //left boundary points
            i = 0;

            for (j = 0; j < Ny; j++) {
                phi_new[j][i] = exp( -1 * pow(j * dy - 2.5 * H + 0.5 * dy ,2) );
            }


            //bottom boundary points except corners
            j = 0;

            for (i = 1; i < Nx - 1; i++) {

                diffu = Diff * (   (phi[j][i+1] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - phi[j][i] )/(pow(dy,2))    );

                convec =  ( u[j][i+1] * (phi[j][i] + phi[j][i+1]) / 2 - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (v[j + 1][i] * (phi[j][i] + phi[j + 1][i])/2 - 0 ) / dy;

                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }
   

            //top right corner

            j = Ny - 1;
            i = Nx - 1;

            diffu = Diff * (   (- phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (- phi[j][i] + phi[j-1][i])/(pow(dy,2))    );
            convec =  ( u[j][i] * (phi[j][i] + phi[j][i]) / 2 - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (  0  - v[j][i] * (phi[j - 1][i] + phi[j][i])/2 ) / dy;
            phi_new[j][i] = phi[j][i] + dt * (diffu - convec);


            //bottom right corner

            j = 0;
            i = Nx - 1;

            diffu = Diff * (   ( - phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - phi[j][i] )/(pow(dy,2))    );
            convec =  ( u[j][i] * (phi[j][i] + phi[j][i]) / 2 - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (v[j + 1][i] * (phi[j][i] + phi[j + 1][i])/2 - 0 ) / dy;
            phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
           


            

            //bottom boundary of void
            j = Bj - 1;
            for (i = Bi; i < Bi + void_len; i++) {
                diffu = Diff * (   (phi[j][i+1] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   ( - phi[j][i] + phi[j-1][i])/(pow(dy,2))    );
                convec =  ( u[j][i+1] * (phi[j][i] + phi[j][i+1]) / 2 - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (0 - v[j][i] * (phi[j - 1][i] + phi[j][i])/2 ) / dy;
                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }

            //top boundary of void
            j = Bj + void_len;
            for (i = Bi; i < Bi + void_len; i++) {
                diffu = Diff * (   (phi[j][i+1] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - phi[j][i] )/(pow(dy,2))    );
                convec =  ( u[j][i+1] * (phi[j][i] + phi[j][i+1]) / 2 - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (v[j + 1][i] * (phi[j][i] + phi[j + 1][i])/2 - 0 ) / dy;
                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }
           
            //left boundary of void
            i = Bi - 1;
            for (j = Bj; j < Bj + void_len; j++) {
                diffu = Diff * (   (0 - phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - 2*phi[j][i] + phi[j-1][i])/(pow(dy,2))    );
                convec =  (   0  - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (v[j + 1][i] * (phi[j][i] + phi[j + 1][i])/2 - v[j][i] * (phi[j - 1][i] + phi[j][i])/2 ) / dy;
                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }

            //right boundary of void
            i = Bi + void_len;
            for (j = Bj; j < Bj + void_len; j++) {
                diffu = Diff * (   (phi[j][i+1] - phi[j][i] ) / (pow(dx,2))   +   (phi[j + 1][i] - 2*phi[j][i] + phi[j-1][i])/(pow(dy,2))    );
                convec =  ( u[j][i+1] * (phi[j][i] + phi[j][i+1]) / 2 - 0 ) / dx   +   (v[j + 1][i] * (phi[j][i] + phi[j + 1][i])/2 - v[j][i] * (phi[j - 1][i] + phi[j][i])/2 ) / dy;
                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }






            ///////////// deal with cells that are not on boundary but are next to cells that form boundaries (just use regular second order differencing)
            
            //top 
            j = Ny - 2;
            for (i = 1; i < Nx - 1; i++) {
                diffu = Diff * (   (phi[j][i+1] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - 2*phi[j][i] + phi[j-1][i])/(pow(dy,2))    );
                convec =  ( u[j][i+1] * (phi[j][i] + phi[j][i+1]) / 2 - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (v[j + 1][i] * (phi[j][i] + phi[j + 1][i])/2 - v[j][i] * (phi[j - 1][i] + phi[j][i])/2 ) / dy;
                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }

            //left
            i = 1;
            for (j = 1; j < Ny - 1; j++) {
                diffu = Diff * (   (phi[j][i+1] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - 2*phi[j][i] + phi[j-1][i])/(pow(dy,2))    );
                convec =  ( u[j][i+1] * (phi[j][i] + phi[j][i+1]) / 2 - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (v[j + 1][i] * (phi[j][i] + phi[j + 1][i])/2 - v[j][i] * (phi[j - 1][i] + phi[j][i])/2 ) / dy;
                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }

            //bottom
            j = 1;
            for (i = 1; i < Nx - 1; i++) {
                diffu = Diff * (   (phi[j][i+1] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - 2*phi[j][i] + phi[j-1][i])/(pow(dy,2))    );
                convec =  ( u[j][i+1] * (phi[j][i] + phi[j][i+1]) / 2 - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (v[j + 1][i] * (phi[j][i] + phi[j + 1][i])/2 - v[j][i] * (phi[j - 1][i] + phi[j][i])/2 ) / dy;
                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }

            //right
            i = Nx - 2;
            for (j = 1; j < Ny - 1; j++) {
                diffu = Diff * (   (phi[j][i+1] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - 2*phi[j][i] + phi[j-1][i])/(pow(dy,2))    );
                convec =  ( u[j][i+1] * (phi[j][i] + phi[j][i+1]) / 2 - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (v[j + 1][i] * (phi[j][i] + phi[j + 1][i])/2 - v[j][i] * (phi[j - 1][i] + phi[j][i])/2 ) / dy;
                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }

            //void left
            i = Bi - 2;
            for (j = Bj; j < Bj + void_len; j++) {
                diffu = Diff * (   (phi[j][i+1] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - 2*phi[j][i] + phi[j-1][i])/(pow(dy,2))    );
                convec =  ( u[j][i+1] * (phi[j][i] + phi[j][i+1]) / 2 - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (v[j + 1][i] * (phi[j][i] + phi[j + 1][i])/2 - v[j][i] * (phi[j - 1][i] + phi[j][i])/2 ) / dy;
                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }

            //void right
            i = Bi + void_len + 1;
            for (j = Bj; j < Bj + void_len; j++) {
                diffu = Diff * (   (phi[j][i+1] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - 2*phi[j][i] + phi[j-1][i])/(pow(dy,2))    );
                convec =  ( u[j][i+1] * (phi[j][i] + phi[j][i+1]) / 2 - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (v[j + 1][i] * (phi[j][i] + phi[j + 1][i])/2 - v[j][i] * (phi[j - 1][i] + phi[j][i])/2 ) / dy;
                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }

            //void top
            j = Bj + void_len + 1;
            for (i = Bi; i < Bi + void_len; i++) {
                diffu = Diff * (   (phi[j][i+1] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - 2*phi[j][i] + phi[j-1][i])/(pow(dy,2))    );
                convec =  ( u[j][i+1] * (phi[j][i] + phi[j][i+1]) / 2 - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (v[j + 1][i] * (phi[j][i] + phi[j + 1][i])/2 - v[j][i] * (phi[j - 1][i] + phi[j][i])/2 ) / dy;
                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }

            //void bottom
            j = Bj - 2;
            for (i = Bi; i < Bi + void_len; i++) {
                diffu = Diff * (   (phi[j][i+1] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - 2*phi[j][i] + phi[j-1][i])/(pow(dy,2))    );
                convec =  ( u[j][i+1] * (phi[j][i] + phi[j][i+1]) / 2 - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (v[j + 1][i] * (phi[j][i] + phi[j + 1][i])/2 - v[j][i] * (phi[j - 1][i] + phi[j][i])/2 ) / dy;
                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }


            //inside square, set phi_new = 0

            for (j = Bj; j < Bj + void_len; j++) {                        
                for (i = Bi; i < Bi + void_len; i++) {
                    phi_new[j][i] = 0;
                }
            }


            // update phi everywhere

            for (j = 0; j < Ny; j++) {
                for (i = 0; i < Nx; i++) {
                    phi[j][i] = phi_new[j][i];
                }
            }



        }







        ///////// print periodocially

        if (iter == 1000) {
            print_current_data(iter, u, v, phi, Nx, Ny, convective_method);
            printf("Data was printed for Gauss Seidel method at step %d", iter);
        }

        if (iter == 2000) {
            print_current_data(iter, u, v, phi, Nx, Ny, convective_method);
            printf("Data was printed for Gauss Seidel method at step %d", iter);
        }

        if (iter == 3000) {
            print_current_data(iter, u, v, phi, Nx, Ny, convective_method);
            printf("Data was printed for Gauss Seidel method at step %d", iter);
        }

        if (iter == 5000) {
            print_current_data(iter, u, v, phi, Nx, Ny, convective_method);
            printf("Data was printed for Gauss Seidel method at step %d", iter);
        }

        if (iter == 7000) {
            print_current_data(iter, u, v, phi, Nx, Ny, convective_method);
            printf("Data was printed for Gauss Seidel method at step %d", iter);
        }

        if (iter == 10000) {
            print_current_data(iter, u, v, phi, Nx, Ny, convective_method);
            printf("Data was printed for Gauss Seidel method at step %d", iter);
        }

        if (iter == 15000) {
            print_current_data(iter, u, v, phi, Nx, Ny, convective_method);
            printf("Data was printed for Gauss Seidel method at step %d", iter);
        }

        if (iter == 20000) {
            print_current_data(iter, u, v, phi, Nx, Ny, convective_method);
            printf("Data was printed for Gauss Seidel method at step %d", iter);
        }

        if (iter == 25000) {
            print_current_data(iter, u, v, phi, Nx, Ny, convective_method);
            printf("Data was printed for Gauss Seidel method at step %d", iter);
        }

        if (iter == 30000) {
            print_current_data(iter, u, v, phi, Nx, Ny, convective_method);
            printf("Data was printed for Gauss Seidel method at step %d", iter);
        }

        if (iter == 35000) {
            print_current_data(iter, u, v, phi, Nx, Ny, convective_method);
            printf("Data was printed for Gauss Seidel method at step %d", iter);
        }

        if (iter == 40000) {
            print_current_data(iter, u, v, phi, Nx, Ny, convective_method);
            printf("Data was printed for Gauss Seidel method at step %d", iter);
        }

        if (iter == 45000) {
            print_current_data(iter, u, v, phi, Nx, Ny, convective_method);
            printf("Data was printed for Gauss Seidel method at step %d", iter);
        }

        if (iter == 50000) {
            print_current_data(iter, u, v, phi, Nx, Ny, convective_method);
            printf("Data was printed for Gauss Seidel method at step %d", iter);
        }

        if (iter == 55000) {
            print_current_data(iter, u, v, phi, Nx, Ny, convective_method);
            printf("Data was printed for Gauss Seidel method at step %d", iter);
        }

        


        

        

    }



    print_current_data(iter, u, v, phi, Nx, Ny, convective_method);
    printf("Data was printed for Gauss Seidel method at step %d", iter);



    /* -------------------------------------------------------------------------
    ---------------------------- Freeing Variables -----------------------------
    ------------------------------------------------------------------------- */
   

   
    for (i = 0; i < Ny; i++) {
        free(u[i]);
        free(v[i]);

        free(Hx_old[i]);
        free(Hy_old[i]);

        free(step1_mat_x[i]);
        free(step1_mat_y[i]);

        free(du_s[i]);
        free(dv_s[i]);
        free(du_ss[i]);
        free(dv_ss[i]);

        free(p[i]);
        free(grad_u_star_over_dt[i]);
        free(phi[i]);
        free(phi_new[i]);

    }
   
    free(u);
    free(v);

    free(Hx_old);
    free(Hy_old);




    free(du_s);
    free(dv_s);
    free(du_ss);
    free(dv_ss);

    free(step1_mat_x);
    free(step1_mat_y);

    free(temp_vec_x_small);
    free(temp_vec_x_small2);
    free(temp_vec_x_medium);
    free(temp_vec_x_medium2);
    free(temp_vec_x_long);
    free(temp_vec_x_long2);

    free(temp_vec_y_small);
    free(temp_vec_y_small2);
    free(temp_vec_y_small3);
    free(temp_vec_y_small4);
    free(temp_vec_y_long);
    free(temp_vec_y_long2);

    free(p);
    free(grad_u_star_over_dt);
    free(phi);
    free(phi_new);

    for (i = 0; i < Ny; i++) {
        free(u_star[i]);
    }

    for (i = 0; i < Ny + 1; i++) {
        free(v_star[i]);
    }

    free(u_star);
    free(v_star);


    printf("\n End of script, congrats!");
   
    return 0;



}