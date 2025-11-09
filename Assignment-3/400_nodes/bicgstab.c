#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define MAX_ITER 10000
#define TOLERANCE 1e-6

// Function to read matrix from Kmat.txt
int read_matrix(const char* filename, double** K, int* n) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        printf("Error: Cannot open %s\n", filename);
        return -1;
    }
    
    // First, count the number of elements to determine matrix size
    int count = 0;
    double temp;
    while (fscanf(file, "%lf", &temp) == 1) {
        count++;
    }
    
    *n = (int)sqrt(count);
    printf("Matrix size detected: %d x %d\n", *n, *n);
    
    // Allocate memory for matrix
    *K = (double*)malloc((*n) * (*n) * sizeof(double));
    if (!*K) {
        printf("Error: Memory allocation failed\n");
        fclose(file);
        return -1;
    }
    
    // Rewind and read the matrix
    rewind(file);
    int i;
    for (i = 0; i < (*n) * (*n); i++) {
        if (fscanf(file, "%lf", &(*K)[i]) != 1) {
            printf("Error reading matrix element %d\n", i);
            free(*K);
            fclose(file);
            return -1;
        }
    }
    
    fclose(file);
    return 0;
}

// Function to read vector from Fvec.txt
int read_vector(const char* filename, double** F, int* n) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        printf("Error: Cannot open %s\n", filename);
        return -1;
    }
    
    // Count elements
    int count = 0;
    double temp;
    while (fscanf(file, "%lf", &temp) == 1) {
        count++;
    }
    
    *n = count;
    printf("Vector size detected: %d\n", *n);
    
    // Allocate memory
    *F = (double*)malloc(*n * sizeof(double));
    if (!*F) {
        printf("Error: Memory allocation failed\n");
        fclose(file);
        return -1;
    }
    
    // Rewind and read vector
    rewind(file);
    int i;
    for (i = 0; i < *n; i++) {
        if (fscanf(file, "%lf", &(*F)[i]) != 1) {
            printf("Error reading vector element %d\n", i);
            free(*F);
            fclose(file);
            return -1;
        }
    }
    
    fclose(file);
    return 0;
}

// Matrix-vector multiplication: y = A*x
void matvec(double* A, double* x, double* y, int n) {
    int i, j;
    for (i = 0; i < n; i++) {
        y[i] = 0.0;
        for (j = 0; j < n; j++) {
            y[i] += A[i*n + j] * x[j];
        }
    }
}

// Vector operations
double dot_product(double* a, double* b, int n) {
    double result = 0.0;
    int i;
    for (i = 0; i < n; i++) {
        result += a[i] * b[i];
    }
    return result;
}

void vector_copy(double* src, double* dest, int n) {
    int i;
    for (i = 0; i < n; i++) {
        dest[i] = src[i];
    }
}

void vector_axpy(double* y, double a, double* x, int n) {
    // y = y + a*x
    int i;
    for (i = 0; i < n; i++) {
        y[i] += a * x[i];
    }
}

void vector_scale(double* x, double a, int n) {
    int i;
    for (i = 0; i < n; i++) {
        x[i] *= a;
    }
}

void vector_add(double* a, double* b, double* c, double s, int n){
    int i;
    for(i = 0; i < n; i++){
        c[i] = a[i] + s * b[i];
    }
}
double vector_norm(double* x, int n) {
    return sqrt(dot_product(x, x, n));
}

// Simple Jacobi preconditioner (diagonal preconditioning)
void jacobi_preconditioner(double* A, double* r, double* z, int n) {
    int i;
    for (i = 0; i < n; i++) {
        if (fabs(A[i*n + i]) > 1e-12) {
            z[i] = r[i] / A[i*n + i];
        } else {
            z[i] = r[i]; // If diagonal is zero, no preconditioning
        }
    }
}

// Conjugate Gradient method with preconditioning (BigSTABC-like algorithm)
int bigstabc_solve(double* A, double* b, double* x, int n) {
    printf("Starting BiCGSTAB solver...\n");
    
    // Allocate working vectors
    double* r = (double*)malloc(n * sizeof(double));  // residual
    double* r_old = (double*)malloc(n*sizeof(double)); // old residual
    double* r0 = (double*)malloc(n*sizeof(double)); // initial residual   
    double* s = (double*)malloc(n * sizeof(double));  // preconditioned residual
    double* p = (double*)malloc(n * sizeof(double));  // search direction
    double* Ap = (double*)malloc(n * sizeof(double)); // A*p
    double* As = (double*)malloc(n * sizeof(double)); // A*s
    double* temp = (double*)malloc(n * sizeof(double)); // temporary vector for safe operations
    
    if (!r || !r_old || !r0 || !s || !p || !Ap || !As || !temp) {
        printf("Error: Memory allocation failed in solver\n");
        return -1;
    }
    
    // Initialize: r = b - A*x (x starts as zero vector)
    matvec(A, x, temp, n);
    vector_add(b,temp,r,-1.0,n);

    for(int i = 0; i < n; i++){
        r0[i] = (rand() % 100)/100.0;  // Random values between 0 and 99
    }
    // int i;
    // for (i = 0; i < n; i++) {
    //     r[i] = b[i] - r[i];
    // }
    // vector_copy(r0, r, n);
    vector_copy(r, p, n);  // Initialize p = r0
    
    double initial_residual = vector_norm(r, n);
    printf("Initial residual norm: %e\n", initial_residual);
    printf("Initial r0 norm: %e\n", vector_norm(r0, n));
    printf("Initial p norm: %e\n", vector_norm(p, n));
    
    if (initial_residual < TOLERANCE) {
        printf("Initial guess is already converged!\n");
        free(r); free(r_old); free(r0); free(s); free(p); free(Ap); free(As); free(temp);
        return 0;
    }
    
    // Main BiCGSTAB iteration loop
    double time_start  = clock();

    int iter;
    for (iter = 0; iter < MAX_ITER; iter++) {
        // Ap = A * p
        matvec(A, p, Ap, n);
        
        // Debug information
        double p_norm = vector_norm(p, n);
        double ap_norm = vector_norm(Ap, n);
        double r0_norm = vector_norm(r0, n);
        
        if (iter % 100 == 0) {
            printf("Iter %d: ||p|| = %e, ||Ap|| = %e, ||r0|| = %e\n", 
                   iter, p_norm, ap_norm, r0_norm);
        }

        double  ap_dot_r0= dot_product(Ap, r0, n);
        
        if (fabs(ap_dot_r0) < 1e-15) {
            printf("BiCGSTAB breakdown: Ap dot r0 = %e at iteration %d\n", ap_dot_r0, iter);
            printf("This indicates the algorithm has broken down.\n");
            printf("Trying to restart with current residual as new r0...\n");
            
            printf("DEBUG: Before restart - r norm: %e\n", vector_norm(r, n));
            printf("DEBUG: Before restart - r0 norm: %e\n", vector_norm(r0, n));
            printf("DEBUG: Before restart - p norm: %e\n", vector_norm(p, n));
            printf("DEBUG: Before restart - Ap norm: %e\n", vector_norm(Ap, n));
            
            // Restart: use current residual as new r0
            vector_copy(r, r0, n);
            vector_copy(r, p, n);
            
            // CRUCIAL FIX: We must recalculate Ap since p has changed!
            matvec(A, p, Ap, n);
            ap_dot_r0 = dot_product(Ap, r0, n);
            
            printf("DEBUG: After restart - new r0 norm: %e\n", vector_norm(r0, n));
            printf("DEBUG: After restart - new p norm: %e\n", vector_norm(p, n));
            printf("DEBUG: After restart - new Ap norm: %e\n", vector_norm(Ap, n));
            printf("DEBUG: After restart - new ap_dot_r0: %e\n", ap_dot_r0);
            
            if (fabs(ap_dot_r0) < 1e-15) {
                printf("Restart failed. Terminating.\n");
                break;
            }
            printf("Restart successful. Continuing...\n");
        }
        
        double r_dot_r0 = dot_product(r, r0, n);
        double alpha = r_dot_r0 / ap_dot_r0;

        vector_add(r, Ap, s, -alpha, n);  // s = r - alpha * Ap

        matvec(A, s, As, n);
        double as_dot_as = dot_product(As, As, n);
        if (fabs(as_dot_as) < 1e-15) {
            printf("Warning: As dot As is nearly zero at iteration %d\n", iter);
            break;
        }
        
        double omega = dot_product(As, s, n) / as_dot_as;

        // Update x safely using temporary vector
        vector_add(x, p, temp, alpha, n);   // temp = x + alpha * p
        vector_add(temp, s, x, omega, n);   // x = temp + omega * s

        vector_copy(r, r_old, n);
        vector_add(s, As, r, -omega, n);    // r = s - omega * As

        double r_new_dot_r0 = dot_product(r, r0, n);
        double r_old_dot_r0 = dot_product(r_old, r0, n);
        
        if (fabs(r_old_dot_r0) < 1e-15 || fabs(omega) < 1e-15) {
            printf("Warning: Division by nearly zero at iteration %d\n", iter);
            printf("r_old_dot_r0 = %e, omega = %e\n", r_old_dot_r0, omega);
            break;
        }
        
        double beta = (r_new_dot_r0 / r_old_dot_r0) * (alpha / omega);
        
        double residual_norm = vector_norm(r, n);
        double relative_residual = residual_norm / initial_residual;
        
        // Check for NaN
        if (isnan(residual_norm) || isnan(relative_residual)) {
            printf("Error: NaN detected at iteration %d\n", iter);
            printf("alpha = %e, omega = %e, beta = %e\n", alpha, omega, beta);
            printf("r_dot_r0 = %e, ap_dot_r0 = %e\n", r_dot_r0, ap_dot_r0);
            break;
        }
        
        // if (iter % 100 == 0 || relative_residual < TOLERANCE) {
        //     printf("Iteration %d: residual = %e, relative = %e\n", 
        //            iter, residual_norm, relative_residual);
        // }
        
        if (relative_residual < TOLERANCE) {
            double total_time = (clock() - time_start) / CLOCKS_PER_SEC;
            printf("Total time taken: %f seconds\n", total_time);
            printf("Converged in %d iterations\n", iter + 1);
            free(r); free(r_old); free(r0); free(s); free(p); free(Ap); free(As); free(temp);
            return iter + 1;
        }
        
        if (fabs(residual_norm) < 1e-16) {
            printf("Warning: residual is nearly zero at iteration %d\n", iter);
            break;
        }
        
        // Update p safely: p = r + beta * (p - omega * Ap)
        vector_add(p, Ap, temp, -omega, n);  // temp = p - omega * Ap
        vector_add(r, temp, p, beta, n);     // p = r + beta * temp
    }
    double total_time = (clock() - time_start) / CLOCKS_PER_SEC;
    printf("Total time taken: %f seconds\n", total_time);
    printf("Warning: Maximum iterations reached without convergence\n");
    free(r); free(r_old); free(r0); free(s); free(p); free(Ap); free(As); free(temp);
    return -1;
}

// Write solution to file
void write_solution(double* x, int n, const char* filename) {
    FILE* file = fopen(filename, "w");
    if (!file) {
        printf("Error: Cannot create output file %s\n", filename);
        return;
    }
    
    fprintf(file, "# Solution vector (node displacements)\n");
    int i;
    for (i = 0; i < n; i++) {
        fprintf(file, "Node %d: %15.8e\n", i+1, x[i]);
    }
    
    fclose(file);
    printf("Solution written to %s\n", filename);
}

// Verify solution by computing residual
void verify_solution(double* A, double* b, double* x, int n) {
    double* Ax = (double*)malloc(n * sizeof(double));
    if (!Ax) return;
    
    matvec(A, x, Ax, n);
    
    double residual_norm = 0.0;
    double b_norm = 0.0;
    int i;
    
    for (i = 0; i < n; i++) {
        double residual = b[i] - Ax[i];
        residual_norm += residual * residual;
        b_norm += b[i] * b[i];
    }
    
    residual_norm = sqrt(residual_norm);
    b_norm = sqrt(b_norm);
    
    printf("\nSolution verification:\n");
    printf("||b - A*x|| = %e\n", residual_norm);
    printf("||b||       = %e\n", b_norm);
    printf("Relative error = %e\n", residual_norm / b_norm);
    
    free(Ax);
}

int main() {
    printf("BigSTABC Linear System Solver\n");
    printf("=============================\n\n");
    
    double* K = NULL;  // Stiffness matrix
    double* F = NULL;  // Force vector
    double* x = NULL;  // Solution vector
    int n_matrix, n_vector;
    
    // Read matrix and vector files
    if (read_matrix("Kmat.txt", &K, &n_matrix) != 0) {
        printf("Failed to read matrix file\n");
        return 1;
    }
    
    if (read_vector("Fvec.txt", &F, &n_vector) != 0) {
        printf("Failed to read vector file\n");
        free(K);
        return 1;
    }
    
    // Check dimensions match
    if (n_matrix != n_vector) {
        printf("Error: Matrix size (%d) doesn't match vector size (%d)\n", 
               n_matrix, n_vector);
        free(K);
        free(F);
        return 1;
    }
    
    int n = n_matrix;
    printf("System size: %d x %d\n\n", n, n);
    
    // Allocate solution vector (initialize to zero)
    x = (double*)calloc(n, sizeof(double));
    if (!x) {
        printf("Error: Failed to allocate solution vector\n");
        free(K);
        free(F);
        return 1;
    }
    
    // Solve the system K*x = F
    int iterations = bigstabc_solve(K, F, x, n);
    
    if (iterations > 0) {
        printf("\nSolution completed successfully!\n");
        
        // Verify the solution
        verify_solution(K, F, x, n);
        
        // Write solution to file
        write_solution(x, n, "solution.txt");
        
        // Print some solution statistics
        double max_displacement = 0.0;
        int max_node = 0;
        int i;
        for (i = 0; i < n; i++) {
            if (fabs(x[i]) > max_displacement) {
                max_displacement = fabs(x[i]);
                max_node = i + 1;
            }
        }
        
        printf("\nSolution statistics:\n");
        printf("Maximum displacement: %e at node %d\n", max_displacement, max_node);
        printf("Solution vector norm: %e\n", vector_norm(x, n));
        
    } else {
        printf("Solution failed to converge\n");
    }
    
    // Cleanup
    free(K);
    free(F);
    free(x);
    
    return (iterations > 0) ? 0 : 1;
}