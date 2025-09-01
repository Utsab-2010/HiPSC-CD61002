#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_ITER 1000
#define TOL 1e-6

void mat_vec_mul(int n, double** A, double* x, double* y) {
    for (int i = 0; i < n; i++) {
        y[i] = 0.0;
        for (int j = 0; j < n; j++) {
            y[i] += A[i][j] * x[j];
        }
    }
}

double norm_diff(int n, double* x_new, double* x_old) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        double diff = x_new[i] - x_old[i];
        sum += diff * diff;
    }
    return sqrt(sum);
}

void jacobi_solver(int n, double** A, double* B, double* x) {
    double *x_new = (double*)malloc(n * sizeof(double));

    for (int iter = 0; iter < MAX_ITER; iter++) {
        for (int i = 0; i < n; i++) {
            double s = 0.0;
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    s += A[i][j] * x[j];
                }
            }
            x_new[i] = (B[i] - s) / A[i][i];
        }

        if (norm_diff(n, x_new, x) < TOL) {
            printf("Jacobi method converged in %d iterations.\n", iter);
            break;
        }

        // Copy x_new to x
        for (int i = 0; i < n; i++) {
            x[i] = x_new[i];
        }
    }

    free(x_new);
}

int main() {
    // Example size
    int n = 3;

    // Allocate memory for matrix A
    double** A = (double**)malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++)
        A[i] = (double*)malloc(n * sizeof(double));

    // Example matrix A
    A[0][0] = 4; A[0][1] = 1; A[0][2] = 2;
    A[1][0] = 3; A[1][1] = 5; A[1][2] = 1;
    A[2][0] = 1; A[2][1] = 1; A[2][2] = 3;

    // Example vector B
    double B[3] = {4,7,3};
    double x[3] = {0,0,0};  // Initial guess

    jacobi_solver(n, A, B, x);

    printf("Solution:\n");
    for (int i = 0; i < n; i++) {
        printf("%f\n", x[i]);
    }

    // Free memory
    for (int i = 0; i < n; i++)
        free(A[i]);
    free(A);

    return 0;
}
