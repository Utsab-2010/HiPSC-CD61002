#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

int p(int j, int x, int number_nodes_x) {
    return j * number_nodes_x + x;
}

int main(int argc, char* argv[]) {

    // MPI init
    MPI_Init(&argc, &argv);
    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    const double domain_len_x = 1.0;
    const double domain_len_y = 1.0;
    int number_steps_x;
    int number_steps_y;
    double omega = 1.25; // SOR relaxation factor
    if (myrank == 0) {
        printf("Enter number of steps in x and y directions and omega:\n");
        if (scanf("%d %d %lf", &number_steps_x, &number_steps_y, &omega) != 3) {
            fprintf(stderr, "Rank 0: Invalid input\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    // Broadcast number of steps to all processes
    MPI_Bcast(&number_steps_x, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&number_steps_y, 1, MPI_INT, 0, MPI_COMM_WORLD);

    double step_size_x = domain_len_x / (double)number_steps_x;
    double step_size_y = domain_len_y / (double)number_steps_y;

    double coeff_x = 1.0 / (step_size_x * step_size_x);
    double coeff_y = 1.0 / (step_size_y * step_size_y);

    // int Nx_global = number_steps_x;
    // int Ny_global = number_steps_y;

    // double omega = 2.0 / (1.0 + sin(M_PI / (double)( (Nx_global > Ny_global) ?
    //                                                 Nx_global : Ny_global )));
    // if (myrank == 0) {
    //     printf("Using optimal omega = %f\n", omega);
    // }

    // 1D domain decomposition along x (in terms of steps)
    int number_steps_x_subdomain = number_steps_x / nproc;
    if (myrank < number_steps_x % nproc) {
        number_steps_x_subdomain += 1;
    }

    // Local number of nodes in x (steps + 1), global in y
    int number_nodes_x = number_steps_x_subdomain + 1;
    int number_nodes_y = number_steps_y + 1;

    // (node_x_start / node_x_end are computed but not used further)
    int node_x_start, node_x_end;
    if (myrank < number_steps_x % nproc) {
        node_x_start = myrank * number_nodes_x;
        node_x_end   = (myrank + 1) * number_nodes_x;
    } else {
        node_x_start = myrank * number_nodes_x + number_steps_x % nproc;
        node_x_end   = (myrank + 1) * number_nodes_x + number_steps_x % nproc;
    }

    // Allocate and initialize matrices and vectors
    int i, j, k;
    int npoints = number_nodes_x * number_nodes_y;

    double **A    = (double**) malloc(npoints * sizeof(double*));
    double *b     = (double*)  malloc(npoints * sizeof(double));
    double *x     = (double*)  malloc(npoints * sizeof(double));
    double *x_old = (double*)  malloc(npoints * sizeof(double));

    if (A == NULL || b == NULL || x == NULL || x_old == NULL) {
        printf("Rank %d: Memory allocation failed (A/b/x/x_old)\n", myrank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    for (j = 0; j < npoints; j++) {
        A[j] = (double*) malloc(npoints * sizeof(double));
        if (A[j] == NULL) {
            printf("Rank %d: Memory allocation failed for A[%d]\n", myrank, j);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        b[j]     = 0.0;
        x[j]     = 0.0;
        x_old[j] = 0.0;
        for (k = 0; k < npoints; k++) {
            A[j][k] = 0.0;
        }
    }

    // ----------------------------------------------------------------------
    // Build the matrix A and RHS b
    // Top boundary = 1, other boundaries = 0
    // ----------------------------------------------------------------------

    // Left boundary (global x = 0) only on rank 0
    if (myrank == 0) {
        for (j = 0; j < number_nodes_y; j++) { // leftmost local node is global left boundary
            int idx = p(j, 0, number_nodes_x);
            A[idx][idx] = 1.0;
            b[idx]      = 0.0;
        }
    }

    // Right boundary (global x = domain_len_x) only on last rank
    if (myrank == nproc - 1) {
        for (j = 0; j < number_nodes_y; j++) { // rightmost local node is global right boundary
            int idx = p(j, number_nodes_x - 1, number_nodes_x);
            A[idx][idx] = 1.0;
            b[idx]      = 0.0;
        }
    }

    // Bottom boundary (y = 0) – local row j = 0 on every rank
    for (i = 0; i < number_nodes_x; i++) {
        int idx = p(0, i, number_nodes_x);
        A[idx][idx] = 1.0;
        b[idx]      = 0.0;
    }

    // Top boundary (y = domain_len_y) – local row j = number_nodes_y - 1 on every rank
    for (i = 0; i < number_nodes_x; i++) {
        int idx = p(number_nodes_y - 1, i, number_nodes_x);
        A[idx][idx] = 1.0;
        b[idx]      = 1.0;  // top boundary = 1
    }

    // Interior points
    for (j = 1; j < number_nodes_y - 1; j++) {
        for (i = 1; i < number_nodes_x - 1; i++) {
            int row = p(j, i, number_nodes_x);
            A[row][row]                         = -2.0 * (coeff_x + coeff_y);
            A[row][p(j, i - 1, number_nodes_x)] =  coeff_x;
            A[row][p(j, i + 1, number_nodes_x)] =  coeff_x;
            A[row][p(j - 1, i, number_nodes_x)] =  coeff_y;
            A[row][p(j + 1, i, number_nodes_x)] =  coeff_y;
        }
    }

    // ----------------------------------------------------------------------
    // Buffers for halo exchange (vertical columns)
    // ----------------------------------------------------------------------
    double *Send_to_R   = (double*) malloc(number_nodes_y * sizeof(double));
    double *Send_to_L   = (double*) malloc(number_nodes_y * sizeof(double));
    double *Recv_from_R = (double*) malloc(number_nodes_y * sizeof(double));
    double *Recv_from_L = (double*) malloc(number_nodes_y * sizeof(double));

    if (!Send_to_R || !Send_to_L || !Recv_from_R || !Recv_from_L) {
        printf("Rank %d: Memory allocation failed for halo buffers\n", myrank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int max_iterations      = 10000;
    double point_error;
    double global_error     = 0.0;
    double global_error_max = 1e-5;
    double local_error;
    double sum;

    double mytime;
    // double omega = 1.92; // SOR relaxation factor

    MPI_Barrier(MPI_COMM_WORLD);
    mytime = MPI_Wtime();

    // Neighbor ranks for halo exchange
    int left  = (myrank == 0)        ? MPI_PROC_NULL : myrank - 1;
    int right = (myrank == nproc-1)  ? MPI_PROC_NULL : myrank + 1;

    for (k = 1; k <= max_iterations; k++) {

        // Prepare halo data from current iterate x_old
        if (right != MPI_PROC_NULL) {
            for (j = 0; j < number_nodes_y; j++) {
                // send near-right interior column
                Send_to_R[j] = x_old[p(j, number_nodes_x - 3, number_nodes_x)];
            }
        }

        if (left != MPI_PROC_NULL) {
            for (j = 0; j < number_nodes_y; j++) {
                // send near-left interior column
                Send_to_L[j] = x_old[p(j, 2, number_nodes_x)];
            }
        }

        // ------------------------------------------------------------------
        // Halo exchange using MPI_Sendrecv (safe, no deadlock, no bad ranks)
        // ------------------------------------------------------------------

        // Exchange data with left neighbor (Recv_from_L gets neighbor's right-side data)
        MPI_Sendrecv(
            Send_to_R,             // sendbuf
            number_nodes_y,
            MPI_DOUBLE,
            right,                 // dest
            0,                     // tag
            Recv_from_L,           // recvbuf
            number_nodes_y,
            MPI_DOUBLE,
            left,                  // source
            0,                     // tag
            MPI_COMM_WORLD,
            MPI_STATUS_IGNORE
        );

        // Exchange data with right neighbor (Recv_from_R gets neighbor's left-side data)
        MPI_Sendrecv(
            Send_to_L,
            number_nodes_y,
            MPI_DOUBLE,
            left,
            1,
            Recv_from_R,
            number_nodes_y,
            MPI_DOUBLE,
            right,
            1,
            MPI_COMM_WORLD,
            MPI_STATUS_IGNORE
        );

        // ------------------------------------------------------------------
        // Use received halo values to update boundary-like RHS terms (if desired)
        // ------------------------------------------------------------------
        if (myrank != 0 && myrank != nproc - 1) {
            for (j = 1; j < number_nodes_y - 1; j++) {
                b[p(j, 0, number_nodes_x)]                  = Recv_from_L[j];
                b[p(j, number_nodes_x - 1, number_nodes_x)] = Recv_from_R[j];
            }
        } else if (myrank == 0 && right != MPI_PROC_NULL) {
            for (j = 1; j < number_nodes_y - 1; j++) {
                b[p(j, number_nodes_x - 1, number_nodes_x)] = Recv_from_R[j];
            }
        } else if (myrank == nproc - 1 && left != MPI_PROC_NULL) {
            for (j = 1; j < number_nodes_y - 1; j++) {
                b[p(j, 0, number_nodes_x)]                  = Recv_from_L[j];
            }
        }

        // ------------------------------------------------------------------
        // SOR iteration on local system
        // ------------------------------------------------------------------
        local_error = 0.0;

        for (j = 0; j < npoints; j++) {
            if (fabs(A[j][j]) < 1e-14) {
                // skip if diagonal is zero (should not happen for proper interior/boundary)
                continue;
            }

            sum = b[j];

            for (i = 0; i < npoints; i++) {
                if (A[j][i] == 0.0 || i == j) continue;

                if (i < j) {
                    // use the new value in this iteration (Gauss-Seidel part)
                    sum -= A[j][i] * x[i];
                } else {
                    // use the old value
                    sum -= A[j][i] * x_old[i];
                }
            }

            // Gauss–Seidel value
            sum = sum / A[j][j];

            // SOR relaxation
            sum = (1.0 - omega) * x_old[j] + omega * sum;

            point_error = fabs(sum - x_old[j]);
            if (point_error > local_error) {
                local_error = point_error;
            }

            x[j] = sum;
        }

        // Global max error across all processes
        MPI_Allreduce(&local_error, &global_error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        if (myrank == 0) {
            printf("Iteration %d: Global Error = %e\n", k, global_error);
        }

        if (global_error < global_error_max) {
            if (myrank == 0) {
                printf("Converged in %d iterations with error %e\n", k, global_error);
            }
            break;
        }

        // Update old solution
        for (j = 0; j < npoints; j++) {
            x_old[j] = x[j];
        }
    }

    mytime = MPI_Wtime() - mytime;
    if (myrank == 0) {
        printf("Total time taken: %f seconds\n", mytime);
    }

    // Free memory
    for (j = 0; j < npoints; j++) {
        free(A[j]);
    }
    free(A);
    free(b);
    free(x);
    free(x_old);
    free(Send_to_R);
    free(Send_to_L);
    free(Recv_from_R);
    free(Recv_from_L);

    MPI_Finalize();
    return 0;
}
