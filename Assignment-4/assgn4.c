/* mpi_jacobi_sor.c
   Parallel weighted-Jacobi (SOR-like) solver for 2D Laplace using row-wise domain decomposition.
   Compile: mpicc -O3 -std=c99 -lm mpi_jacobi_sor.c -o mpi_jacobi_sor
   Run: mpirun -np <P> ./mpi_jacobi_sor <n> <omega> <tol> <max_iter>
     where n = number of interior points per row (so total unknowns = n*n)
     omega = relaxation factor (0 < omega <= 2), e.g. 1.8818 from your single-process run
     tol = stopping tolerance on max update (e.g. 1e-6)
     max_iter = max iterations (e.g. 10000)
*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int min(int a,int b){ return a<b? a:b; }

int main(int argc, char *argv[]){
    MPI_Init(&argc,&argv);
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

    if(argc < 5){
        if(rank==0) printf("Usage: %s <n> <omega> <tol> <max_iter>\n  n = interior points per row (total unknowns = n*n)\n", argv[0]);
        MPI_Finalize();
        return 0;
    }

    int n = atoi(argv[1]);             // interior points per row
    double omega = atof(argv[2]);      // relaxation factor for weighted Jacobi
    double tol = atof(argv[3]);
    int max_iter = atoi(argv[4]);

    int N = n;                         // keep name consistent
    // Partition rows among processes (each process gets local_n interior rows)
    int base = N / nprocs;
    int rem  = N % nprocs;

    int local_n = base + (rank < rem ? 1 : 0); // number of *interior* rows assigned to this rank
    // global start row index (0..N-1) for this rank
    int start_row = rank * base + (rank < rem ? rank : rem);
    int end_row = start_row + local_n - 1;

    // We'll store grid with ghost rows top and bottom: (local_n + 2) x N
    int rows_with_ghosts = local_n + 2;
    // allocate 2 arrays: old and new
    double *u_old = (double*)calloc(rows_with_ghosts * N, sizeof(double));
    double *u_new = (double*)calloc(rows_with_ghosts * N, sizeof(double));
    if(!u_old || !u_new){ printf("Rank %d: allocation failed\n",rank); MPI_Abort(MPI_COMM_WORLD,1); }

    // RHS b is zero except top BC contributes to some rows. For simplicity, incorporate BC into "b_term"
    // Top boundary (global row N corresponds to top Dirichlet = 1.0). Note: interior rows are 0..N-1 with top at row N-1.
    // Here we assume top boundary T=1 (y=1); other boundaries = 0.

    // initialize u_old to 0 already by calloc
    // Prepare neighbor ranks for communication
    int rank_up = (rank == nprocs-1) ? MPI_PROC_NULL : rank + 1;   // rank holding rows with larger global row index
    int rank_down = (rank == 0) ? MPI_PROC_NULL : rank - 1;       // rank holding rows with smaller global row index

    MPI_Barrier(MPI_COMM_WORLD);
    double t_start = MPI_Wtime();

    int iter;
    double global_max_diff = 0.0;
    for(iter=0; iter < max_iter; iter++){
        // Exchange ghost rows: send my first interior row upward? careful with indexing
        // Our local array layout: row 0 = ghost-top, rows 1..local_n = real, row local_n+1 = ghost-bottom
        // Global row indexing: 0 (bottom) ... N-1 (top). We used start_row..end_row accordingly.
        // Send my top real row (row local_n) to rank_up as their ghost-bottom, receive their top real row into my ghost-top? Let's use Sendrecv for both directions.

        // Send my top real row to rank_up -> it will become their ghost-bottom.
        MPI_Status st;
        // send top real row (index local_n) to rank_up, receive neighbor top ghost (my ghost-top at index 0) from rank_up
        MPI_Sendrecv(
            &u_old[ (local_n)*N ],  // send my top real row
            N, MPI_DOUBLE, rank_up, 0,
            &u_old[ 0 ], N, MPI_DOUBLE, rank_up, 1,
            MPI_COMM_WORLD, &st
        );
        // send bottom real row (index 1) to rank_down, receive neighbor bottom ghost into my last index (local_n+1)
        MPI_Sendrecv(
            &u_old[ 1*N ], N, MPI_DOUBLE, rank_down, 1,
            &u_old[ (local_n+1)*N ], N, MPI_DOUBLE, rank_down, 0,
            MPI_COMM_WORLD, &st
        );

        // compute local updates using weighted Jacobi stencil
        double local_max_diff = 0.0;
        for(int i_local = 1; i_local <= local_n; ++i_local){
            int global_row = start_row + (i_local - 1); // which global interior row this is (0 .. N-1)
            for(int j=0; j < N; ++j){
                // neighbors:
                double left  = (j==0)   ? 0.0 : u_old[ i_local*N + (j-1) ];
                double right = (j==N-1) ? 0.0 : u_old[ i_local*N + (j+1) ];
                double down  = u_old[ (i_local-1)*N + j ];  // bottom neighbor (could be ghost row or interior)
                double up    = u_old[ (i_local+1)*N + j ];  // top neighbor

                // contribution from boundary top = 1.0 -> only matters for points adjacent to top boundary (global_row == N-1)
                double b_term = 0.0;
                if(global_row == N-1){
                    // top boundary is at global row N (outside interior), with T=1
                    // discretization for Laplace: average of neighbors => RHS term for top neighbor becomes +1 (we treat boundaries as known)
                    b_term = 1.0;
                }

                double jacobi = 0.25 * (left + right + up + down + b_term);
                // Weighted Jacobi:
                double newval = (1.0 - omega) * u_old[ i_local*N + j ] + omega * jacobi;

                u_new[ i_local*N + j ] = newval;
                double diff = fabs(newval - u_old[ i_local*N + j ]);
                if(diff > local_max_diff) local_max_diff = diff;
            }
        }

        // compute global max
        MPI_Allreduce(&local_max_diff, &global_max_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        // swap pointers
        double *tmp = u_old; u_old = u_new; u_new = tmp;

        if(global_max_diff <= tol) break;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double t_end = MPI_Wtime();
    double elapsed = t_end - t_start;

    // gather timings to rank 0: we will collect the max elapsed (worst-case)
    double max_elapsed;
    MPI_Reduce(&elapsed, &max_elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if(rank == 0){
        printf("n=%d (N*N=%d unknowns), procs=%d, omega=%.6f, tol=%.1e, iter=%d, time=%.6f s\n",
               N, N*N, nprocs, omega, tol, iter+1, max_elapsed);
    }

    // optionally write the solution slice from rank 0 to file (omitted here for brevity)
    free(u_old); free(u_new);
    MPI_Finalize();
    return 0;
}