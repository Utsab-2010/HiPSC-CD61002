#include <math.h>
#include <stdlib.h>
#include <stdio.h>

int p(int j, int x, int number_nodes_x){
    return j * number_nodes_x + x;
}

int main(int argc, char* argv[]){
    const double domain_len_x = 1;
    const double domain_len_y = 1;
    int number_steps_x ;
    int number_steps_y ;
    // printf("Enter number of steps in x and y directions:\n");
    scanf("%d %d",&number_steps_x,&number_steps_y);
    double step_size_x = domain_len_x/number_steps_x;
    double step_size_y = domain_len_y/number_steps_y;

    double coeff_x = 1/(step_size_x*step_size_x);
    double coeff_y = 1/(step_size_y*step_size_y);

    int number_nodes_x = number_steps_x + 1;
    int number_nodes_y = number_steps_y + 1;
    

    //allocate and initialize matrices and vectors
    int i,j,k;
    int npoints = number_nodes_x * number_nodes_y;
    double** A = (double**) malloc(npoints * sizeof(double*));
    double* b = (double*) malloc(npoints * sizeof(double));
    double* x = (double*) malloc(npoints * sizeof(double));
    double* x_old = (double*) malloc(npoints * sizeof(double));

    for(int j=0;j<npoints;j++){
        A[j] = (double*) malloc(npoints * sizeof(double));
        b[j] = 0.0;
        x[j] = 0.0;
        x_old[j] = 0.0;
        for(int k=0;k<npoints;k++){
            A[j][k] = 0.0;
        }
    }

    //build the matrices top boundary = 1, other boundaries = 0

    for(j=0;j<number_nodes_y;j++){ //leftmost boundary
        A[p(j,0,number_nodes_x)][p(j,0,number_nodes_x)] = 1.0;
        b[p(j,0,number_nodes_x)] = 0.0;
    }
    for(j=0;j<number_nodes_y;j++){ //rightmost boundary
        A[p(j,number_nodes_x-1,number_nodes_x)][p(j,number_nodes_x-1,number_nodes_x)] = 1.0;
        b[p(j,number_nodes_x-1,number_nodes_x)] = 0.0;
    }
    for(i=1;i<number_nodes_x-1;i++){ //bottom boundary
        A[p(0,i,number_nodes_x)][p(0,i,number_nodes_x)] = 1.0;
        b[p(0,i,number_nodes_x)] = 0.0;
    }
    for(i=1;i<number_nodes_x-1;i++){ //top boundary
        A[p(number_nodes_y-1,i,number_nodes_x)][p(number_nodes_y-1,i,number_nodes_x)] = 1.0;
        b[p(number_nodes_y-1,i,number_nodes_x)] = 1.0;
    }

    for(j=1;j<number_nodes_y-1;j++){
        for(i=1;i<number_nodes_x-1;i++){
            A[p(j,i,number_nodes_x)][p(j,i,number_nodes_x)] = -2.0*(coeff_x + coeff_y);
            A[p(j,i,number_nodes_x)][p(j,i-1,number_nodes_x)] = coeff_x;
            A[p(j,i,number_nodes_x)][p(j,i+1,number_nodes_x)] = coeff_x;
            A[p(j,i,number_nodes_x)][p(j-1,i,number_nodes_x)] = coeff_y;
            A[p(j,i,number_nodes_x)][p(j+1,i,number_nodes_x)] = coeff_y;
        }
    }

    int max_iterations = 10000;
    double row_error,point_error,global_error,global_error_ma=1e-8,sum ;
    double omega = 1.25; //relaxation factor
    for(k=1;k<=max_iterations;k++){
        global_error = 0.0;
        printf("Iteration %d\n",k);
        for(j=0;j<npoints;j++){
            sum = 0.0;

            for(i=0;i<npoints;i++){
                
                if (A[j][i]== 0.0) continue;
                
                if (j<i){
                    sum -= A[j][i] * x[i];
                }
                else if (j>i){
                    sum -= A[j][i] * x_old[i];
                }
                sum+=b[j];
                sum *= omega/A[j][j];
                sum += (1.0 - omega) * x_old[j];
            }
            point_error = fabs(sum - x_old[j]);
            if (point_error > global_error){
                global_error = point_error;    
            }   
            x[j] = sum;
            
        }
        if (global_error < global_error_ma){
            printf("Converged in %d iterations with error %e\n",k,global_error);
            break;
        }
        for(j=0;j<npoints;j++){
            x_old[j] = x[j];    
        }
    } 
}   