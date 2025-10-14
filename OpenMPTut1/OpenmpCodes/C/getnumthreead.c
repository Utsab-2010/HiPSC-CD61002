#include <omp.h> 
#include <stdio.h> 
#include <stdlib.h> 
  
int main(int argc, char* argv[]) 
{
    #pragma omp parallel 
    { 
       int id = omp_get_thread_num();
       int nthread = omp_get_num_threads(); 
       printf("This thread id %d from total %d threads \n",id,nthread);
    } 
} 
