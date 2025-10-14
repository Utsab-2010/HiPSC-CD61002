#include <omp.h> 
#include <stdio.h> 
#include <stdlib.h> 
  
int main(int argc, char* argv[]) 
{ 
    #pragma omp parallel num_threads(4)
    { 
       int id = omp_get_thread_num(); 
       printf("Hello World from thread id %d \n",id);
    }
    #pragma omp parallel
    {
       int id = omp_get_thread_num();
       printf(" Again Hello World from thread id %d \n",id);
    }
} 
