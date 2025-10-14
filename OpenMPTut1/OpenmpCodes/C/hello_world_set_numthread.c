#include <omp.h> 
#include <stdio.h> 
#include <stdlib.h> 
  
int main(int argc, char* argv[]) 
{
  omp_set_num_threads(4); 
    #pragma omp parallel 
    { 
       int id = omp_get_thread_num(); 
       printf("Hello World from thread id %d \n",id);
    } 
} 
