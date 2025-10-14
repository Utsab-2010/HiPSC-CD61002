#include <omp.h> 
#include <stdio.h> 
#include <stdlib.h> 
  
int main(int argc, char* argv[]) 
{
int id; 
int x = 5;  // shared memory x
#pragma omp parallel private(x,id)
{
id=omp_get_thread_num();
x = x+id; // updating private x
printf("private: x is %d from thread: %d\n",x,id); 
} // writing out private x by each thread
printf("after: x is %d\n",x);
} 
