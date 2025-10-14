#include <omp.h> 
#include <stdio.h> 
#include <stdlib.h> 
  
int main(int argc, char* argv[]) 
{
int i,id,a; 
#pragma omp parallel for private(a,id) ordered
for (i=0;i<8;i++)
{
#pragma omp ordered
id=omp_get_thread_num();
a = id+1; // updating private a
printf("a is %d from thread: %d\n",a,id); 
} }                                          
