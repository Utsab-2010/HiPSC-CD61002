#include <omp.h> 
#include <stdio.h> 
#include <stdlib.h> 
  
int main(int argc, char* argv[]) 
{
int i,id;
int a[10]; 
omp_set_num_threads(3);
#pragma omp parallel private(i,id) shared(a)
{
for (i=0;i<10;i++)
{
id=omp_get_thread_num();
a[i] = i+1; // updating private a
printf("a[%d] is %d from thread: %d\n",i,a[i],id); 
}}}
