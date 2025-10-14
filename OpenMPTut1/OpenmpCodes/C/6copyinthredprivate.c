#include <omp.h> 
#include <stdio.h> 
#include <stdlib.h> 
  
int main(int argc, char* argv[]) 
{
int id;
static int x; 
#pragma omp threadprivate(x)
#pragma omp parallel private(id)
{
id=omp_get_thread_num();
x = x+id; // updating threadprivate x
printf("private: x is %d from thread: %d\n",x,id); 
} 
#pragma omp parallel copyin(x)
    //copying copies masterthread's threadprivate
    //values to other threads
id=0;
#pragma omp parallel private(id) 
{
id=omp_get_thread_num();
x = x+10; // updating threadprivate x
printf("2nd parallel zone x is %d from thread: %d\n",x,id);
} // writing threadprivate variables 
} // second parallel loop 
