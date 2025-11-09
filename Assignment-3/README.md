To run the omp script:
- 1. Compile using ```gcc bicgstab_omp -o bicgstab_omp -lm -fopenmp```
- 2. run it using `bicgstab_omp 100_nodes 4`. Instead of "100_nodes" you can use the folder path having hte Kmat.txt and Fvec.txt files.
