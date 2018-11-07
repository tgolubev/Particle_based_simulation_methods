# Parallel 2D Molecular Dynamics

## Description
This is a parallelized version of the 2D MD code using OpenMPI (Message Passing Interface).
* The parallelization is done via domain decomposition in 1D along the x direction. 
* Each processor (except the master) has a subsystem for which it performs all calculations. 
* When atoms move between domains, information is exchanged using Isend and Irecv MPI calls. 
* Mirror boundary conditions are applied for exterior system boundaries. 
* A cutoff radius where force calculations byond the radius are ignored, is used to significantly speed up calculations. 
* MPI Gather and Gatherv is used for the master to gather system information from the other processors.

## Compiling and running:

The included make file requires having OpenMPI installed.

To compile in Linux, just use `make`

To run in Linux: `mpirun -np [number of processors] MPI-md` 

Note: Processor 0 is used as the master for gathering data from other processors and output. It does not perform any MD calculations. # of subdomains = # of processors - 1 

**Note: For stability, each processor's domain must be larger than the cutoff radius. This means that for the current parameters in main.cpp, using <= 5 processors works best.**
