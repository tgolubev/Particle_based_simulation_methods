#include "velocityverlet.h"
#include "system.h"
#include "atom.h"
#include "send_atoms.h"
#include <mpi.h>
#include "global.h"
#include "unitconverter.h"
#include <omp.h>
#include <stdio.h>


void VelocityVerlet::integrate(System &system, double dt) //passing by reference &system
{


    int nprocs, rank, namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);   //find ID
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);  //find # of processors
    MPI_Get_processor_name(processor_name, &namelen);

    //for OpenMP
    int iam = 0, np = 1;

    double half_dt = 0.5*dt;
    if(m_firstStep) {
        system.calculateForces();
        m_firstStep = false;
    }

    //std::cout <<"in VV" <<std::endl;

//#pragma omp parallel default(shared) private(iam, np)

       // np = omp_get_num_threads();
        //iam = omp_get_thread_num();
        //printf("This is thread %d out of %d from process %d out of %d on %s\n", iam, np, rank, nprocs, processor_name);



    //#pragma omp parallel for
    for(int i = 0;i<system.num_atoms();i++){
    //for(Atom *atom : system.atoms()) {  //compiler doesn't like this version when used together with "#pragma for"
        Atom* atom = system.atoms(i);

        //this operates on the vectors directly using vec2 class--> try doing with for loop and see if timing difference
        atom->velocity += atom->force*half_dt/mass;
        atom->position += atom->velocity*dt;  //NOTE: since v is computed 1st, this is actually x(t+dt) = x(t) + vt+0.5at^2 since v = v+0.5at
    }



    system.applyMirrorBCs(dt);
    //system.applyPeriodicBoundaryConditions();

    if(nprocs >1){
       send_atoms(&system);  //send and recieve atoms which have left their processor's domain
    }

    system.calculateForces(); // New positions, recompute forces

    //#pragma omp parallel for
    for(int i = 0;i<system.num_atoms();i++){
        Atom* atom = system.atoms(i);
    //for(Atom *atom : system.atoms()) {
        atom->velocity += atom->force*half_dt/mass;  //calculate new velocities
    }
}
