//-------------------------------------------------------------------------------------------------------------
/*                   Parallelized 2D Molecular Dynamics

          Author: Timofey Golubev with debugging help of Xukun Xiang.

 Acknowledgements:
 Initial inspiration (with permission) for the molecular dynamics code is from
 the code found at: https://github.com/andeplane/molecular-dynamics-fys3150v.

 Description:

 This is a parallelized (using MPI) object-oriented molecular dynamics code.
 Parallelization is done using domain decomposition in 1D. The code can simulate
 the effects on atomic-like particles due to a moving external electric potential.
 This code is the "Model" portion of a group project on a parallelized interactive
 real time simulation of molecular dynamics using the Lennard Jones potential,
 meant for use on a Raspberry Pi supercomputer with an XBox One controller available at:
 https://github.com/MrHelloBye/Parallel-2D-LJMD

 Technical details:
 This file allows to control most of the program parameters.
 Output is a movie.xyz file with atom positions which can allow to create animations in
 i.e Ovito. Additional system properties such as energies, temperature, pressure, density,
 diffusion coefficient can be output if desired by using StatisticsSampler class. Output
 can be made in SI of Lennard Jones units by using UnitConverter class.
 The frequency of output can be varied as desired. Of course, less frequent output, results
 in lower CPU time.

 For a description of the molecular dynamics algorithms, see my report at:
 https://github.com/tgolubev/Computational-Physics-PHY905MSU/blob/master/Project4/molecular_dynamics/Report/Project4.pdf

*/
//-------------------------------------------------------------------------------------------------------------

#include "global.h"
#include "random.h"
#include "lennardjones.h"
#include "velocityverlet.h"
#include "system.h"
#include "statisticssampler.h"
#include "atom.h"
#include "io.h"
#include "unitconverter.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <time.h>
#include <stdio.h>
#include <mpi.h>

using namespace std;
using namespace chrono;

int main(int argc, char **argv)
{
    int my_id, nprocs;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // General simulation parameters
    double initialTemperature = 600.; //in K
    double currentTemperature;
    double latticeConstant =3.8;      // in Angstrom
    double sigma = 3.4;               //atom/particle diameter in Angstrom for LJ potential
    double epsilon = 1.0318e-2;       // epsilon from LJ in eV
    double side_length;
    double total_dt_time= 0.0;

    vec2 Total_systemSize(90, 30); // since using SC lattice--> just gives # of atoms in each dimension
    vec2 subsystemSize;
    subsystemSize[0] = Total_systemSize[0]/(nprocs-1);  //1D domain decomposition along x
    subsystemSize[1] = Total_systemSize[1];
    vec2 subsystemOrigin; //bottom left corner of each subsystem, this will be defined in each processor seperately

    int StatSample_freq = 10; // statistics sample frequency
    int N_time_steps = 10000; // number of time steps

    //for NVT ensemble
    int N_rescale_steps = 1;  //number of steps over which to gradually rescale velocities: has to be large enough to prevent instability

    //-----------------------------------------------------------------------------------------------------------------------------------------------
    // Variables for keeping track of atoms in each domain
    int *NumInBox_array;
    int NumInBox = 0;
    double *positions;  //stores positions of each domain
    double *Allpositions;
    MPI_Datatype stype;

    //Initialize MD units
    UnitConverter::initializeMDUnits(sigma, epsilon);
    initialTemperature = UnitConverter::temperatureFromSI(initialTemperature);
    latticeConstant = UnitConverter::lengthFromAngstroms(latticeConstant);
    double dt = UnitConverter::timeFromSI(2e-15);  // argument is in seconds

    cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
    cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
    cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
    cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K" << endl;
    cout <<"Epsilon is " << epsilon <<" eV and sigma is " <<sigma <<" in Angstrom" << endl;
    cout << "mass in LJ units" << mass <<endl;
    cout << "timestep in LJ units" <<dt <<endl;

    System system;
    if (my_id != 0) {
        system.createSCLattice(Total_systemSize,subsystemSize, latticeConstant, initialTemperature, mass, subsystemOrigin);
        system.potential().setEpsilon(1.0);                                          //if don't set to 1.0, must modify LJ eqn.
        system.potential().setSigma(UnitConverter::lengthFromAngstroms(sigma));      //i.e. LJ atom/particle diameter,
        system.m_sample_freq=100;                                                    //statistics sampler freq.
        system.removeTotalMomentum();
    }

    //setup external potential--> object extPotential is declare as public member of system
    system.extPotential.position.set(5.,system.systemSize(1)/2);
    system.extPotential.setMax(200.);  //200 is stable when use 1sigma for stdev
    system.extPotential.setStdev(UnitConverter::lengthFromAngstroms(1*sigma));

    StatisticsSampler statisticsSampler;

    FILE *movie;
    //IO movie("movie.xyz");

    if (my_id != 0) {
        cout << setw(20) << "Timestep" <<
                setw(20) << "Time" <<
                setw(20) << "Temperature(not K!)" <<
                setw(20) << "KineticEnergy" <<
                setw(20) << "Processor Rank" <<
                setw(20) << "PotentialEnergy" <<
                setw(20) << "TotalEnergy" << endl;
    } else {
        movie = fopen("movie.xyz","w+");
        NumInBox_array = new int[nprocs]; //allocate memery only in proc 1, array for recieving number of particles in each procs domain
    }

    high_resolution_clock::time_point start2 = high_resolution_clock::now();

    for (int timestep = 0; timestep < N_time_steps; timestep++) {

        //move around the ext Potential to show example
        // in the full project, the positions are recieved from the X-box controller
        if(timestep % 10 == 0){
               if(system.extPotential.position[0] < 85) system.extPotential.position[0] += 0.1;
               //if(my_id == 0) cout << "potential position" << system.extPotential.position[0] <<std::endl;
        }

        if (my_id != 0) {

            system.step(dt);  //only do timestepping for non-root procs

            //use sampler to calculate system parameters periodically
            if (timestep % system.m_sample_freq ==0)
                statisticsSampler.sample(system);

            //periodically rescale Velocities to keep T constant (NVT ensemble)
            if (timestep % 100 == 0) {
                statisticsSampler.sampleTemperature(system); //sample temperature, so can rescale as often as we want
                currentTemperature = statisticsSampler.temperature();  //this gets the value of temperature member variable
                //Note: initial temperature is the desired temperature here
                system.rescaleVelocities(statisticsSampler, currentTemperature, initialTemperature, N_rescale_steps);
            }

            if (timestep % 100 == 0 ) {
                // Print the timestep and system properties every 100 timesteps
                cout << setw(20) << system.steps() <<
                        setw(20) << system.time() <<
                        setw(20) << statisticsSampler.temperature() <<
                        setw(20) << statisticsSampler.kineticEnergy() <<
                        setw(20) << "proc" << my_id <<
                        setw(20) << statisticsSampler.potentialEnergy() <<
                        setw(20) << statisticsSampler.totalEnergy() << endl;
            }
            NumInBox = system.num_atoms();

            positions = new double[2*system.num_atoms()];
            int i = 0, j=0;
            for (Atom *atom : system.atoms()) {
                positions[j] = system.atoms(i)->position[0];
                positions[j+1] = system.atoms(i)->position[1];
                i++;
                j+=2;
            }
        }

        //gather number of particles in each processor 1st
        //MPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
        MPI_Gather(&NumInBox, 1, MPI_INT, &NumInBox_array[0], 1, MPI_INT, 0, MPI_COMM_WORLD);
        int offset=0;
        int displs[nprocs];
        int rcounts[nprocs];
        int NumGlobal = 0;

        if (my_id == 0) {
            //cout << "result of gather " << NumInBox_array[0] << " " << NumInBox_array[1] << " " << NumInBox_array[2] << endl;

            //make array of displacements for MPI_Gatherv
            for (int i = 0; i < nprocs; i++) {
                displs[i] = offset;
                offset += 2*NumInBox_array[i];
                rcounts[i] = 2*NumInBox_array[i];
                NumGlobal += NumInBox_array[i];
            }
            Allpositions = new double [2*NumGlobal];
        }

        //create MPI datatype to allow passing of varying sizes of arrays
        MPI_Type_vector(NumInBox, 2, 2, MPI_DOUBLE, &stype);  //make vector for atom positions, 2 elements/each atom
        MPI_Type_commit(&stype);

        //Gatherv allows to recieve varying sizes of data
        //MPI_Gatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,void *recvbuf, const int recvcounts[], const int displs[], MPI_Datatype recvtype, int root, MPI_Comm comm)
        MPI_Gatherv(&positions[0], 1, stype, Allpositions, rcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if (timestep % 100 == 0) {
            if (my_id == 0) {
                fprintf(movie,"%d \n", NumGlobal);
                fprintf(movie,"%11.3f \n",timestep);  //comment line

                int num=0;
                for (int i = 0; i < nprocs; i++) {
                    for (int j = 0; j < NumInBox_array[i]; j++) {
                        fprintf(movie, "H \t %11.3f \t %11.3f \t %d \n",Allpositions[num+2*j],Allpositions[num+2*j+1],i);  //this i will output the proc which each atom is a part of
                    }
                    num += 2*NumInBox_array[i];  //keeps track of start locations for elements belonging to next proc
                }

                delete [] Allpositions;

                // print position of external potential
                double x_pos = system.extPotential.position[0];
                double y_pos = system.extPotential.position[1];
                //std::cout << system.extPotential.position[0] << system.extPotential.position[1];
                fprintf(movie, "Ar \t %11.3f \t %11.3f \n", x_pos, y_pos);
            }
        }
    }

    if (my_id == 0) fclose(movie);

    //timing
    high_resolution_clock::time_point finish2 = high_resolution_clock::now();
    duration<double> time2 = duration_cast<duration<double>>(finish2-start2);
    cout << "Total CPU time for proc " << my_id <<" is " << time2.count() << endl;

    MPI_Finalize();

    return 0;
}
