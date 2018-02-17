//This is an object oriented molecular dynamics code which relies on many classes. This file allows to control most of the
//program parameters and choose which simulation approach (NVT, NVE, or gradual heating) to use. Certain blocks need
//to be commented/uncommented to use different approaches. This file can take input arguments, but they are not required.

//Outputs are a movie.xyz file with atom positionitions which can allow to create animations in i.e Ovito and statistics.txt
//which calculates system energetics, temperature, and diffusion coefficient. The frequency of output can be varied as
//desired. Of course, less frequent output, makes the program run faster.

//Author: Timofey Golubev based on the code found at https://github.com/andeplane/molecular-dynamics-fys3150.


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
#include <stdio.h> //for printf
#include <mpi.h>  //get the mpi functions from here

using namespace std;
using namespace chrono;


//void decompositione1D_system(int system_length, int my_id,int nproc, int*  subdomain_start,int* subdomain_size) { //int* means pointer to int
//  *subsystem_length = system_length / nproc; //subsystem_length is a pointer, *subsystem_length changes the value of variable that pointer points to
//}


int main(int argc, char **argv)
{
   int my_id, nproc;
    
    MPI_Init(&argc, &argv);  //initialize MPI environment: takes in arguments from the mpirun command
    //everything past this point is executed on EACH processor in parallel

    //each processor finds out what it's ID (AKA rank) is and how many processors there are
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);   //find ID
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);  //find # of processors

    // Find your part of the domain--> determines subdomain location and size for each processor
    //int subsystem_length, subsystem_start;  //these will gain values after calling the decompositione_system fnc.
    //decompositione1D_system(system_length, my_id, nproc, &subsystem_start, &subdomain_size);
    //cout << "length of each subsystem" << subsystem_length <<endl;
    

    //all variables will be defined in EACH processor
    vec2 Total_systemSize(20,10); //TOTAL system dimensions
    vec2 subsystemSize; //this will be defined in each processor seperately
    vec2 subsystemOrigin; //bottom left corner of each subsystem, this will be defined in each processor seperately
    double initialTemperature = 10; //in K
    double currentTemperature;
    double latticeConstant =10.256;  //in angstroms  //FOR 2D seems need larger latticeConstant to prevent blowup
    double sigma = 3.4; //atom/particle diameter in Angstrom for LJ potential
    double epsilon = 1.0318e-2; // epsilon from LJ in eV
    double side_length;
    double mass = 6.63352088e-26; // mass in kg
    double total_dt_time= 0.0;

    int StatSample_freq = 10;

    int N_time_steps = 50000; //number of time steps

    //for NVT ensemble
    int N_steps = 100;  //number of steps over which to gradually rescale velocities: has to be large enough to prevent instability

    //-----------------------------------------------------------------------------------------------------------------------------------------------
    //Initialize MD units
     UnitConverter::initializeMDUnits(sigma, epsilon);
     initialTemperature = UnitConverter::temperatureFromSI(initialTemperature);
     latticeConstant = UnitConverter::lengthFromAngstroms(latticeConstant);
     double dt = UnitConverter::timeFromSI(2e-14); //

    cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
    cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
    cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
    cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K" << endl;
    cout <<"Epsilon is " << epsilon <<" eV and sigma is " <<sigma <<" in Angstrom" << endl;


    System system;
    system.createSCLattice(Total_systemSize,subsystemSize, latticeConstant, initialTemperature, mass, subsystemOrigin);
    //system.createFCCLattice(numberOfUnitCellsEachDimension, latticeConstant, initialTemperature, mass);
    //system.createRandompositionitions(num_particles, side_length, initialTemperature, mass);
    system.potential().setEpsilon(1.0); //if don't set to 1.0, must modify LJ eqn.
    system.potential().setSigma(UnitConverter::lengthFromAngstroms(sigma));      //i.e. LJ atom/particle diameter,
    system.m_sample_freq=10; //statistics sampler freq.
    system.removeTotalMomentum();

    create_MPI_ATOM();

    StatisticsSampler statisticsSampler;

    //record left half of system
    IO movie("movie.xyz"); // To write the positionitions to file, create IO object called "movie"
    IO movie2("movie2.xyz");
    if( my_id == 0 ) {
        movie.saveState(system, statisticsSampler);  //save the initial particle positionitions to file. pass statisticsSampler object too, so can use density function...
    }
    //record right half of system
    if(my_id == 1) {
        movie.saveState(system, statisticsSampler);  //save the initial particle positionitions to file. pass statisticsSampler object too, so can use density function...
    }


    //         ofstream velocities;
    //         velocities.open ("initial_velocities.txt");
    //         for(Atom *atom : system.atoms()) {
    //                 velocities <<UnitConverter::velocityToSI(atom->velocity.x()) <<" "<<
    //                 UnitConverter::velocityToSI(atom->velocity.y()) <<" "<< "\n";
    //         }
    //         velocities.close();


    cout << setw(20) << "Timestep" <<
            setw(20) << "Time" <<
            setw(20) << "Temperature(not K!)" <<
            setw(20) << "KineticEnergy" <<
            setw(20) << "PotentialEnergy" <<
            setw(20) << "TotalEnergy" << endl;

    high_resolution_clock::time_point start2 = high_resolution_clock::now();

    for(int timestep=0; timestep< N_time_steps; timestep++) {
        high_resolution_clock::time_point startdt = high_resolution_clock::now();
        system.step(dt);
        high_resolution_clock::time_point finishdt = high_resolution_clock::now();
        duration<double> time_dt = duration_cast<duration<double>>(finishdt-startdt);
        // cout << "dt time" << time_dt.count() <<endl;
        total_dt_time += time_dt.count();



        //Uncomment the below block to implement gradual heating of the system.
        /*
            //heat system gradually
            statisticsSampler.sampleKineticEnergy(system);      //can't sample temperature w/o sampling KE!
            statisticsSampler.sampleTemperature(system);  //sample temp at every timestep
            system.increaseTemperature(statisticsSampler, UnitConverter::temperatureFromSI(0.0001));  //Increase T by this increment (in K)
            */

        //use sampler to calculate system parameters and save to statistics.txt file
        if(timestep % system.m_sample_freq ==0){
            //to save CPU, can choose to sample only periodically
            statisticsSampler.sample(system);
        }

        //Uncoment the below block to use NVT ensemble.

        /*

            //periodically rescale Velocities to keep T constant (NVT ensemble)
            //if(timestep % 100 == 0){
                //rescale the velocities at every time step
            //sample temperature, so can rescale as often as we want
                statisticsSampler.sampleTemperature(system);
                currentTemperature = statisticsSampler.temperature();  //this gets the value of temperature member variable
            //Note: initial temperature is the desired temperature here
            system.rescaleVelocities(statisticsSampler, currentTemperature, initialTemperature, N_steps);
            //}
            */


        if( timestep % 1000 == 0 ) {
            // Print the timestep and system properties every 1000 timesteps
            cout << setw(20) << system.steps() <<
                    setw(20) << system.time() <<
                    setw(20) << statisticsSampler.temperature() <<
                    setw(20) << statisticsSampler.kineticEnergy() <<
                    setw(20) << statisticsSampler.potentialEnergy() <<
                    setw(20) << statisticsSampler.totalEnergy() << endl;
        }


    }

    //timing
    high_resolution_clock::time_point finish2 = high_resolution_clock::now();
    duration<double> time2 = duration_cast<duration<double>>(finish2-start2);
    cout << "Total CPU time = " << time2.count() << endl;

    cout<<"Total dt time = " <<total_dt_time <<endl;


    MPI_Finalize();



    return 0;
}
