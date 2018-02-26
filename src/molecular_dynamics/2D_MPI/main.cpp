//This is an object oriented molecular dynamics code which relies on many classes. This file allows to control most of the
//program parameters and choose which simulation approach (NVT, NVE) to use. Certain blocks need
//to be commented/uncommented to use different approaches.

//Outputs are a movie.xyz file with atom positions which can allow to create animations in i.e Ovito and statistics.txt
//which calculates system energetics, temperature, and diffusion coefficient. The frequency of output can be varied as
//desired. Of course, less frequent output, makes the program run faster.

//Author: Timofey Golubev based on the code found at https://github.com/andeplane/molecular-dynamics-fys3150.

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
#include <stdio.h> //for printf
#include <mpi.h>  //get the mpi functions from here


using namespace std;
using namespace chrono;

int main(int argc, char **argv)
{
   int my_id, nprocs;
    
    MPI_Init(&argc, &argv);  //initialize MPI environment: takes in arguments from the mpirun command
    //everything past this point is executed on EACH processor in parallel

    //each processor finds out what it's ID (aka rank) is and how many processors there are
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);   //find ID
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);  //find # of processors

    double initialTemperature = 600.;//in K
    double currentTemperature;
    double latticeConstant =3.8;    //in angstroms
    double sigma = 3.4;             //atom/particle diameter in Angstrom for LJ potential
    double epsilon = 1.0318e-2;     // epsilon from LJ in eV
    double side_length;

    double total_dt_time= 0.0;

    //all variables will be defined in EACH processor
    vec2 Total_systemSize(90,30); //since using SC lattice--> just gives # of atoms in each dimension
    vec2 subsystemSize;
    subsystemSize[0] = Total_systemSize[0]/(nprocs-1);  //1D parallelization along x
    subsystemSize[1] = Total_systemSize[1]; //MUST CHANGE THIS TO /nprocs when do 2D parallelization

    vec2 subsystemOrigin; //bottom left corner of each subsystem, this will be defined in each processor seperately

    int StatSample_freq = 10;

    int N_time_steps = 10000; //number of time steps

    //for NVT ensemble
    int N_steps = 1;  //number of steps over which to gradually rescale velocities: has to be large enough to prevent instability
    //using the more basic thermostat of sharp rescaling every 100 steps

    //-----------------------------------------------------------------------------------------------------------------------------------------------
    int *NumInBox_array;  //declare
    int NumInBox = 0;
    double * positions;  //stores positions of each domain
    double * Allpositions;
    MPI_Datatype stype;



    //Initialize MD units
    if (my_id ==0){
        NumInBox_array = new int[nprocs]; //allocate memery only in proc 1, array for recieving number of particles in each procs domain
    }

    UnitConverter::initializeMDUnits(sigma, epsilon);
    //double mass = UnitConverter::massFromSI(6.63352088e-26); // mass in kg
    initialTemperature = UnitConverter::temperatureFromSI(initialTemperature);
    latticeConstant = UnitConverter::lengthFromAngstroms(latticeConstant);
    double dt = UnitConverter::timeFromSI(2e-15); //

    cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
    cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
    cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
    cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K" << endl;
    cout <<"Epsilon is " << epsilon <<" eV and sigma is " <<sigma <<" in Angstrom" << endl;
    cout << "mass in LJ units" << mass <<endl;
    cout << "timestep in LJ units" <<dt <<endl;

    System system;
    if(my_id !=0){
        system.createSCLattice(Total_systemSize,subsystemSize, latticeConstant, initialTemperature, mass, subsystemOrigin);
        system.potential().setEpsilon(1.0); //if don't set to 1.0, must modify LJ eqn.
        system.potential().setSigma(UnitConverter::lengthFromAngstroms(sigma));      //i.e. LJ atom/particle diameter,
        system.m_sample_freq=100; //statistics sampler freq.
        system.removeTotalMomentum();
    }

    //create_MPI_ATOM();

    StatisticsSampler statisticsSampler;

    FILE *movie;
    //IO movie("movie.xyz");

    //record left half of system

    //IO movie("movie.xyz"); // To write the positionitions to file, create IO object called "movie"

    /*

    if( my_id == 0 ) {
        // movie.saveState(system, statisticsSampler);  //save the initial particle positionitions to file. pass statisticsSampler object too, so can use density function...
    }
    //record right half of system
    IO movie2("movie_proc1.xyz");

    if(my_id == 2) {
        movie2.saveState(system, statisticsSampler);  //save the initial particle positionitions to file. pass statisticsSampler object too, so can use density function...
    }
    */

    if(my_id != 0){
        cout << setw(20) << "Timestep" <<
                setw(20) << "Time" <<
                setw(20) << "Temperature(not K!)" <<
                setw(20) << "KineticEnergy" <<

                setw(20) << "PotentialEnergy" <<
                setw(20) << "TotalEnergy" << endl;
    }else{
        movie = fopen("movie.xyz","w+");
    }

    high_resolution_clock::time_point start2 = high_resolution_clock::now();

    for(int timestep=0; timestep< N_time_steps; timestep++) {
        if(my_id != 0){
            //cout <<"line 140" <<endl;

            system.step(dt);  //only do timestepping for non-root procs

            //cout <<"timestep" <<timestep <<endl;

            //doesn't get here...

            //use sampler to calculate system parameters periodically
            if(timestep % system.m_sample_freq ==0){
                statisticsSampler.sample(system);
            }

            //periodically rescale Velocities to keep T constant (NVT ensemble)
            if(timestep % 100 == 0){
                statisticsSampler.sampleTemperature(system); //sample temperature, so can rescale as often as we want
                currentTemperature = statisticsSampler.temperature();  //this gets the value of temperature member variable
                //Note: initial temperature is the desired temperature here
                system.rescaleVelocities(statisticsSampler, currentTemperature, initialTemperature, N_steps);
            }

            if( timestep % 100 == 0 ) {
                // Print the timestep and system properties every 1000 timesteps
                cout << setw(20) << system.steps() <<
                        setw(20) << system.time() <<
                        setw(20) << statisticsSampler.temperature() <<
                        setw(20) << statisticsSampler.kineticEnergy() <<
                        setw(20) << "proc" << my_id <<  //this results in std::bad cast error
                        setw(20) << statisticsSampler.potentialEnergy() <<
                        setw(20) << statisticsSampler.totalEnergy() << endl;
            }
            NumInBox = system.num_atoms();

           //positions = getPositions(system);

            positions = new double[2*system.num_atoms()];
           //  cout <<"line 177" << system.atoms(0)->position[0] << endl;
            int i = 0, j=0;
            for(Atom *atom : system.atoms()) {
                positions[j] = system.atoms(i)->position[0];
                positions[j+1] = system.atoms(i)->position[1];
                i++;
                j+=2;
            }

           // cout <<"posotins" << positions[0] <<endl;

        }



        //gather
        //MPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)

        //gather number of particles in each processor 1st
       //right now proc 0 is BOTH the proc where data is gathered and still doing simulation...
        MPI_Gather(&NumInBox,1,MPI_INT,&NumInBox_array[0],1,MPI_INT,0,MPI_COMM_WORLD);
        int offset=0;
        int displs[nprocs];
        int rcounts[nprocs];
        int NumGlobal = 0;

        if(my_id == 0){
            //cout << "result of gather " << NumInBox_array[0] << " " << NumInBox_array[1] << " " << NumInBox_array[2] << endl;

            //make array of displacements for MPI_Gatherv
            for (int i=0;i<nprocs;i++){
                displs[i] = offset;
                offset += 2*NumInBox_array[i];
                rcounts[i] = 2*NumInBox_array[i];
                NumGlobal += NumInBox_array[i];
            }
            Allpositions = new double [2*NumGlobal];
        }




        //create MPI datatype to allow passing of varying sizes of arrays
         MPI_Type_vector(NumInBox,2,2,MPI_DOUBLE,&stype);  //make vector for atom positions, 2 elements/each atom
         MPI_Type_commit(&stype);

        //Gatherv allows to recieve varying sizes of data
        //MPI_Gatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,void *recvbuf, const int recvcounts[], const int displs[], MPI_Datatype recvtype, int root, MPI_Comm comm)


        MPI_Gatherv(&positions[0], 1, stype, Allpositions, rcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

       if(timestep % 100 ==0){
         if(my_id ==0){
             fprintf(movie,"%d \n", NumGlobal);
             fprintf(movie,"%11.3f \n",timestep);
             for(int i =0;i<NumGlobal;i++){
                 fprintf(movie, "H \t %11.3f \t %11.3f \n",Allpositions[2*i],Allpositions[2*i+1]);
             }

             delete [] Allpositions;
         }
       }






        //MPI_Gatherv(&positions[0], 2*NumInBox_array[my_id], MPI_INT, Allpositions, 2*NumInBox_array[my_id], displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

         //cout << "positions" << positions <<endl;

         /*
         if(my_id ==0) {
             for(int i = 0;i<2*NumInBox;i++){
                Allpositions[i] = positions[i];
             }
         }
         */

         //cout << "allpositions" << Allpositions[0] << " " <<Allpositions[1] <<endl;

        /*
         if(my_id !=0){
            MPI_Send(positions,1,stype,0,my_id,MPI_COMM_WORLD);
        }else{
            MPI_Recv(&Allpositions[displs[1]],rcounts[1],MPI_DOUBLE,1,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        */




    }//end of for loop

    if(my_id ==0) fclose(movie);


    //timing
    high_resolution_clock::time_point finish2 = high_resolution_clock::now();
    duration<double> time2 = duration_cast<duration<double>>(finish2-start2);
    cout << "Total CPU time for proc " << my_id <<" is " << time2.count() << endl;

    MPI_Finalize();

    return 0;
}
