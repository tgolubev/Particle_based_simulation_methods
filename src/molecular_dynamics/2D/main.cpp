//This is an object oriented molecular dynamics code which relies on many classes. This file allows to control most of the
//program parameters and choose which simulation approach (NVT, NVE, or gradual heating) to use. Certain blocks need
//to be commented/uncommented to use different approaches. This file can take input arguments, but they are not required.

//Outputs are a movie.xyz file with atom positions which can allow to create animations in i.e Ovito and statistics.txt
//which calculates system energetics, temperature, and diffusion coefficient. The frequency of output can be varied as
//desired. Of course, less frequent output, makes the program run faster.

//Author: Timofey Golubev based on the code found at https://github.com/andeplane/molecular-dynamics-fys3150.


#include "math/random.h"
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

using namespace std;
using namespace chrono;


int main()
{
    int system_select;
    int num_particles;
    bool input_parameters;
    vec2 numberOfUnitCellsEachDimension(6,6);
    double initialTemperature = 10; //in K
    double currentTemperature;
    double latticeConstant =10.256;  //in angstroms  //FOR 2D seems need larger latticeConstant to prevent blowup
    double sigma = 3.4; //atom/particle diameter in Angstrom for LJ potential
    double epsilon = 1.0318e-2; // epsilon from LJ in eV
    double side_length;
    double variance;
    double mass = 6.63352088e-26; // mass in kg
    bool input_variance;
    double total_dt_time= 0.0;

    double N_time_steps = 15000; //number of time steps

    //for NVT ensemble
    double N_steps = 100;  //number of steps over which to gradually rescale velocities: has to be large enough to prevent instability

    //input options interface
    cout<<"Do you want to input parameters, overiding default values?"<<endl;
    cin >> input_parameters;
    cout <<"Select System: 1 Random positions or 2 fcc lattice" <<endl;
    cin >> system_select;

    //---------------------------------------------------------------------------------------------------------------------------------
    //System dependent parameters
    if(system_select ==1){
        cout<<"Number of particles" <<endl;
        cin >>num_particles;
        cout<<"Side length of simulation cube" <<endl;
        cin >>side_length;
        if(input_parameters ==1){
            cout <<"particle mass (in kg)" <<endl;
            cin >> mass;
            cout <<"Initial temperature (in K)" <<endl;
            cin >> initialTemperature;
            cout<<"Sigma for LJ (atom diameter)" <<endl;
            cin >> sigma;
            cout<<"Epsilon for LJ (in eV)" <<endl;
            cin >> epsilon;
            cout<<"Would you like to input variance for Maxwellian velocity distr. (default = kT/m)? 1 Yes, 0 No" <<endl;
            cin >> input_variance;
            if(input_variance ==1){
                cout<<"Variance"<<endl;
                cin>> variance;
            }
        }
    }

    if(system_select == 2 && input_parameters ==1){
        cout <<"Number of Unit cells in each dimension" <<endl;
        cin>> numberOfUnitCellsEachDimension[0];
        numberOfUnitCellsEachDimension[1]  = numberOfUnitCellsEachDimension[0];
        cout <<"Initial temperature (in K)" <<endl;
       // cin >> intialTemperature;
        //initialTemperature = UnitConverter::temperatureFromSI(intialTemperature);
        cout <<"lattice constant (in Angstroms)" <<endl;
        cin >>latticeConstant;  
        cout <<"atom mass (in kg)" <<endl;
        cin >> mass;
        cout<<"Sigma for LJ (atom diameter)" <<endl;
        cin >> sigma;
        cout<<"Epsilon for LJ (in eV)" <<endl;
        cin >> epsilon;
        cout<<"Would you like to input variance for Maxwellian velocity distr. (default = kT/m)? 1 Yes, 0 No" <<endl;
        cin >> input_variance;
        if(input_variance ==1)
            cout<<"Variance. Note: k_b = 1, mass in amu, T "<<endl;
            cin>> variance;
    }

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
    if(system_select ==2) system.createFCCLattice(numberOfUnitCellsEachDimension, latticeConstant, initialTemperature, variance, input_variance, mass);
    else system.createRandomPositions(num_particles, side_length, initialTemperature, variance, input_variance, mass);
    system.potential().setEpsilon(1.0); //if don't set to 1.0, must modify LJ eqn.
    system.potential().setSigma(UnitConverter::lengthFromAngstroms(sigma));      //i.e. LJ atom/particle diameter,

    system.m_sample_freq=10; //statistics sampler freq.


    system.removeTotalMomentum();

    StatisticsSampler statisticsSampler;
    IO movie("movie.xyz"); // To write the positions to file, create IO object called "movie"
    movie.saveState(system, statisticsSampler);  //save the initial particle positions to file. pass statisticsSampler object too, so can use density function...

    ofstream velocities;
     velocities.open ("initial_velocities.txt");
     for(Atom *atom : system.atoms()) {
            velocities <<UnitConverter::velocityToSI(atom->velocity.x()) <<" "<<
            UnitConverter::velocityToSI(atom->velocity.y()) <<" "<< "\n";
     }
     velocities.close();


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
        if(timestep % 1000 ==0){
          //save atom coordinates only periodically to save CPU and file size
          movie.saveState(system, statisticsSampler);  //calls saveState fnc in io.cpp which saves the state to the movie.xyz file
        }
    }

    //timing
    high_resolution_clock::time_point finish2 = high_resolution_clock::now();
    duration<double> time2 = duration_cast<duration<double>>(finish2-start2);
    cout << "Total CPU time = " << time2.count() << endl;

    cout<<"Total dt time = " <<total_dt_time <<endl;

    movie.close();

    return 0;
}
