#include "system.h"
#include "statisticssampler.h"
#include "lennardjones.h"
#include <iostream>
#include <iomanip>
#include "unitconverter.h"

using std::ofstream; using std::cout; using std::endl;

StatisticsSampler::StatisticsSampler()
{

}

void StatisticsSampler::saveToFile(System &system)
{
    // Save the statistical properties for each timestep for plotting etc.
    // First, open the file if it's not open already

    if(!m_file.good()) {   //m_file is an ofstream
        m_file.open("statistics.txt", ofstream::out);
        // If it's still not open, something bad happened...
        if(!m_file.good()) {
            cout << "Error, could not open statistics.txt" << endl;
            exit(1);
        }
    }

    // Print out values here
    //Using SI units

   /* m_file << std::setw(15) <<UnitConverter::timeToSI(system.time()) << //UnitConverter:: allows to access fncs in that class
              //std::setw(10) << m_totalMomentum <<
              std::setw(15) << UnitConverter::energyToSI(m_kineticEnergy)<<
              std::setw(15) << UnitConverter::energyToSI(m_potentialEnergy) <<
              std::setw(15) << UnitConverter::temperatureToSI(m_temperature)<<
              std::setw(15) << UnitConverter::diffusionToSI(m_diffusion_coeff)<<
              std::setw(15) << UnitConverter::pressureToSI(m_pressure)<<endl;   //SI units are N/m^3 which is Pascal, 1atm = 1.01325*10^5 Pa
              */

    //Using MD units

    m_file << std::setw(10) <<system.time() <<
              //std::setw(10) << m_totalMomentum <<
              std::setw(10) << m_kineticEnergy<<
              std::setw(10) << m_potentialEnergy <<
              std::setw(12) << m_temperature<<
              std::setw(12) << m_diffusion_coeff<<
              std::setw(12) << m_density<<
              std::setw(12) << m_pressure_ideal <<
              std::setw(12) << m_pressure_ext <<
              std::setw(12) << m_pressure<<endl;
}


void StatisticsSampler::sample(System &system)
{
    // Here you should measure different kinds of statistical properties and save it to a file.
    //sampleMomentum(system);
    sampleKineticEnergy(system);
    samplePotentialEnergy(system);
    sampleTemperature(system);
    //sampleDiffusionCoeff(system);
    //sampleDensity(system);
    samplePressure(system);
    saveToFile(system);
}

void StatisticsSampler::sampleMomentum(System &system)
{
    m_totalMomentum.set(0,0,0);  //reset total-momentum to 0
    for(Atom *atom : system.atoms()) { //c++11 way of iterating through  entire vector or array
          m_totalMomentum += atom->mass()*atom->velocity;  //mass() returns value of m_mass (atom's mass)
     }
}

void StatisticsSampler::sampleKineticEnergy(System &system)
{
    m_kineticEnergy = 0.; // Remember to reset the value from the previous timestep
    for(Atom *atom : system.atoms()) {
        m_kineticEnergy += 0.5*atom->mass()*atom->velocity.lengthSquared();
    }
}

void StatisticsSampler::samplePotentialEnergy(System &system)
{
    m_potentialEnergy = system.potential().potentialEnergy();
}

void StatisticsSampler::sampleTemperature(System &system)
{
    //Reuse the kinetic energy that we already calculated
    m_temperature = (2./3.)*m_kineticEnergy/system.num_atoms();
    //cout <<"num_atoms= " << system.num_atoms() << endl;
}

void StatisticsSampler::sampleDensity(System &system)
{
  //reset quantities
    m_totalMass = 0;
  //total mass over volume
    for(Atom *atom : system.atoms()) {
        m_totalMass += atom->mass();

        m_density = m_totalMass/system.volume();  //already have this function
    }

}



void StatisticsSampler::samplePressure(System &system){             //ideal gas pressure: P = (n/3N)sum(mv_i^2)
    LennardJones m_potential;
    //m_pressure_ideal = 2*m_kineticEnergy*m_density/(3*system.num_atoms()); //NOTE: density must be sampled 1st to be able to get pressure
    m_pressure_ideal = (system.num_atoms()/system.volume())*m_temperature;  //NOTE: P_ideal = nT  , if rewrite the eqns: use KE = (3/2)NT
    //std::cout <<"pressure virial in statasamaple" <<system.m_potential.pressureVirial()<<std::endl;
    m_pressure_ext = (1/(3*system.volume()))*system.m_potential.pressureVirial(); //NOTE: m_potential is a LJ object IN system class pressure_virial is public variable so is accessible
    m_pressure = m_pressure_ideal + m_pressure_ext;

    //external pressure needs to be found in LJ since involves displacements
}


void StatisticsSampler::sumPressures(){
    //this is for averaging
    pressure_ideal_sum += m_pressure_ideal;
    pressure_ext_sum += m_pressure_ext;
    pressure_sum += m_pressure;
    //std::cout <<"pressure sum" <<pressure_sum <<std::endl;
    num_samples++;
}

void StatisticsSampler::sampleDiffusionCoeff(System &system)
{
    double displacements_sqrd_sum = 0.0;  //reset displacements sum
    for(Atom *atom : system.atoms()) {
        vec3 total_displacement;
        for(int j=0;j<3;j++){
            //takes into account displacement within 1 cell plus displacement due to crossing boundaries into neighboring image cells (PBCs)
            total_displacement[j] = (atom->position[j] - atom->initial_position(j)) + atom->num_bndry_crossings[j]*system.systemSize(j);
        }
        double total_displacement_sqrd = total_displacement.lengthSquared();
        displacements_sqrd_sum += total_displacement_sqrd;
        }
    m_diffusion_coeff = displacements_sqrd_sum/(6*system.time()*system.num_atoms());  //Einstein relation: D = (mean sqr displacement)/6t
}


