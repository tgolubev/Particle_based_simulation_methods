//------------------------------------------------------------------------------------------------------
// Definition of StatisticsSampler class which contains all functions for calculating physical
// quantities about the system such as temperature, pressure, etc.
//
//
// Small modifications by Timofey Golubev
//------------------------------------------------------------------------------------------------------

#ifndef STATISTICSSAMPLER_H
#define STATISTICSSAMPLER_H
#include <fstream>
#include "vec2.h"

class System; // Promise the compiler that this is a class even though we haven't included system.h here

class StatisticsSampler
{
private:
    std::ofstream m_file;
    vec2 m_totalMomentum;
    double m_kineticEnergy = 0;
    double m_potentialEnergy = 0;
    double m_temperature = 0;
    double m_density = 0;
    double m_diffusion_coeff = 0;
    double m_totalMass = 0;
public:
    StatisticsSampler();
    void saveToFile(System &system);
    void sample(System &system);
    void sampleKineticEnergy(System &system);
    void samplePotentialEnergy(System &system);
    void sampleTemperature(System &system);
    void sampleDensity(System &system);
    void sampleDiffusionCoeff(System &system);
    void sampleMomentum(System &system);
    vec2 totalMomentum(){return m_totalMomentum;}
    double kineticEnergy() { return m_kineticEnergy; }
    double potentialEnergy() { return m_potentialEnergy; }
    double totalEnergy() { return m_kineticEnergy+m_potentialEnergy; }
    double temperature() { return m_temperature; }
    double density() { return m_density; }
};
#endif
