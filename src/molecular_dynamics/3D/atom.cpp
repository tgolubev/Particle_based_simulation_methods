#include "atom.h"
#include "math/random.h"
#include <cmath>

Atom::Atom(double mass) :
    m_mass(mass)
{

}

void Atom::setInitialPosition(double x, double y, double z)
{
    m_initial_position.set(x,y,z);
    position.set(x,y,z);
}

void Atom::resetForce()
{
    force.zeros();
}

void Atom::resetVelocityMaxwellian(double temperature, double variance, bool input_variance)
{
    // Resetting the velocity according to a Maxwell-Boltzmann distribution (see http://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution )
    double boltzmannConstant = 1.0; // In these units, the boltzmann constant equals 1
    double standardDeviation;
    if(input_variance ==1) standardDeviation = sqrt(variance);
    else standardDeviation = sqrt(boltzmannConstant*temperature/m_mass);
    velocity.randomGaussian(0, standardDeviation);  //note: randomGaussian is defined in vec3.cpp. arguments(mean, stdev), so this is Gaussian centered about 0
}
