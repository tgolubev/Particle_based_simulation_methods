#define _USE_MATH_DEFINES
#include "unitconverter.h"
#include <iostream>
//#include <cmath>
#include <math.h>

using namespace std;

double UnitConverter::m0 = 0;
double UnitConverter::q0 = 0;
double UnitConverter::hbar0 = 0;
double UnitConverter::electricConstant0 = 0;
double UnitConverter::a0 = 0;
double UnitConverter::a0_angstrom = 0;
double UnitConverter::E0 = 0;
double UnitConverter::E0ev = 0;
double UnitConverter::kb = 0;
double UnitConverter::t0 = 0;
double UnitConverter::F0 = 0;
double UnitConverter::T0 = 0;
double UnitConverter::P0 = 0;
double UnitConverter::v0 = 0;
double UnitConverter::visc0 = 0;
double UnitConverter::diff0 = 0;
std::string UnitConverter::currentUnits = "No units chosen.";
bool UnitConverter::initialized = false;


void UnitConverter::initializeMDUnits(double sigma, double epsilon) {
    UnitConverter::initialized = true;
    UnitConverter::currentUnits = "MD units";
    // Molecular Dynamics units

    // Fundamental units
    double m0 = 1.66053892e-27;         // 1 amu in SI [kg]
    //double L0 = 1e-10;
    double L0 = sigma*1.e-10;                // 1 particle diameter in SI [m]
    double kb = 1.3806488e-23;          // SI [J/K]
    double E0eV = epsilon;            //epsilon (from coeff. in front of LJ) in in eV
    double E0 = 1.60217657e-19*E0eV;    // was Ar value  convert epsilon to SI [J]

    UnitConverter::m0 = m0;
    UnitConverter::kb = kb;

    // Derived units
    UnitConverter::a0 = L0;
    UnitConverter::E0 = E0;
    UnitConverter::t0 = L0*sqrt(m0/E0);
    UnitConverter::v0 = L0/t0;
    UnitConverter::F0 = E0/a0;
    UnitConverter::T0 = E0/kb;  //epsilon/kb  --> epsilon is from LJ
    UnitConverter::P0 = E0/(a0*a0*a0);
    UnitConverter::visc0 = P0*t0;
    UnitConverter::diff0 = a0*a0/t0;
    UnitConverter::E0ev = E0eV;
    UnitConverter::hbar0 = INFINITY; // Not really used, give INF as a warning-ish?
    UnitConverter::q0 = INFINITY; // Not really used, give INF as a warning-ish?
    UnitConverter::electricConstant0 = INFINITY; // Not really used, give INF as a warning-ish?
}

void UnitConverter::makeSureInitialized() {
    if(!UnitConverter::initialized) std::cout <<"ERROR: Unit converter not initialized" << std::endl;//UnitConverter::initialize(MDUnits);
}

/*void UnitConverter::initialize(Units type) {
    //if(type == MDUnits) UnitConverter::initializeMDUnits(sigma);
}
*/

double UnitConverter::pressureToSI(double P) {UnitConverter::makeSureInitialized(); return UnitConverter::P0*P; }
double UnitConverter::pressureFromSI(double P) {UnitConverter::makeSureInitialized(); return P/UnitConverter::P0; }

double UnitConverter::temperatureToSI(double T) {UnitConverter::makeSureInitialized(); return UnitConverter::T0*T; }
double UnitConverter::temperatureFromSI(double T) {UnitConverter::makeSureInitialized(); return T/UnitConverter::T0; }

double UnitConverter::massToSI(double m) {UnitConverter::makeSureInitialized(); return UnitConverter::m0*m; }
double UnitConverter::massFromSI(double m) {UnitConverter::makeSureInitialized(); return m/UnitConverter::m0; }

double UnitConverter::lengthToSI(double L) {UnitConverter::makeSureInitialized(); return UnitConverter::a0*L; }
double UnitConverter::lengthFromSI(double L) {UnitConverter::makeSureInitialized(); return L/UnitConverter::a0; }

double UnitConverter::lengthToAngstroms(double L) {UnitConverter::makeSureInitialized(); return UnitConverter::a0*L*1e10; }
double UnitConverter::lengthFromAngstroms(double L) {UnitConverter::makeSureInitialized(); return L/(UnitConverter::a0*1e10); }


vec2 UnitConverter::lengthToSI(vec2 position)
{
    return vec2(UnitConverter::lengthToSI(position.x()), UnitConverter::lengthToSI(position.y()));
}

vec2 UnitConverter::lengthFromSI(vec2 position)
{
    return vec2(UnitConverter::lengthFromSI(position.x()), UnitConverter::lengthFromSI(position.y()));
}

vec2 UnitConverter::lengthToAngstroms(vec2 position)
{
    return vec2(UnitConverter::lengthToAngstroms(position.x()), UnitConverter::lengthToAngstroms(position.y()));
}

vec2 UnitConverter::lengthFromAngstroms(vec2 position)
{
    return vec2(UnitConverter::lengthFromAngstroms(position.x()), UnitConverter::lengthFromAngstroms(position.y()));
}

vec2 UnitConverter::velocityToSI(vec2 velocity)
{
    return vec2(UnitConverter::velocityToSI(velocity.x()), UnitConverter::velocityToSI(velocity.y()));
}

vec2 UnitConverter::velocityFromSI(vec2 velocity)
{
    return vec2(UnitConverter::velocityFromSI(velocity.x()), UnitConverter::velocityFromSI(velocity.y()));
}

double UnitConverter::forceToSI(double F) {UnitConverter::makeSureInitialized(); return UnitConverter::F0*F; }
double UnitConverter::forceFromSI(double F) {UnitConverter::makeSureInitialized(); return F/UnitConverter::F0; }

double UnitConverter::energyToSI(double E) {UnitConverter::makeSureInitialized(); return UnitConverter::E0*E; }  //this converts back correctly by multiplying by epsilon
double UnitConverter::energyFromSI(double E) {UnitConverter::makeSureInitialized(); return E/UnitConverter::E0; }

double UnitConverter::energyToEv(double E) {UnitConverter::makeSureInitialized(); return UnitConverter::E0ev*E; }
double UnitConverter::energyFromEv(double E) {UnitConverter::makeSureInitialized(); return E/UnitConverter::E0ev; }

double UnitConverter::degreesToRadians(double v) {UnitConverter::makeSureInitialized(); return M_PI/180*v; }
double UnitConverter::radiansToDegrees(double v) {UnitConverter::makeSureInitialized(); return 180/M_PI*v; }

double UnitConverter::timeToSI(double t) {UnitConverter::makeSureInitialized(); return UnitConverter::t0*t; }
double UnitConverter::timeFromSI(double t) {UnitConverter::makeSureInitialized(); return t/UnitConverter::t0; }

double UnitConverter::velocityToSI(double v) {UnitConverter::makeSureInitialized(); return v*UnitConverter::v0; }
double UnitConverter::velocityFromSI(double v) {UnitConverter::makeSureInitialized(); return v/UnitConverter::v0; }

double UnitConverter::diffusionToSI(double d) {UnitConverter::makeSureInitialized(); return d*UnitConverter::diff0; }
double UnitConverter::diffusionFromSI(double d) {UnitConverter::makeSureInitialized(); return d/UnitConverter::diff0; }
