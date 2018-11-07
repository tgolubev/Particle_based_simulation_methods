#ifndef EXTPOTENTIAL_H
#define EXTPOTENTIAL_H

#include "vec2.h"

//! Contains properties and functions related to the external electric potential
class ExtPotential{
 public:
    vec2 position;
    double stdev;
    double potentialMean;
    double max;       //max of Gaussian
    double twostdevSqrd;

    ExtPotential() { }
    void setPosition(double x, double y);
    void setStdev(double input_stdev);
    void setMax(double inputmax);

    vec2 getForcefromPotential(vec2 atomPositions, double cutoffSqrd);
};


#endif // EXTPOTENTIAL_H
