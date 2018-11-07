//------------------------------------------------------------------------------------------------------
// Definition of VelocityVerlet class which controls the integration of
// the equations of motion in order to time-step the simulation.
//
// By Timofey Golubev
//------------------------------------------------------------------------------------------------------

#ifndef VELOCITYVERLET_H
#define VELOCITYVERLET_H

class VelocityVerlet
{
private:
    bool m_firstStep = true;
public:
    VelocityVerlet() {}
    void integrate(class System &system, double dt);
};
#endif
