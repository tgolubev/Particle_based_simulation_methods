#include "velocityverlet.h"
#include "system.h"
#include "atom.h"
#include "send_atoms.h"


void VelocityVerlet::integrate(System &system, double dt) //passing by reference &system, passes using the address but doesn't create a variable (pointer) whose value = that address
{


    double half_dt = 0.5*dt;
    if(m_firstStep) {
        system.calculateForces();
        m_firstStep = false;
    }

    for(Atom *atom : system.atoms()) {
        //this operates on the vectors directly using vec3 class
        atom->velocity += atom->force*half_dt/atom->mass();
        atom->position += atom->velocity*dt;  //NOTE: since v is computed 1st, this is actually v = vt+0.5at^2 since v = v+0.5at
    }

    system.applyPeriodicBoundaryConditions();

    //send_atoms accepts a pointer, so declare a pointer(variable whose value = the memory address)
    System  * pt_system = &system; //assign it the address of system

    send_atoms(&system);  //send and recieve atoms which have left their processor's domain

    system.calculateForces(); // New positions, recompute forces

    for(Atom *atom : system.atoms()) {
        atom->velocity += atom->force*half_dt/atom->mass();  //calculate new velocities
    }
}
