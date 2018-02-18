#include "system.h"
#include "velocityverlet.h"
#include "lennardjones.h"
#include "statisticssampler.h"
#include "unitconverter.h"
#include "random.h"
#include <map>
#include <mpi.h>
#include <algorithm> //need for sort etc.

System::System()   //constructor
{
  m_num_atoms = 0;  //make num_atoms = 0 upon initialization
}

System::~System()  //desctructor: remove all the atoms
{
    for(Atom *atom : m_atoms) {
        delete atom;
    }
    m_atoms.clear();    //vector containing all atoms objects
}

void System::applyPeriodicBoundaryConditions() {
    // Read here: http://en.wikipedia.org/wiki/Periodic_boundary_conditions#Practical_implementation:_continuity_and_the_minimum_image_convention
    //using version A: where fold the atoms back into the simulation cell. The cell has its lower left corner at origin of coord. system.

    for(Atom *atom : atoms()) {

        for(int j=0;j<2;j++){
            //fold atoms back into the box if they escape the box
            //Note: if atom in 1 step moves further than the neighboring image of the simulation cell, it is not fully brought back.
            if (atom->position[j] <  0.) {
                atom->position[j] += m_systemSize[j];
                atom->num_bndry_crossings[j] -= 1;  //crossing left or bottom boundary is counted as -1 crossing.
            }
            if (atom->position[j] >  m_systemSize[j]){
                atom->position[j] -= m_systemSize[j];
                atom->num_bndry_crossings[j] += 1;    //crossing right or top boundary is counted as +1 crossing
            }
        }

    }
}


// remove atoms that escape the SUB system (corresponding to current processor)
/*
void System::removeEscapedAtoms() {
    for(Atom *atom : atoms()) {
        if(atom->position[0]>m_subsystemSize[0])
            //m_atoms.erase(atom);
            delete atom;
    }
}
*/



void System::rescaleVelocities(StatisticsSampler &statisticsSampler, double currentTemperature, double desiredTemperature, int N_steps){
    //rescale velocities using equipartition theorem: v_desired = sqrt(T_desired/T_actual)*v_actual
    //double rescaling_factor = sqrt(desiredTemperature/statisticsSampler.temperature()); //sqrt(T_desired/T_actual)

    double rescaling_factor = sqrt(1+ (1/N_steps)*(desiredTemperature/currentTemperature-1)); //sqrt(1+(1/N_steps)(T_desired/T_current - 1))
    //Note: if N_steps is set to 1, then just gives the old scaling: sqrt(T_desired/T_actual)
    for(Atom *atom : atoms()) {
        atom->velocity *= rescaling_factor;  //a*=b means a = a*b
    }
    removeTotalMomentum();  //If don't do this, eventually will have drifting issue!
}

void System::removeTotalMomentum() {
    // Find the total momentum and remove momentum equally on each atom so the total momentum becomes zero.
    vec2 total_momentum;
    for(Atom *atom : atoms()) { //c++11 way of iterating through  entire vector or array
        total_momentum += atom->mass()*atom->velocity;  //mass() returns value of m_mass (atom's mass)
    }
    vec2 Mom_per_atom;   //2D momentum components
    Mom_per_atom = total_momentum/num_atoms(); //this is amount of abs(momentum) per atom that need to remove to get total to 0
    for(Atom *atom : atoms()) {
        for(int j=0;j<2;j++){
            //evenly modify components of velocity of each atom to yield total system momentum of 0
            atom->velocity[j] -= Mom_per_atom[j]/atom->mass();
        }
    }

    //test if total_momentum was rezeroed 0: (a unit test)
    /*
    total_momentum.print("Total Momentum before removeMomentum");
    total_momentum.set(0,0,0);  //reset total-momentum to 0
    for(Atom *atom : atoms()) { //c++11 way of iterating through  entire vector or array
          total_momentum += atom->mass()*atom->velocity;  //mass() returns value of m_mass (atom's mass)
     }
    total_momentum.print("Total Momentum after removeMomentum");   //print() is fnc in vec2 class, takes in a string input
    */
}

//note: currently upon initialization, each processor will create the ENTIRE global system. Then will send atoms and delete uneecesary atoms when send_atoms
//is called in velocityverlet.cpp. Could optimize to have only the subdomains created here, but not sure how much CPU this saves since just corresponds to 1 dt send commands.
void System::createSCLattice(vec2 Total_systemSize, vec2 subsystemSize, double latticeConstant, double temperature, double mass, vec2 subsystemOrigin) {

    double x = 0;
    //each processor finds out what it's ID (AKA rank) is and how many processors there are
    int nprocs, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);   //find ID
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);  //find # of processors

    x =rank*subsystemSize[0]*latticeConstant;  //initializes x at the right position for the processor's domain

    //std::cout <<"x" <<x <<std::endl;  //works correctly


    for(int i=0;i<subsystemSize[0];i++){  //< b/c i starts w/ 0
        //i.e. i = 0,1...N_x-1

        double y = 0;

        for(int j=0;j<subsystemSize[1];j++){

                Atom *atom = new Atom(UnitConverter::massFromSI(mass)); //uses mass in kg: mass is correct

                atom->setInitialPosition(x,y);
                atom->num_bndry_crossings.set(0.,0.);   //make sure initial # of bndry crossings is 0
                atom->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom);     //add element to vector m_atoms 1 element (atom object)
                y += latticeConstant;
        }
        x +=latticeConstant;
    }

    //this sets the TOTAL system size--> is used for PBCs which are at the outer boundaries
    setSystemSize(latticeConstant*Total_systemSize); //system size set by multiply vec2 # of unit cells by latticeConstant
    vec2 include_bndry(0.00001,0.00001);  //will add onto SystemSize so include bndry for domain distributions among processors in SimSize
    setSimSize(Total_systemSize*latticeConstant + include_bndry);
    std::cout<<"system size = " << m_systemSize <<std::endl;
    std::cout<<"num_atoms = " << num_atoms() <<std::endl;

    m_num_atoms = num_atoms();

    setSubSystemSize(latticeConstant*subsystemSize);

}




void System::createRandomPositions(int num_particles, double side_length, double temperature, double mass){

    //Places particles randomly into a cube
    for(int i=0; i<num_particles; i++) {      //for all atoms (100 right now)
        Atom *atom = new Atom(UnitConverter::massFromSI(mass)); //uses mass in kg
        //Choose random x,y,z positions
        double x = Random::nextDouble(0, side_length); // random number in the interval [0,10]. nextDouble is defined in random.h
        double y = Random::nextDouble(0, side_length);

        atom->position.set(x,y);    //set atom to position x,y,z

        atom->resetVelocityMaxwellian(temperature);
        m_atoms.push_back(atom);     //add element to vector m_atoms 1 element (atom object)

        vec2 system_size;                                           //define a vector object containing system size
        system_size.set(side_length,side_length);       //set the values of elements of the vector
        setSystemSize(system_size); //system size defines size of system in EACH dimension
        m_num_atoms = num_atoms();
    }
}

void System::increaseTemperature(StatisticsSampler &statisticsSampler, double increment){
    //increases system temperature by factor: meant for use at each MD step to slowly heat system.
    double velocity_rescaling_factor = sqrt((statisticsSampler.temperature()+increment)/statisticsSampler.temperature()); //sqrt(T_new/T_current)=sqrt(factor)
    for(Atom *atom : atoms()) {
        atom->velocity *= velocity_rescaling_factor;  //a*=b means a = a*b
    }
}

void System::calculateForces() {
    for(Atom *atom : m_atoms) {
        atom->resetForce();
    }
    m_potential.calculateForces(*this); // this is a pointer, *this is a reference to this object
}

void System::step(double dt) {
    m_integrator.integrate(*this, dt);  //this calls velocityverlet.cpp
    m_steps++;
    m_time += dt;

}

//---------------------------------------------------------------------
/*!
 Remove atoms from the system.  Needs to sort indices because erase() operation reorders things; also, because of this it is fastest to pop from lowest to highest index.
 Returns the number of atoms deleted.
 \param [in] indices Vector of local indices of atoms to delete from the system
*/
int System::delete_atoms (std::vector <int> indices) {
        std::vector <Atom*>::iterator it = m_atoms.begin();
        std::map <int, int>::iterator map_it;

        // Sort indices from lowest to highest
        sort (indices.begin(), indices.end());

        // Pop in this order
        int shift = 0;
        int upper;
        for (unsigned int i = 0; i < indices.size(); ++i) {
                // Update glob_to_loc for all atoms following the erased atoms
                if (i < indices.size()-1) {
                        upper = indices[i+1]-shift;
                } else {
                        upper = m_atoms.size();
                }
                for (int j = indices[i]+1-shift; j < upper; ++j) {
                        glob_to_loc_id_[m_atoms[j]->atom_index] -= (shift+1);
                }

                map_it = glob_to_loc_id_.find(m_atoms[indices[i]-shift]->atom_index);
                m_atoms.erase(it - shift + indices[i]); //must pass it a accept a const iterator: note: 'it' is an iterator
                glob_to_loc_id_.erase(map_it);
                ++shift;
        }
        m_num_atoms -= shift;
        return shift;
}
/* I think this version is not needed --> is for 1 atom only it seems...
/*
 Attempt to push ONE atom into the system.  This assigns the map automatically to link the atoms global index to the local storage location.
 This reallocates the internal vector that stores the atoms;
 \param [in] natoms Length of the array of atoms to add to the system.
 \param [in] \*new_atoms Pointer to an array of atoms the user has created elsewhere.

std::vector <int> System::add_atoms (const int natoms, Atom &new_atoms) {
        int index = m_atoms.size();
        std::vector <int> update_proc(natoms);  //make vector of size natoms
        for (int i = 0; i < natoms; ++i) {
                m_atoms.push_back(new_atoms[i]);
                glob_to_loc_id_[new_atoms[i].atom_index] = index;
                update_proc[i] = new_atoms[i].atom_index; //filled with indices of the newly added atoms
                index++;
        }
        m_num_atoms += natoms;
        return update_proc;
}
*/

/*!
 Attempt to push an atom(s) into the system.  This assigns the map automatically to link the atoms global index to the local storage location.
 This reallocates the internal vector that stores the atoms.
 \param [in] natoms Length of the array of atoms to add to the system.
 \param [in] \*new_atoms Pointer to an array of atoms the user has created elsewhere.
*/
std::vector <int> System::add_atoms (std::vector <Atom*> new_atoms) {
        int index = m_atoms.size(), natoms = new_atoms.size();  //index is correct --> gives 200 = total number of atoms
        std::vector <int> update_proc(natoms);
        //gets through here

        std::cout <<"index at line 248"<< index <<std::endl;
        std::cout <<"natoms at line 248"<< natoms <<std::endl;  //finds that there are 100 new atoms --> ok.. make sense for 2 processors...
        //issue is that it might try to add atoms that already are in the processor AGAIN... b/c intially didn't split up the domain...
        for (int i = 0; i < natoms; ++i) {
               m_atoms.push_back(new_atoms.at(i)); //note: at() is member fnc of c++ vector class: Returns a reference to the element at position n in the vector. is almost same as [] operator but here checks whether i is within bounds of the vector.
               //gets through here fine
               std::cout << "line 252 system" <<std::endl;
               std::cout <<"new_atoms[0]" << new_atoms[0] <<std::endl;  //new atoms[0] is a pointer --> good
               std::cout <<"new_atoms[0]->atom_index" <<new_atoms[0]->m_mass << std::endl; //it can't find the new_atoms.at(i).atom_index!  //IT SEEMS THE POINTER, BUT CAN'T FIND atom_index!!
               glob_to_loc_id_[new_atoms.at(i)->atom_index] = index;  //BREAKS HERE!  //it has atom_index = 0...
               std::cout <<"line 252 in addatoms" <<std::endl;
               update_proc[i] = new_atoms.at(i)->atom_index; //note: if accessing member property through a pointer, must use -> instead of .
               index++;


        }
        m_num_atoms += natoms;
        return update_proc;
}

/*!
 Attempt to push ghost atom(s) into the system.  This assigns the map automatically to link the atoms global index to the local storage location.
 This reallocates the internal vector that stores the atoms; if a memory error occurs during such reallocation, an error is given and the system exits.
 Does not change m_num_atoms (the number of atoms a processor is responsible for)
 \param [in] natoms Length of the array of atoms to add to the system.
 \param [in] \*new_atoms Pointer to an array of atoms the user has created elsewhere.
 */
void System::add_ghost_atoms (const int natoms, std::vector <Atom*> new_atoms) {
        for (int i = 0; i < natoms; ++i) {
            /* Only add the atom to the system if it is not already contained in the system.
                        Do this by going through the atoms and comparing atom_index of the atoms already in the system and the atoms to be added */
            bool found_atom = false;
            for (unsigned int j = 0; j < m_atoms.size(); ++j) {
                if (m_atoms[j]->atom_index == new_atoms.at(i)->atom_index) {
                    found_atom = true;
                    break;
                }
            }
            if (!found_atom) {
                m_atoms.push_back(new_atoms.at(i));
            }

        }
        return;
}

/*!
 Clears the atoms communicated from neighbouring domains from the list of atoms stored in the system leaving only the atoms the system is responsible for.
 */
void System::clear_ghost_atoms () {
        m_atoms.erase(m_atoms.begin()+m_num_atoms, m_atoms.end());
}
