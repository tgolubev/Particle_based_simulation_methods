#include "lennardjones.h"
#include "system.h"
#include "math.h"

double LennardJones::potentialEnergy() const
{
    return m_potentialEnergy;
}

double LennardJones::sigma() const
{
    return m_sigma;
}

void LennardJones::setSigma(double sigma)
{
    m_sigma = sigma;
    m_sigma_sqrd = sigma*sigma;
}

double LennardJones::epsilon() const
{
    return m_epsilon;
}


void LennardJones::setEpsilon(double epsilon)
{
    m_epsilon = epsilon;
    m_four_epsilon = 4.0*m_epsilon;
    m_twntyfour_epsilon = 24.0*m_epsilon;  //must reset this too
}

void LennardJones::calculateForces(System &system)  //object system is passed by reference (allows changing)
{
 m_potentialEnergy = 0;  //reset potential energy

 double skin_cutoff = 3.*m_sigma;

 double too_close = 0.7*m_sigma;
 double too_close_sqrd = too_close*too_close;

 const double skin_cutoff_sqrd = skin_cutoff*skin_cutoff;

 for(int current_index=0; current_index<system.num_atoms()-1; current_index++){  //-1 b/c don't need to calculate pairs when get to last atom
    Atom *current_atom = system.atoms(current_index);  //system.atoms(index) returns the pointer to the atom corresponding to index


    //temporarily test calculating all forces explicitely for pairs  THIS WORKS!, but newton's 3rd law doesn't!
    //for(int other_index=1;other_index<system.num_atoms();other_index++){
      // if(current_index == other_index) continue;

    for(int other_index=current_index+1;other_index<system.num_atoms();other_index++){   //to avoid double counting
        Atom *other_atom = system.atoms(other_index);


        //distance and vector btw objects dx should obey min. image convention: only closest distance to particle or its image is considered
        vec2 displacement(0.,0.);
        for(int j=0;j<2;j++){
             displacement[j] = current_atom->position[j] - other_atom->position[j];
             //for cases where the folded back particle will be closer than its image to a given particle
             if (displacement[j] >  system.halfsystemSize(j)) displacement[j] -= system.systemSize(j);   //systemSize(j) returns m_systemSize[j] from system class
             if (displacement[j] <= -system.halfsystemSize(j)) displacement[j] += system.systemSize(j);
         }
        double radiusSqrd = displacement.lengthSquared();


        //apply interaction range cutoff
      if(radiusSqrd > skin_cutoff_sqrd) {
          //  std::cout << "radius when used cutoff" << sqrt(radiusSqrd) <<std::endl;
         continue;
       }



        // if(radiusSqrd < too_close_sqrd ) {
               // std::cout << "reached too close radius of" << sqrt(radiusSqrd) <<std::endl;
          //      radiusSqrd = too_close_sqrd;   //if gets too close, just change radius to be larger so force doesn't blow up (don't actually change the forces)
         //}



        double radius = sqrt(radiusSqrd);
        double sigma_over_radius = m_sigma/radius;


       //double total_force_over_r = (m_twntyfour_epsilon/pow(m_sigma,2.))*(2.0*pow(sigma_over_radius,14.)-pow(sigma_over_radius,8.));

        //below version is if make in unit converter L0 = sigma*1e-10;
        //note: epsilon was set to 1, so coeff. in front is just 24
        double total_force_over_r = 24.*(2.0*pow(radius,-14.)-pow(radius,-8.));
        //double total_force_over_r = m_twntyfour_epsilon*(2.0*pow(radius,-14.)-pow(radius,-8.));
        //ATTRACTIVE FORCE SHOULD POINT TOWARDS OTHER ATOM. REPULSIVE AWAY FROM OTHER ATOM!!!

        //std::cout << "individual force/r calc" <<total_force_over_r <<std::endl;


        //find and set force components
        //double total_force_over_r = total_force/radius; //precalculate to save 2 FLOPS
        for(int j=0;j<2;j++) {
            current_atom->force[j] += total_force_over_r*displacement[j]; //i.e. Fx = (F/r)*x
            other_atom->force[j] -= total_force_over_r*displacement[j]; //using Newton's 3rd law  NOTE CANNOT TAKE THE NEGATIVE OF THE CURRENT ATOM FORCE B/C THAT IS AN ADDITIVE FORCE!

        }

       // std::cout <<"current atom force within loop" <<current_atom->force[0] <<std::endl;

        if(system.steps() % system.m_sample_freq ==0){
            //calculate potential energy every m_sample_freq steps
           //m_potentialEnergy += m_four_epsilon*(pow(sigma_over_radius,12.)-pow(sigma_over_radius,6));
            //below version is if make in unit converter L0 = sigma*1e-10;
            //since epsilon is set to 1
            m_potentialEnergy += 4.*(pow(radius,-12.)-pow(radius,-6));
           //m_potentialEnergy += m_four_epsilon*(pow(radius,-12.)-pow(radius,-6));
        }
    }
 }


 //output forces

  for(Atom *atom : system.atoms()) {
     // std::cout<< "atom force after Lj LOOP" << atom->force[0] <<" " <<atom->force[1] <<"timestep" <<system.steps() << std::endl;
  }


}



