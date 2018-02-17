#include "lennardjones.h"
#include "system.h"
#include <cmath>
#include <mpi.h>

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

    const double skin_cutoff = 3.5*m_sigma;
    const double skin_cutoff_sqrd = skin_cutoff*skin_cutoff;

    vec2 sys_size = system.subsystemSize(); //returns size of LOCAL processors system box
    int decomp_dim = 0;  // 0 or 1, x or y direction of decomposition
    m_potentialEnergy = 0;  //reset potential energy
    int nprocs, rank;
    MPI_Comm_size (MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

    // Store number of atoms to send to and receive from the processor on the left and on the right
    int num_to_left = 0;
    int num_to_right = 0;
    int num_from_left, num_from_right;

    // Store atoms to send and receive
    std::vector<Atom> to_left;
    std::vector<Atom> to_right;
    std::vector<Atom> from_left;
    std::vector<Atom> from_right;

    MPI_Request req[4], req2[4];
    MPI_Status stat[4], stat2[4];
    
    
   //First calculate interactions btw. all atoms in the system

    for(int current_index=0; current_index<system.num_atoms()-1; current_index++){  //-1 b/c don't need to calculate pairs when get to last atom
        Atom *current_atom = system.atoms(current_index);  //system.atoms(index) returns the pointer to the atom corresponding to index

        for(int other_index=current_index+1;other_index<system.num_atoms();other_index++){   //to avoid double counting
            Atom *other_atom = system.atoms(other_index);


            //distance and vector btw objects dx
              vec2 displacement(0.,0.);

             for(int j=0;j<2;j++){
                 displacement[j] = current_atom->position[j] - other_atom->position[j];
                   //no min. image convention here for parallel version
                //for cases where the folded back particle will be closer than its image to a given particle
                // if (displacement[j] >  system.halfsystemSize(j)) displacement[j] -= system.systemSize(j);   //systemSize(j) returns m_systemSize[j] from system class
                //if (displacement[j] <= -system.halfsystemSize(j)) displacement[j] += system.systemSize(j);
           }


            double radiusSqrd = displacement.lengthSquared();
            //apply interaction range cutoff
            //std::cout<<radius<<std::endl;
            // double radius = sqrt(radiusSqrd);
            //if(radius > 5.0*m_sigma) continue;


            if(radiusSqrd > skin_cutoff_sqrd) continue;

            double radius = sqrt(radiusSqrd);
            double sigma_over_radius = m_sigma/radius;


            double total_force_over_r = 24.*(2.0*pow(radius,-14.)-pow(radius,-8.));
            //ATTRACTIVE FORCE SHOULD POINT TOWARDS OTHER ATOM. REPULSIVE AWAY FROM OTHER ATOM!!!

            //find and set force components
            //double total_force_over_r = total_force/radius; //precalculate to save 2 FLOPS
            for(int j=0;j<2;j++) {
                current_atom->force[j] += total_force_over_r*displacement[j]; //i.e. Fx = (F/r)*x
                other_atom->force[j] -= current_atom->force[j]; //using Newton's 3rd law
            }
            

            
             /*
             if(system.steps() % system.m_sample_freq ==0){
                 //calculate potential energy every m_sample_freq steps
             //m_potentialEnergy += m_four_epsilon*(pow(sigma_over_radius,12.)-pow(sigma_over_radius,6));
                 //below version is if make in unit converter L0 = sigma*1e-10;
                 //since epsilon is set to 1
                 m_potentialEnergy += 4.*(pow(radius,-12.)-pow(radius,-6));
             //m_potentialEnergy += m_four_epsilon*(pow(radius,-12.)-pow(radius,-6));
             }
              */

        }//end of inner loop
        
        //remember if atom need to be sent to neighboring processor
        // special case for nprocs==2; don't want to send same atom twise
        if (nprocs == 2) {                                //SHOULD BE USING SUBSYSTEM SIZE HERE!
            if (system.atoms(current_index)->position[decomp_dim] < rank * sys_size[decomp_dim] / nprocs + skin_cutoff) {
                to_left.push_back(*(system.atoms(current_index)));
                num_to_left++;
            }
            else if  (system.atoms(current_index)->position[decomp_dim] > (rank + 1) * sys_size[decomp_dim] / nprocs - skin_cutoff) {
                to_right.push_back(*(system.atoms(current_index)));
                num_to_right++;
            }
        } else {
            if (system.atoms(current_index)->position[decomp_dim] < rank * sys_size[decomp_dim] / nprocs + skin_cutoff) {
                to_left.push_back(*(system.atoms(current_index)));
                num_to_left++;
            }
            if  (system.atoms(current_index)->position[decomp_dim] > (rank + 1) * sys_size[decomp_dim] / nprocs - skin_cutoff) {
                to_right.push_back(*(system.atoms(current_index)));
                num_to_right++;
            }
        }
            
    }//end of outer loop
    
    //first just test w/o worrying about ghost atoms
    
    
    //send ghost atoms
//     if (nprocs > 1) {
// 		// Send ghost atoms to neighboring processor
// 		MPI_Isend(&num_to_left, 1, MPI_INT, (rank - 1 + nprocs) % nprocs, 1, MPI_COMM_WORLD, req);
// 		MPI_Irecv(&num_from_left, 1, MPI_INT, (rank - 1 + nprocs) % nprocs, 1, MPI_COMM_WORLD, req+1);
// 		MPI_Isend(&num_to_right, 1, MPI_INT, (rank + 1) % nprocs, 1, MPI_COMM_WORLD, req+2);
// 		MPI_Irecv(&num_from_right, 1, MPI_INT, (rank + 1) % nprocs, 1, MPI_COMM_WORLD, req+3);
// 		MPI_Waitall (4, req, stat);
// 
// 		from_left.resize(num_from_left);
// 		from_right.resize(num_from_right);
// 		
// 		MPI_Isend(&to_left[0], num_to_left, MPI_ATOM, (rank - 1 + nprocs) % nprocs, 1, MPI_COMM_WORLD, req2);
// 		MPI_Irecv(&from_left[0], num_from_left, MPI_ATOM, (rank - 1 + nprocs) % nprocs, 1, MPI_COMM_WORLD, req2+1);
// 		MPI_Isend(&to_right[0], num_to_right, MPI_ATOM, (rank + 1) % nprocs, 1, MPI_COMM_WORLD, req2+2);
// 		MPI_Irecv(&from_right[0], num_from_right, MPI_ATOM, (rank + 1) % nprocs, 1, MPI_COMM_WORLD, req2+3);
// 		MPI_Waitall (4, req2, stat2);
        
        
    
    
    
    
    
    
    
    
    


}



