//This is a code for testing the tridiagonal solver in the context of solving the 1D Poisson equation


//March 30, 2018
//Timofey Golubev


#include <iostream>
#include "Set_AV_diags.h"
#include "tridiag_solver.h"


int main()
{

    int num_elements = 100;  //number of nodes in the system
    double *epsilon = new double[num_elements+2];
    double *a = new double[num_elements]; //main diag
    double *b = new double[num_elements-1]; //upper diag
    double *c = new double[num_elements-1]; //lower diag
    double *V = new double[num_elements];  //array for solution, electric potential
    double *rhs = new double[num_elements];

    double V_leftBC = 0.;
    double V_rightBC= 0.9;

    for(int i=0;i<=num_elements+1;i++){
        epsilon[i] = 3.8;
    }

    a = set_main_diag(epsilon, a,  num_elements);
    b = set_upper_diag(epsilon, b, num_elements);
    c = set_lower_diag(epsilon, c, num_elements);

    //setup rhs of Poisson eqn.
    for (int i = 1;i<= num_elements;i++){
         rhs[i] = 0;  //for now
    }

    //test i.e. having a dipole at accross 24-25th node
    rhs[24] = 0.681; //corresponds to having 5 holes/(75nm)^2 plane at the interface, with z-mesh size of 1nm
    rhs[25] = -0.681;

    //for bndrys
    rhs[1] = rhs[1] - epsilon[1]*V_leftBC;
    rhs[num_elements] = rhs[num_elements] - epsilon[num_elements]*V_rightBC;

    V = TriCRSSolver(a, b, c, rhs,  num_elements);

    //output the result to terminal
    for(int i = 1;i<= num_elements; i++){
        std::cout << V[i] << " " <<std::endl;
    }

    return 0;
}
