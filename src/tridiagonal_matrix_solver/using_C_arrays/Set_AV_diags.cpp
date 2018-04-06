#include <iostream>


//this is a in tridiag_solver
double* set_main_diag(double *epsilon, double *main_diag, int num_elements){
    //int num_elements = sizeof(main_diag)/sizeof(main_diag[1]);

    //std::cout<<"numelements " <<num_elements<< std::endl;

    for (int i=1; i<= num_elements;i++){
        main_diag[i] = -2.*epsilon[i+1];     //NOTE:  NEED TO CONSIDER IF THIS SHOULD BE i+1 or i....
    }
    return main_diag;
}

//this is b in tridiag_solver
double* set_upper_diag(double *epsilon, double *upper_diag, int num_elements){
    //int num_elements = sizeof(upper_diag)/sizeof(upper_diag[1]);
    
    for (int i = 1; i<=num_elements-1; i++){
        upper_diag[i] = epsilon[i+1];    //1st element here corresponds to 2nd row of matrix
    }
    return upper_diag;
}

//this is c in tridiag_solver
double* set_lower_diag(double *epsilon, double *lower_diag, int num_elements){
    //int num_elements = sizeof(lower_diag)/sizeof(lower_diag[1]);
    
    for (int i = 1; i<=num_elements-1; i++){
        lower_diag[i] = epsilon[i+1];    //1st element here corresponds to the 1st row of matrix: so need to use i.e. epsilon corresponding to fullV(3) 
    }
    return lower_diag;
    
}
