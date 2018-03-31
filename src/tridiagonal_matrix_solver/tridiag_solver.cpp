// % !==================================================================================
// % ! CRSSolver = Compressed Row Storage (CRS] + Linear Equations Solver
// % !
// % ! This is a solver for a set of linear equations in the form of matrix
// % ! Ax = b  where A = nxn matrix, x and b is vector.
// % ! where matrix A is a "TRIDIAGONAL matrix".
// % ! 
// % !
//    This uses the technique from "Inversion of Jacobi's Tridiagonal Matrix" R. Usmani, 1993
//   The variable notation is almost the same, except, here zeta is theta in Usmani.
//
//% ! NOTE: Inverse matrix will be calculated row by row, then used to calculate the
//% ! answer x element by element. I'll not keep the whole inverse matrix. 


//TriCRSSolver will return a pointer to the answer array x. Takes pointers to the arrays as input
// a = array containing elements of main diagonal. indices: (a1.....an)
// b = array containing elements of upper diagonal. indices (b1....b_n-1)
//c = array containing elements of lower diagonal. indices (c1...c_n-1)

// % !==================================================================================


#include <iostream>
#include <numeric>  //allows to use inner product fnc.


double* TriCRSSolver(double *a, double *b, double *c, double *rhs, int num_elements){  //later should figure out how it could find the num_elements itself, using length of main diag: a

    
    //In order to get a dynamic array need to use a pointer and new, can't define an array with [variable] for it's lenght...., will get "expression did not evaluate to a constant" error
    double *zeta = new double[num_elements+1];
    double *phi= new double[num_elements+1];
    double *A_inverse = new double[num_elements];
    double temp;
    
    double* x = new double[num_elements];  //x is defined as a pointer to an array of num_elements doubles

    //[a, b, c] = Get_abc(A_val, num_elements); //INSTEAD OF DOING THIS UNPACKING, JUST DIRECTLY PASS a,b,c to this fnc.!!!--> this makes it less general, b/c
    //in code calling this fnc will need to specify a,b,c (the 3 diagonals of tridiag) matrix, instead of using 1 array A_val, but really isn't much differeence.

    //std::cout <<"c[num)el] " <<c[num_elements-1] <<std::endl;

    zeta[0] = 1.0;
    zeta[1] = a[1];
    for(int i = 2; i<=num_elements; i++){
        zeta[i] = a[i]*zeta[i-1] - b[i-1]*c[i-1]*zeta[i-2];
    }
   
    phi[num_elements+1] = 1.0;
    phi[num_elements] = a[num_elements];
    
    for(int i = num_elements-1; i>0; i--){  //decreasing iterations.., for recursion relation
        phi[i] = a[i]*phi[i+1] - b[i]*c[i]*phi[i+2];
    }
    
    //std::cout <<"= " <<phi[num_elements-1] <<std::endl;


   //calculate inverse matrix 
    for(int i = 1; i<=num_elements; i++){
        for (int j = 1; j<= i-1; j++){
            temp = 1.0;
            for(int k = j; k<=i-1; k++){  //k<= i-1 b/c there's 1 less element in off-diagonals than main diagonal (which is indexed by i)
                temp = temp*c[k];   //multiplies all the c's together--> also is what sometimes causes Usmani method to blow up
            }
            A_inverse[j] = pow(-1., i+j)*temp*zeta[j-1]*phi[i+1]/zeta[num_elements];
            //std::cout <<"a inverse " << A_inverse[j] <<std::endl;
            //here seems fine
        }
        for(int j = i; j<=num_elements;j++){
            temp = 1.0;
            for(int k = i; k<=j-1;k++){
                temp = temp * b[k];
            }
            A_inverse[j] = pow(-1., i+j)*temp*zeta[i-1]*phi[j+1]/zeta[num_elements];
            //std::cout <<"a inverse " << A_inverse[j] <<std::endl;
        }
        
        double init = 0;

        //hardcode the dot product
        x[i] = 0;
        for(int m = 1; m<= num_elements;m++){
            x[i] += A_inverse[m]*rhs[m];
        }

        //x[i] = std::inner_product(A_inverse,A_inverse+num_elements,rhs,init); //SEEMS THIS INNER PRODUCT DOESN'T WORK!
        //inner_product (InputIterator1 first1, InputIterator1 last1, InputIterator2 first2, T init); //init is initial value of scalar product, which we want to be 0
    }

    //%test[i] = A_inverse*rhs
return x;
}

