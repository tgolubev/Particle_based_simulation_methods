%Script for testing my C++ tridiagonal solver for 1D Poisson equation by
%solving a tridiagonal matrix using Matlab's solver

num_elements = 100;

V_leftBC = 0.;
V_rightBC= 0.2;

epsilon(1:num_elements+2) = 3.8;

%Setup the tridiagonal matrix
AV = zeros(num_elements,3);
for i=1:num_elements     
    AV(i,2) = -2.*epsilon(i+1);    
end

for i = 1:num_elements-1
    AV(i,1) = epsilon(i+1);    %1st element here corresponds to 2nd row of matrix
end
AV(num_elements, 1) = 0;     %last element is unused

for i = 2:num_elements
    AV(i,3) = epsilon(i+1);      %1st element here corresponds to the 1st row of matrix: so need to use i.e. epsilon corresponding to fullV(3) 
end
AV(1,3) = 0;                 %1st element of Ap(:,3) is unused


AV_val = spdiags(AV,-1:1,num_elements,num_elements);  %A = spdiags(B,d,m,n) creates an m-by-n sparse matrix by taking the columns of B and placing them along the diagonals specified by d.

%Setup the right-hand side of Poisson's equation
for i = 1: num_elements
    rhs(i,1) = 1;   
end

%for bndrys 
 rhs(1) = rhs(1) - epsilon(1)*V_leftBC;
 rhs(num_elements) = rhs(num_elements) - epsilon(num_elements)*V_rightBC;
        
 V = AV_val\rhs