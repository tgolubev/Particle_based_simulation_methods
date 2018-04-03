%1D Particle-in-Cell code with periodic boundary conditions.
clear vars;

%Timofey Golubev
%% Parameters
num_electrons = 100;  
num_ions = num_electrons;
num_cells = 10;  
L = 10;  %length of the system
dx = L/num_cells;
epsilon_0 = 8.854*10^-12;
e = 1.6*10^19; 
initial_spacing = L/(num_electrons);  %initially space all electrons equidistant
perturb = 0.05*initial_spacing;  %a small perturbation of electrons positions from the ions positions in order to generate oscillations
m_e = 9.10938356*10^31;  %electron mass (kg)
dt = 10^-10;
N_steps = 1000;

initial_vel = 0;  %initial velocity for electrons
cross_section_A =  L^2;  %make it like a cubic system, but only  treat 1D, this is for rho in Poisson
cell_V = cross_section_A*dx;  %volume of each cell


%% Initializations

%allocate poisson matrix
A = zeros(num_cells,num_cells);


% spacing diagram: | .  .  .  .  . | .  .  .  .  . |  %these are showing
% repeating versions of the entire system!, not individual cells!


%initialize charge in each cell to zero
for i = 1:num_cells
    cell(i).charge = 0;
end

%Initialize ions in each cell
% for i = 1:num_ions
%     ion(i).position= (i-1)*initial_spacing + initial_spacing/2;  %dx/2 to shift off of the edges
%     %ion(i).position
%     %position is being calculated correctly
%     cell_index = floor(ion(i).position/dx) +1;  %+1 b/c we can't index from 0
%     %cell_index
%     %this finds indices correctly
%     cell(cell_index).charge = cell(cell_index).charge +1;
%     %cell(cell_index).charge;
% end

%Initialize electrons evenly distributed in the cells
for i = 1:num_electrons
    electron(i).name = 'name';
    electron(i).position= (i-1)*initial_spacing + initial_spacing/2 + perturb;  %dx/2 to shift off of the edges, perturb off of ions positions
    %electron(i).position
    %position is being calculated correctly
    electron(i).force = 0;  %initialize with 0 force felt
    electron(i).velocity = initial_vel;
    
    %assign each electron to its cell
    cell_index = floor(electron(i).position/dx) +1;  %+1 b/c we can't index from 0
    electron(i).cell = cell_index;  %keep track of which cell the electron is in
    %cell_index
    %this finds indices correctly
    cell(cell_index).charge = cell(cell_index).charge -1;  %NOTE: -1 b/c electrons!
    cell(cell_index).charge;
end

%testing
%intially charge in each cell will be 0
%for i = 1:num_cells
%    cell(cell_index).charge
%end

%calculate charge density on each cell
CV = dx^2/epsilon_0;
for cell_index = 1:num_cells
    cell(cell_index).rho = e*cell(cell_index).charge./(cell_V);
    %cell(cell_index).rho;
    %NOTE: initial charge densities are not equal for some reason!
    
    %setup rhs of Poisson matrix equation
    rhs(cell_index) = CV*cell(cell_index).rho;  %NOTE:  check the constants CV!
    
    %setup Poisson matrix main diagonal
    A(cell_index,cell_index) = -2;
    
end

%setup Poisson matrix off-diagonals
for cell_index = 2:num_cells -1
    A(cell_index,cell_index+1) = 1;
    A(cell_index,cell_index-1) = 1;
end
%explicitely  define the elements which are missed in the loop
A(1,2) = 1;
A(num_cells, num_cells-1)= 1;
%define the 2 non-tridiag elements which impose the PBCs
A(1,num_cells) = 1;
A(num_cells,1) = 1;

V = A\(rhs');

%create file
movie = fopen(fullfile('C:\Users\Tim\Documents\CMSE890\particle_in_cell_methods\results\movie.xyz'),'w'); % w: open file for writing

%% Main loop
step = 0;
while step<= N_steps

    %% Calculate Electric field in each cell
    for i = 2:num_cells-1
        E(i) = (V(i+1)-V(i-1))/(2*dx);
    end
    %for boundary cells use PBCs
    E(1) = (V(2) - V(num_cells))/(2*dx);
    E(num_cells) = (V(1) - V(num_cells-1))/(2*dx);


    %% Calculate Forces felt by each electron
    %F = -eE
    for i = 1:num_electrons
        electron(i).force = -e*E(electron(i).cell);  %based on E at the location fo electron
        electron(i).force
    end
    %ISSUE IS THAT THE ELECTRONS ARE FEELING NO FORCES IF I HAVE IONS, B/C THERE IS NO E
    %field!--> b/c I have no net charge in each cell!!
    %if i DON'T INITIALIZE FORCES, THEN ELECTRONS FEEL FORCES!
    %but if I do electrons only, everything just blows appart!, all
    %repulsive

    %% Leapfrog integrator
    %for 1st timestep, need to use 1/2 Euler step to find velocity
    if(step == 0)
        for i = 1:num_electrons
         electron(i).velocity = initial_vel + (dt/2)*electron(i).force./m_e;  %this is v(+delta t/2)
        end
    end
    
    for i = 1:num_electrons
        electron(i).velocity = electron(i).velocity + dt*electron(i).force./m_e; 
        electron(i).position = electron(i).position + dt*electron(i).velocity;
    end
    
    %save to file
    fprintf(movie,'%.0f \r\n', num_electrons);
    fprintf(movie,'%.0f \r\n', step);      
    for i = 1:num_electrons
        fprintf(movie,'H %.4e %.4e  \r\n ',electron(i).position, electron(i).velocity);
        %electron(i).velocity
        %electron(i).position
    end
    
    step = step +1;
end

fclose(movie);

