%1D Particle-in-Cell code with periodic boundary conditions.
%Set up are xy planes of charge seperated in the z direction.
clear all;

%Timofey Golubev
%% Parameters

num_cells = 10;  
num_particles_x = 10;  %number of particles in the x direction OF THE SHEET
num_particles_y = 10;
particles_per_cell = num_particles_x*num_particles_y;
num_particles = particles_per_cell*num_cells; %these are COMPUTATIONAL particles 

distance_unit = 10^-9;  %distance units in meters
L = 10*distance_unit;  %length of the system
dz = L/num_cells;
epsilon_0 = 8.854*10^-12;
e = 1.6*10^19; 
initial_spacing = dz;  %INITIAL SPACING FOR THE y and z directions initially space all electrons equidistant
%perturb = 0.2*initial_spacing;  %a small perturbation of electrons positions from the ions positions in order to generate oscillations
m_e = 9.10938356*10^31;  %electron mass (kg)

dt = 10^-15;
N_steps = 1000;
cell_V = L^2*dz;  %volume of each cell

electron_density = 10^16;
rho_background = e*electron_density;       %this specifies the density of plasma we have, since neutralizing background + integral over rho_particles = 0 net charge.
% so I'm assumign that the plasma has ON AVERAGE: 10^16 electron/m^3, which are
% neutralized by the rho_background
%used source: https://arxiv.org/ftp/arxiv/papers/1404/1404.0509.pdf
particle_charge = -rho_background*cell_V;   %particle charge = averaged charge density*volume which particle takes up: 
m_particle = electron_density*m_e*cell_V;

initial_vel = 0;%10^6;  %initial velocity for electrons



%% Initializations

%allocate poisson matrix
A = zeros(num_cells,num_cells);

% spacing diagram: | .  .  .  .  . | .  .  .  .  . |  %these are showing
% repeating versions of the entire system!, not individual cells!

%pre-allocate num_particles in each cell to 0 (for the count to work)
for i = 1:num_cells
    cell(i).num_particles = 0;
end

%Initialize computational particles evenly distributed in the cells
i = 1;
for current_cell = 1:num_cells
    pos_z = (current_cell-1)*dz + dz/4;      %make x position be same for all particles within each cell -> plane of charge%starting position of the slabs is 1 particle at the dz/4 position, inside each cell
    for m = 1:num_particles_x     
        for n = 1:num_particles_y
            particle(i).name = 'name';
            particle(i).pos_z = pos_z;  
            %particle(i).pos_z
            particle(i).pos_x = (m-1)*initial_spacing + initial_spacing/2;  %+dz/2 to offshift away from the boundaries
            particle(i).pos_y = (n-1)*initial_spacing + initial_spacing/2;
            %particle(i).position
            %position is being calculated correctly
            particle(i).force = 0;  %initialize with 0 force felt
            particle(i).velocity = initial_vel;
            i = i+1;
        end
    end
end

for i = 1:num_particles
    %assign each particle to its cell
    cell_index = floor(particle(i).pos_z/dz)+1;  %+1 b/c we can't index from 0
    %cell_index
    particle(i).cell = cell_index;  %keep track of which cell the particle is in
    %this finds indices correctly
    cell(cell_index).num_particles = cell(cell_index).num_particles +1;  
    %cell(cell_index).num_particles;
end

%calculate charge density on each cell
CV = -dz^2/epsilon_0;
%sum over Wz's (the shape/weighting functions)
%use Wz = 1 when in left half of cell, and 0 when in right half of cell
% 1 when z<(delta z)/2
% 0 when z>(delta z)/2
for cell_index = 1:num_cells
    Wz_sum = 0;
    for i = 1:num_particles
        if (particle(i).pos_z - (cell_index-1)*dz) < dz/2  %subtract off the start of each cell's position, to see if particle is located in left 1/2 or right 1/2 of the cell
            Wz = 1;
        else
            Wz = 0;
        end
        Wz_sum = Wz_sum + Wz;
    end
    cell(cell_index).rho = rho_background*(1- (num_cells/num_particles)*Wz_sum);
    
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
        E(i) = (V(i+1)-V(i-1))/(2*dz);
    end
    %for boundary cells use PBC
    E(1) = (V(2) - V(num_cells))/(2*dz);
    E(num_cells) = (V(1) - V(num_cells-1))/(2*dz);


    %% Calculate Forces felt by each electron
    %F = -eE
    for i = 1:num_particles
        particle(i).force = particle_charge*E(particle(i).cell);  %based on E at the location fo electron
        %particle(i).force
    end
    %ISSUE IS THAT THE ELECTRONS ARE FEELING NO FORCES IF I HAVE IONS, B/C THERE IS NO E
    %field!--> b/c I have no net charge in each cell!!
    %if i DON'T INITIALIZE FORCES, THEN ELECTRONS FEEL FORCES!
    %but if I do electrons only, everything just blows appart!, all
    %repulsive
    
    
   %save to file BEFORE proceed to next iter--> so can save the initial
   %state
    fprintf(movie,'%.0f \r\n', num_particles);
    fprintf(movie,'%.0f \r\n', step);      
    for i = 1:num_particles
        fprintf(movie,'H %.4e %.4e %.4e %.4e  \r\n ',particle(i).pos_x, particle(i).pos_y, particle(i).pos_z,  particle(i).velocity);
        %electron(i).velocity
        %electron(i).position
    end
    

    %% Leapfrog integrator
    %for 1st timestep, need to use 1/2 Euler step to find velocity
    if(step == 0)
        for i = 1:num_particles
         particle(i).velocity = initial_vel + (dt/2)*particle(i).force./m_particle;  %this is v(+delta t/2)
        end
    end
    
    for i = 1:num_particles
        particle(i).velocity = particle(i).velocity + dt*particle(i).force./m_particle; 
        particle(i).pos_z = particle(i).pos_z + dt*particle(i).velocity;
        
        %Apply PBC'S
        if particle(i).pos_z > L 
            particle(i).pos_z = particle(i).pos_z -L;
        elseif particle(i).pos_z < L
           particle(i).pos_z = particle(i).pos_z + L;
        end
            
    end
    
 
    
    step = step +1;
end

fclose(movie);

