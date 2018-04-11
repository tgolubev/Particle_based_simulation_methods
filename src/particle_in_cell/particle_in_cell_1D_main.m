%1D Particle-in-Cell code with periodic boundary conditions.
%Set up are xy planes of charge seperated in the z direction.
clear all;

%Timofey Golubev
%% Parameters

num_cells = 30;  
num_particles_x = 1;  %number of particles in the x direction OF THE SHEET
num_particles_y = 1;
particles_per_cell = num_particles_x*num_particles_y;
num_particles = particles_per_cell*num_cells; %these are COMPUTATIONAL particles 

conversion = 10^-10;  %to convert to Ansgtroms for movie file--> OVITO uses Angstroms

distance_unit = 10^-9;  %distance units in meters
L = 10000000*distance_unit;  %length of the system  %NOTE: seems that this # must be very high for it to oscillate (BUT ONLY 1 PLANE OSCILLATES!)..--> 
%NOTE: ALL THIS IS DOING IS SETTING THE SIZE OF THE COMPUTATIONAL
%PARTICLES, NOT THE ACTUAL ELECTRON SPACINGS!!!-->  B/C THE ELECTRON
%DENSITY IS SET BY THE "electron_density" below
% just is saying that the density of computational particles should be
% smaller than electrons (i.e. many electrons / 1 comp particle--> 
%makes perfect sense, comp particle should be more macroscopic quantity

dz = L/num_cells;
epsilon_0 = 8.854*10^-12;
e = 1.6*10^-19; 
initial_spacing = 0.5*dz;  %INITIAL SPACING FOR THE y and z directions initially space all electrons equidistant
%perturb = 0.2*initial_spacing;  %a small perturbation of electrons positions from the ions positions in order to generate oscillations
m_e = 9.10938356*10^-31;  %electron mass (kg)

system_volume = L*(particles_per_cell*initial_spacing^2)

dt = 10^-16;  %10^-16 time step allows it to run with a lower 
N_steps = 10000;
particle_V = initial_spacing^2*dz;  %volume of each computational particle

%NOTE: THE COMOPUTATIONAL PARTICLE DENSITY DOES NOT NEED TO CORRESPOND THE
%THE REAL ELECTRON DENSITY--> CAN HAVE MORE electrons per 1 computational
%particle.
electron_density = 10^16    %changed from 10^16 %density of typical plasma
rho_background = e*electron_density;       %this specifies the density of plasma we have, since neutralizing background + integral over rho_particles = 0 net charge.
% so I'm assumign that the plasma has ON AVERAGE: 10^16 electron/m^3, which are
% neutralized by the rho_background
%used source: https://arxiv.org/ftp/arxiv/papers/1404/1404.0509.pdf
Q = -rho_background*particle_V;   %particle charge: Q = averaged charge density*volume which particle takes up: 
sigma =  Q/(initial_spacing)^2;         %surface charge density (along the xy planes) = total particle charge/size of each particle's xy surface. This is used for force calculation
m_particle = electron_density*m_e*particle_V;

initial_vel = 0;%10^6;  %initial velocity for electrons



%% Initializations

%allocate poisson matrix
A = zeros(num_cells,num_cells);
CV = -dz^2/epsilon_0;  %coefficient for Poisson eqn RHS

% spacing diagram: | .  .  .  .  . | .  .  .  .  . |  %these are showing
% repeating versions of the entire system!, not individual cells!

%pre-allocate num_particles in each cell to 0 (for the count to work)
for i = 1:num_cells
    cell(i).num_particles = 0;
end

%Initialize computational particles evenly distributed in the cells
i = 1;
prev_random_num = 0;
for current_cell = 1:num_cells
    random_num = rand(1);
    while(abs(prev_random_num -random_num) < 0.5*dz)
        random_num = rand(1);          %make sure initialization is not too close
    end
    pos_z = (current_cell-1)*dz + rand(1)*dz; %use random # to get initial positions (rand(1) generates random number btw 0 and 1).
    %pos_z = (current_cell-1)*dz + dz/4;     %make x position be same for all particles within each cell -> plane of charge%starting position of the slabs is 1 particle at the dz/4 position, inside each cell
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
            %particle(i).velocity = 10^-10*rand(1);         %TRY TO USE A RANDOM VELOCITY
            i = i+1;
        end
    end
%     prev_random_num = random_num;
end

%create file
movie = fopen(fullfile('C:\Users\Tim\Documents\CMSE890\particle_in_cell_methods\results\movie.xyz'),'w'); % w: open file for writing

%% Main loop
step = 0;
while step<= N_steps
    added_cell = 0;
    
    %MUST rezero the num_particles in each time step!!
    for i = 1:num_cells
        cell(i).num_particles = 0;
    end
    
     
    for i = 1:num_particles
        %assign each particle to its cell
        cell_index = floor(particle(i).pos_z/dz)+1;  %+1 b/c we can't index from 0
        
        %to prevent the negative indices problem BUT NOT SURE IF THIS WILL
        %MESS UP THE PHYSICS
       % if(cell_index <1) 
        %    cell_index = 1;
        %end
        
        if cell_index > num_cells   %NOTE: once we update this, num_cells gets +1, so next time, if need another cell, will still enter this loop
            cell(cell_index).num_particles = 0;
            num_cells = num_cells + 1;
        end

        particle(i).cell = cell_index;  %keep track of which cell the particle is in
        %this finds indices correctly
        cell(cell_index).num_particles = cell(cell_index).num_particles +1;  
        %cell(cell_index).num_particles;
    end
    step
    
    %calculate charge density on each cell

    %sum over Wz's (the shape/weighting functions)
    %use Wz = 1 when in left half of cell, and 0 when in right half of cell
    % 1 when z<(delta z)/2
    % 0 when z>(delta z)/2
    for cell_index = 1:num_cells
        Wz_sum = 0;
        for i = 1:num_particles
            if abs(particle(i).pos_z - (cell_index-1)*dz) < dz/2  %subtract off the start of each cell's position, to see if particle is located in left 1/2 or right 1/2 of the cell
                Wz = 1;
            else
                Wz = 0;
            end
            Wz_sum = Wz_sum + Wz;
        end
        
        %compare the Wz sum with the # of particles in each cell
        %Wz_sum
        %cell(cell_index).num_particles
        
        cell(cell_index).rho = rho_background*(1- (num_cells/num_particles)*Wz_sum); %this gives the total charge on the grid (in each cell)

        %setup rhs of Poisson matrix equation
        rhs(cell_index) = CV*cell(cell_index).rho;  %NOTE:  check the constants CV    
    end
    
    %Re-setup Poisson matrix b/c the # of cells could have changed which
    %changes matrix size
    %setup Poisson matrix off-diagonals
    for cell_index = 2:num_cells -1
        A(cell_index,cell_index+1) = 1;
        A(cell_index,cell_index-1) = 1;
    end
    %setup Poisson matrix main diagonal
    for cell_index = 1:num_cells

        A(cell_index,cell_index) = -2;
    end
    %explicitely  define the elements which are missed in the loop
    A(1,2) = 1;
    A(num_cells, num_cells-1)= 1;
    %define the 2 non-tridiag elements which impose the PBCs
    A(1,num_cells) = 1;
    A(num_cells,1) = 1;

    V = A\(rhs');

    %% Calculate Electric field in each cell
    for i = 2:num_cells-1
        E(i) = (V(i+1)-V(i-1))/(2*dz);
    end
    %for boundary cells use PBC
    %E(1) = 0;
    %E(num_cells) = 0;
    E(1) = (V(2) - V(num_cells))/(2*dz);
    E(num_cells) = (V(1) - V(num_cells-1))/(2*dz);


    %% Calculate Forces felt by each electron
    %F = sigma* sum(E(grid pts)* Wz        %need to use the shape fnc's again
    for i = 1:num_particles
        EWz_sum = 0;
        for cell_index = 1:num_cells
            if abs(particle(i).pos_z - (cell_index-1)*dz) < dz/2  %subtract off the start of each cell's position, to see if particle is located in left 1/2 or right 1/2 of the cell
                EWz = E(cell_index)*1;
            else
                EWz = 0;  %E*Wz = 0 b/c Wz= 0
            end
            EWz_sum = EWz_sum + EWz;
        end
        particle(i).force = sigma*EWz_sum;  %based on E at the location fo electron
        
        
        %for testing
        %OUTPUT ALL FORCES in 1st step felt by each plane..--> do for case
        %of 1 particle per plane
        for i = 1:num_cells
            particle(i).force
        end
        
    end
    
    
   %save to file BEFORE proceed to next iter--> so can save the initial
   %state
   if(step == 0)
    fprintf(movie,'%.0f \r\n', num_particles);
    fprintf(movie,'%.0f \r\n', step);      
    for i = 1:num_particles
        fprintf(movie,'H %.4e %.4e %.4e %.4e  \r\n ',particle(i).pos_x/conversion, particle(i).pos_y/conversion, particle(i).pos_z/conversion,  particle(i).velocity);
        %electron(i).velocity
        %electron(i).position
    end
   end
    

    %% Leapfrog integrator
    %for 1st timestep, need to use 1/2 Euler step to find velocity
    if(step == 0)
        for i = 1:num_particles
         particle(i).velocity = initial_vel - (dt/2)*particle(i).force./m_particle;  %this is v(-delta t/2)
        end
    end
    
    for i = 1:num_particles
        particle(i).velocity = particle(i).velocity + dt*particle(i).force./m_particle; 
        
        particle(i).pos_z = particle(i).pos_z + dt*particle(i).velocity;
        
        %Apply PBC'S
        if particle(i).pos_z > L 
            particle(i).pos_z = particle(i).pos_z -L;
        elseif particle(i).pos_z < 0
           particle(i).pos_z = particle(i).pos_z + L;
        end
            
    end
    
    %OUTPUT FOR TESTING
%    particle(100).velocity
   % particle(100).pos_z
   % particle(100).force
    
    step = step +1;
    
    fprintf(movie,'%.0f \r\n', num_particles);
    fprintf(movie,'%.0f \r\n', step);      
    for i = 1:num_particles
        fprintf(movie,'H %.4e %.4e %.4e %.4e  \r\n ',particle(i).pos_x/conversion, particle(i).pos_y/conversion, particle(i).pos_z/conversion,  particle(i).velocity);
        %electron(i).velocity
        %electron(i).position
    end
    
    %for plotting
    particle_1_z(step) = particle(1).pos_z;
%    particle_100_z(step) = particle(100).pos_z;
    %particle_250_z(step) = particle(250).pos_z;
    
    steps(step) = step;
    
end

fclose(movie);

figure;
plot(steps, particle_1_z);

figure;
plot(steps, particle_100_z);

figure;
plot(steps, particle_250_z);




