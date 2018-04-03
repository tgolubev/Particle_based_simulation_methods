This is a 1D particle-in-cell simulation with periodic boundary conditions. This attempts to calculate
the plasma oscillation frequency. The algorithm flow is:

1. Put particles into cells
2. Calculate the charge in each cell and calculate a charge density.
1D finite-difference discretization of the grid
4. Use finite-differences to find electric field
6. Use leap-frog integrator to move the particles under influence of the calculated forces.
7. Update which particles are in which cells
8. Repeat steps 2-7 for each time step
