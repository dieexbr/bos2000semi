# Introduction
Reproduction of the results of Bos and Ruban (2000), regarding steady-state solutions of triple-deck equations for supersonic flow over an adiabatic wall.

The numerical scheme is almost identical to that of Cassel et al. (1995), except for the time derivative and the du / dx term in the convective term, which is now taken implicitly.

Matlab files:
- main.m sets up the scale angle, stretching factors and will record a movie at each iteration.
- start.m will initialise the s struct.
- update.m will progress the solution for one iteration.
- find_* will obtain variable * from previous solution.
- export.m will export the data to a folder called "grid_convergence".

# Contact

Please drop a message at dieexbr17@gmail.com should you have any question about this repository.
