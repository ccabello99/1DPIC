# PHY571
Final Project for PHY571


1D Particle in Cell code written in Julia to simulate the interaction of an intense laser pulse with a plasma. The fields are evolved with 1D Finite-Domain Time-Difference methods. The electric potential is calculated from Poisson equation 
and the charge/current density and fields are interpolated to the grid using simple linear interpolation. The macroparticles are weighted and evolved with the relativistic Boris method.
