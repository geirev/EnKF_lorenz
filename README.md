This is an implementation of the EnKF and EnKS with the 
Lorenz equations.   It still uses the old traditional implementation
of the analysis scheme.

Author: Geir Evensen
e-mail: geir.evensen@gmail.com

to compile: 
    edit the makefile to select f90 compiler and options.
    run make

To run program:
    cd ./Run
    ../lorenz
           
The file `infile'  defines a number of parameters for the
simulation.  The mode defines which method to run, i.e. pure ensemble
integration with no assimilation to create a climatology,  enkf, enks
and a lagged enks where info is only carried a certain distance backward
in time.
40.0             ! Final time
1.0              ! Model error variance
2.0              ! Initial error variance
2.0              ! Data error variance
1                ! Read observations from file mesA.dat (1)
2                ! mode: ens=0, ekf=1, eks=2, lag=3

