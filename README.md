# materials-analysis

Follow a breafly description of each program usage:

********* eig.f ************
This program read the EIGENVAL file from vasp extracting  two
column file(s) with the band dispersion E(k), with E in eV and 
k in 2pi/a.
You should compile the program using gfortran (or other fortran
compiler). Edit the file inp.dat as described within it. After
executing the program one or two files will be generated depending
on the spin degree of freedom.

********* grid-spin.f ************
This program will extract the spin-texture (S_x, S_y, S_z)
as a function of kx and ky from the PROCAR file of VASP when runned in 
a uniform grid with SOC. 
You should compile the program using gfortran (or other fortran
compiler). Edit the file inp-grid.dat as described within it. After
executing the program one file containing the spin-texture will be
generated.

********* SQS.f ************
This program read a POSCAR file containing the sites needed 
for the alloying, creating a new POSCAR with the alloyed structure.
You should compile the program using gfortran (or other fortran
compiler). During execution the program will ask for a superstructure cell
size, two interger defining the alloying fraction, and three real numbers
for first, second and thyrd-neghbors correlation threshold.

********* gen-kp-grid.sh *********
This simple program creats a KPOINTS file for vasp containing an uniform
NxN grid of points, with N and the origin of the grid define by the user
editing the gen-kp-grid.sh file
