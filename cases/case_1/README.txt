Description:

This is a CFD project in which steady incompressible flow in region with many obstructions is numerically solved using Semi Implicit Methods for Pressure Linked Equations (SIMPLE) solver.

Guidelines to run on Unix based system:

$ make
# It will make the run file

$ ./run
# This will start the numerical computation

If gnuplot is installed the results (while the simulation is going on) may be viewed by using the following commands on the command line:

$ gnuplot
gnuplot> set size ratio 1.0/3;p 'F1vm' w image,'F2vm' w image,'F3vm' w image,'F4vm' w image,'G12vm' w image,'G23vm' w image,'G34vm' w image,'R1vm' w image,'R2vm' w image

# The above plots the x-velocity in the domain
# In order to view the y-velocity or pressure distribution in the domain replace 'u' in 'R1u', etc. with 'v' or 'p' respectively.

In order to check the residuals, use the following command in gnuplot.
gnuplot> p './res_u' w l,'./res_v' w l,'./res_p' w l

In order to remove all case files written during a previous run of the executable 'run' and use it again, do the following:

$ make resetf
# This will delete all case file written during a previous run

Contact:
In case of error or difficulty, email anupambiswas85@gmail.com.
