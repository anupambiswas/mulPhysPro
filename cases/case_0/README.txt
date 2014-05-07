Description:

This is a CFD project in which steady incompressible flow over a step is numerically solved using Semi Implicit Methods for Pressure Linked Equations (SIMPLE) solver.

Guidelines to run on Unix-based system:

$ make
# It will make the run file

$ ./run
# This will start the simulation

If gnuplot is installed the results (while the simulation is going on) may be viewed by using the following commands on the command line:

$ gnuplot
gnuplot> p 'Channelu' w image,'Extensionu' w image

# The above plots the x-velocity contour-plot in the domain
# In order to view the y-velocity or pressure distribution in the domain replace 'u' in 'Channelu' and 'Extensionu' with 'v' or 'p' respectively.

In order to check the residuals, use the following command in gnuplot.
gnuplot> p './res_u' w l,'./res_v' w l,'./res_p' w l

In order to remove all case files written during a previous run of the executable 'run' and use it again, do the following:

$ make resetf
# This will delete all case file written during a previous run

Contact:
In case of error or difficulty, email anupambiswas85@gmail.com.
