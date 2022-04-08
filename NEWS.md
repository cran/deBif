# deBif 0.1.5 (4/8/2022)

* Changed the calls to the Lapack routines dgetrf, dgecon and dgesvx to correctly pass string from C to Fortran following ‘Writing R Extensions’ §6.6.1

* Changed computation of derivatives in phaseplane() to solve problems with derivative functions that are not vectorizable. Derivatives are now only calculated when either parameter change or the plot options change.

* Stricter checking of input curves in bifurcation() and phaseplane() to prevent crashes when the parameter list has been changed

* Show steady state details in phaseplane for all tabs >= 3 and improved the checking on the steady state solutions in phaseplane()

* Corrected bug in 3D plotting of bifurcation graph that resulted in wrong axis labeling

* First point localisation in 2 parameter plot automatically swaps variables if unsuccessful

* Bugs solved for full report of computations

* Removed the on.exit() call at the start of bifurcation(), which generated a warning about the absence of a graphic context

* Changed the on-exit() calls in phaseplane to correctly plot the steady state points

* Wrapped the call to deSolve time integration method into tryCatch() to catch errors occurring during the time integration 

# deBif 0.1.0 (3/3/2022)

* First official version to be submitted to CRAN. The intensive computational parts (solving for solution points, solving for limit cycles) 
are now carried out by C-code. The package has been tested intensively by a group of BSc students during a 4-day class 

# deBif 0.0.5 (3/6/2021)

* Compared the computation of the limit cycle condition and the jacobian of the limit cycle condition with the same computations in Matcont. The jacobian and the limit cycle condition differ less than 1.0E-5 from each other, showing the correctness of the code 




