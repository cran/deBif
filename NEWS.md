# deBif 0.1.8 (current, 12/1/2023)

* Corrected the checking whether variables had gone of scale during limit cycle continuation

* Updated the DESCRIPTION file as some of the required packages (e.g. shiny) are no longer available for the earlier R versions

* Built-in protection in bif1parplot.R against absence of "eigvals" from the curve list

* Introduced a function deBifReset() to reload the package in case of computational problems

* Implemented a routine that computes the multipliers of a limit cycle, relying on the built-in function solve() to solve for the inverse of a matrix and for the generalized eigenvalues. These multipliers are subsequently used to plot the limit cycle either with a solid line in case it is stable or with a dashed line when unstable

# deBif 0.1.7 (current, 1/18/2023)

* Changed the use of the function sprintf() to snprintf() in the underlying C code as sprintf() has been deprecated in macOS 13

* Corrected a bug (superfluous checking of the value of 'ycol' in case length(state) == 1) in phaseServerFuncs.R that prevented the preventing of all steady states

* Corrected a missing value TRUE/FALSE error that occurred when phaseplane() was started for a 1-dimensional model for the very first time 

* Corrected a bug in phasePlot2D.R that generated an error when selecting a variable with index > 2 in a model with more than 2 variables

* Added an argument 'verify_fa = FALSE' to various calls to shiny::icon() in phaseUI.R and bifUI.R to prevent error messages from the FontAwesome package

# deBif 0.1.6 (5/16/2022)

* Introduced a new numerical option "Initial step size" to control the initial step along a solution curve

* Changed the default grid dimension in the Portrait tab of phaseplane() to 8 and allowed this setting to be varied between 3 and 20

* Reduced the font sizes in the menus of both phaseplane() and bifurcation()

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




