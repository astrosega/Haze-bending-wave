Code use for the Sega et al 2024 ICARUS, Volume 413, 1 May 2024, 115987.

The main codes are wake_rot24.pro for the dynamics of self-gravity wakes, haze_creation.pro for the dynamics of the haze, and saturn_func.pro for the ray-tracing. 
I have included warpers for the main codes that directly plot each figure in the paper, named figure1.pro, figure2.pro ...

INSTALATION

There are three main pieces of code for this project, the orbital evolution of haze particles, the ray-tracing, and the dynamical simulation of rigid self-gravity wakes. The ray-tracing code uses haze_func.pro and haze_creation.pro from the haze package to compute the shape of the haze for each occultation, ortherwise the codes are independent. The .pro files in the root of the repository are functions and procedures that need to be compiled for multiple packages. Package-specific function and procedures lie in their respective folders.

The ray-tracing code expects the "magnus" folder, phases1.sav, and stars_robust.sav in the root directory of you IDL instantiation (by default your User or Home folder). You may also place the .pro files there so that IDL compiles the required functions and procedures automatically. The .pro files in the root of the repository are need in multiple packages

The "magnus" folder contains the UVIS data used, described in Table 2 of Sega et al 2024. phases1.sav contains the observed phase of the wave shown in the same table.
stars_robust.sav contains geometry parameters for each occultation.

EXTERNAL RESOURCES

A stand-alone from-first-principles exposition of linear bending wave theory is available in the Appendix C of my PhD thesis 
https://scholar.colorado.edu/concern/graduate_thesis_or_dissertations/8910jw225

