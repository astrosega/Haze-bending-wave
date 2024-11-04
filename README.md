Code use for the Sega et al 2024 ICARUS, Volume 413, 1 May 2024, 115987.

The main codes are wake_rot24.pro for the dynamics of self-gravity wakes, haze_creation.pro for the dynamics of the haze, and saturn_func.pro for the ray-tracing. 
I have included warpers for the main codes that directly plot each figure in the paper, named figure1.pro, figure2.pro ...

INSTALATION

The ray-tracing code expects the "magnus" folder, phases1.sav, and stars_robust.sav in the root directory of you IDL instantiation (by default your User or Home folder).
The "magnus" folder contains the UVIS data used, described in Table 2 of Sega et al 2024. phases1.sav contains the observed phase of the wave shown in the same table.
stars_robust.sav contains geometry parameters for each occultation.

EXTERNAL RESOURCES

A stand-alone from-first-principles exposition of linear bending wave theory is available in the Appendix C of my PhD thesis 
https://scholar.colorado.edu/concern/graduate_thesis_or_dissertations/8910jw225

