The code in this folder generates radar rainfall ensembles with a method very similar to the REAL method, described by (Germann et al. 2009*).

The main differences are:

1) the weighting system adopted by Germann et al. is not used.
2) the normal score transform is applied to the errors
3) The decomposition of the covariance matrix is done with a singular value decomposition instead of a Cholewsky algorithm, eliminating the requirement for a positive definite covariance matrix.



What you need:

1) a .dat file containing the rain gauges coordinates, structured as two columns of easting and northing values in meters.
2) a .dat file containing the rain gauges values, structured with the following columns:
	1. year - 4 digits
	2. month - 2 digits
	3. day	- 2 digits
	4. hour - 2 digits
	5. minutes - 2 digits
	6. seconds -  - 2 digits
	7. Rain gauge #1 values in mm/h
	...
	n+7. Rain gauge #n values in mm/h
3) a .dat file containing the radar values averaged on 1 hour, in the locations corresponding to the rain gauges. The structure of the file is the same.
4) Radar files in nimrod format gridded at 1km in mm/h.
5) The .m files contained in this folder




How to run it:

To generate the ensembles, open the Germann.m file and complete the section "INFORMATION TO PROVIDE". Then run the script.




* Germann, U., Berenguer, M., Sempere-torres, D., & Zappa, M. (2009). REAL – Ensemble radar precipitation estimation for hydrology in mountainous region. Q. J. R. Meteorol. Soc, 135(February), 445–456. doi:10.1002/qj.375