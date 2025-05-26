README: AzRM MATLAB PACKAGE
In this document, we provide a detailed description of the AzRM MATLAB package for azimuthal reflectivity method. 
The AzRM package includes functions to implement the azimuthal reflectivity method for two kinds of three-layered model, VTI and OA media.
The explanation of the function used in the AzRM MATLAB package are as follows:
	loadmodel_OA/VTI is as the input model data into the main function main AzRM_OA/VTI.
One row in loadmodel_VTI is as the following format:
thickness(km) vp(km/s) vs(km/s) density(g/cm^3) epsilon delta gamma;
One row in loadmodel_OA is as the following format:
thickness(km) vp(km/s) vs(km/s) density(g/cm^3) epsilon delta gamma e (fracture density);

	parameter_OA/VTI is another input data as some necessary parameters into the main function main AzRM_OA/VTI.
parameter_VTI includes:
T_max: maximum modeling time (sec);
delta_T: time interval (msec);
f1: initial frequency (Hz) to model; 
f2: final frequency (Hz) to model;
theta1: initial incidence (degree) to model;
theta2: final incidence (degree) to model;
fdom: dominant frequency (Hz) of the wavelet.

parameter_OA includes:
T_max: maximum modeling time (sec);
delta_T: time interval (msec);
f1: initial frequency (Hz) to model; 
f2: final frequency (Hz) to model;
theta1: initial incidence (degree) to model;
theta2: final incidence (degree) to model;
phi1: initial azimuth (degree) to model;
phi2: final azimuth (degree) to model;
fdom: dominant frequency (Hz) of the wavelet;
phi0: fracture orientation (degree), a single value.

	main_AzRM_OA/VTI has two functions: one is to display the calculated seismic gather; the other is to perform time-depth conversion, to locate the interface position in the time-domain seismic data, and then to observe AVA/AVAZ response.

	Stiffness_matrix_OA/VTI computes the stiffness matrix of anisotropic media using rock physics model and equivalent medium theory.

	AzRM_OA/VTI computes the overall reflection coefficient matrix for each frequency, each azimuth, and all incidence angles in the stratified model. Then, it multiplies with the Ricker wavelet in the frequency domain and uses 1D FFT algorithm to compute the azimuthal angle gather.

	compute_eig_OA/VTI computes the system matrix, the eigenvalue matrix, and the eigenvector matrix for each layer and each horizontal slowness, which is the most important part of this toolbox.

	compute_Rpp_OA/VTI has two functions: First, it computes the reflection and transmission coefficient matrices for each layer. Second, it computes the overall reflection and transmission coefficient matrices using the iterative algorithm.

	Ruger_VTI and Graebner_VTI compute synthetic angle gather using the reflection coefficient approximation and the exact equation for VIT media, respectively.

	Ruger_OA computes azimuthal angle gather using the reflection coefficient approximation for OA media.

	plotseis_AzRM is a function to do a quick plot of the seismic trace gather, published by the CREWES Project at the Department of Geology and Geophysics of the University of Calgary, Calgary, Alberta, Canada (The copyright and ownership is jointly held by G.F. Margrave and the CREWES Project).
