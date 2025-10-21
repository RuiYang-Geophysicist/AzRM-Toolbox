# AzRM-Toolbox
Azimuthal reflectivity modeling (AzRM) tool
AzRM MATLAB Package

This repository contains the MATLAB implementation of the Azimuthal Reflectivity Method (AzRM), supporting both VTI and orthorhombic (OA) media. It is developed for modeling azimuthally anisotropic seismic responses, and supports AVA/AVAZ analysis and synthetic gather generation.

📁 Contents

The package includes the following components:

loadmodel_OA.txt / loadmodel_VTI.txt: Input model files. Each row defines one layer:

VTI model format: thickness(km) vp(km/s) vs(km/s) density(g/cm³) epsilon delta gamma

OA model format: thickness(km) vp(km/s) vs(km/s) density(g/cm³) epsilon delta gamma e (fracture density)

(For the Orthorhombic case, I use Schoenberg's LSD model to establish the relationship between fracture compliance and stiffness coefficients. The fracture density obtained is based on the equivalent medium theory, which considers all layered media as a whole.)

parameters_OA.txt / parameters_VTI.txt: Input parameter files specifying modeling and wavelet settings.

VTI parameters:

T_max: maximum modeling time (s)

delta_T: time interval (s)

f1/f2: frequency range (Hz)

theta1/theta2: incident angle range (°)

fdom: dominant frequency of Ricker wavelet (Hz)

OA parameters (in addition to above):

phi1 / phi2: azimuthal angle range (°)

phi0: fracture orientation (°)

🔧 Function Descriptions

Main Functions

main_AzRM_VTI.m / main_AzRM_OA.m: Main drivers. These perform: Synthetic seismic gather computation, Time-depth conversion, AVA/AVAZ extraction and plotting, Core Algorithms.

Stiffness_matrix_VTI.m / Stiffness_matrix_OA.m: Compute the stiffness tensor of VTI or OA media using equivalent medium theory.

AzRM_VTI.m / AzRM_OA.m: Core Azimuthal Reflectivity Method. For each frequency and angle: Computes PP reflection coefficients, Applies Ricker wavelet in frequency domain, Performs inverse FFT to time domain.

compute_eig_VTI.m / compute_eig_OA.m: Compute eigenvalues and eigenvectors of elastodynamic system matrices, per layer and slowness. This is the core of anisotropic wave modeling.

compute_Rpp_VTI.m / compute_Rpp_OA.m: Compute reflection and transmission matrices using the Fryer-Frazer iterative method.

Comparison Methods: 

Ruger_VTI.m: Approximated reflection coefficient modeling for VTI media (Rüger, 2002)

Graebner_VTI.m: Exact reflection coefficient computation (Graebner, 1992)

Ruger_OA.m: Approximate AVAZ modeling for OA media using Rüger + Chen formulation

Plotting: 

plotseis_AzRM.m:Visualization utility to display angle gathers.(Provided by the CREWES Project, University of Calgary. Ownership: G.F. Margrave & CREWES)

📚 References

Yang R, Chen H, Guo Z, et al. An effective azimuthal reflectivity modeling (AzRM) tool for generating seismic data in anisotropic shale reservoirs[J]. Geophysics, 2025, 90(5): 1-76.
Fryer G J, Frazer L N. Seismic waves in stratified anisotropic media[J]. Geophysical Journal International, 1984, 78(3): 691-710.
Fryer G J, Frazer L N. Seismic waves in stratified anisotropic media—II. Elastodynamic eigensolutions for some anisotropic systems[J]. Geophysical Journal International, 1987, 91(1): 73-101.
Schoenberg M, Helbig K. Orthorhombic media: Modeling elastic wave behavior in a vertically fractured earth[J]. Geophysics, 1997, 62(6): 1954-1974.

📌 Notes

This code is intended for research and academic use. Please verify all outputs carefully.The author does not guarantee correctness or numerical stability in all edge cases.

