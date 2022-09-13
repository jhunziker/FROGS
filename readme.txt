The FROGS package contains the following files: 
-----------------------------------------------

FROGS_Otemma.m: Runs one profile of the example on the Otemma glacier of the paper "Fast 3D ground penetrating radar simulations for glaciers" by J. Hunziker, E.C. Slob and J. Irving, Computers & Geosciences, 2022.

FROGS_validation.m: Runs the validation example of the above-mentioned paper. 

get_corners.m: Function that calculates the position of the corners of the square scattering elements. This is needed to plot the geometry. It is called by the main FROGS scripts. 

get_glacier_bed.m: Function to calculate the position and orientation of the scattering elements to approximate the glacier bed. As input a set of arbitrarily distributed points that define the surface of the glacier bed is needed. 

get_pipe.m: Function that creates a segment of a pipe or a halfpipe built from scattering elements. Originally created to model channels. 

get_plane.m: Function that creates a plane built from scattering elements. 

GNU_General_Public_License.txt: The license under which FROGS can be used, modified and distributed. 

Otemma_chanmap_dx0.1_chanwidth4.mat: Dataset that contains information about a model of a channel at the glacier bed of the Otemma glacier. Derived from Egli, P.E., Irving, J., Lane, S.N., 2021. Characterization of subglacial marginal channels using 3-D analysis of high-density ground-penetrating radar data. Journal of Glaciology 67, 759-772. 

Otemma_points_err1mat.mat: Set of points that define the glacier bed of the Otemma glacier with a precision of 1 m. 

radiation_pattern.m: Function that calculates the radiation pattern of the GPR antenna as a function of orientation and distance from the antenna. 

readme.txt: This file that you are currently looking at. 

sincinterpol.m: Function that performs sinc-interpolation. Used to reconstruct the subsampled frequency vector of the modeled GPR response in the frequency domain. 

trace_ES.txt: GPR-trace created with the semi-analytical code. Used to validate FROGS. 
