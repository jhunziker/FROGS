function [plotx,ploty,plotz] = get_corners(elempos,ellength)
% Calculates the position of the four corners of a square element oriented
% arbitrarilly in 3D.
% The input elempos contains the x-, y-, and z-coordinate of the center of
% the square as well as the azimuth phi and the dip theta. 
% The second input argument is the length of one side of the square. 
% The output are three 4 x 1 vectors containing the x-, y-, and
% z-coordinates of the four corner points of the square. 
% 
% This file is part of FROGS. FROGS is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the Free Software Foundation, 
% either version 3 of the License, or (at your option) any later version. FROGS is distributed 
% in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
% of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
% more details. You should have received a copy of the GNU General Public License along with FROGS. 
% If not, see <https://www.gnu.org/licenses/>. 
% 
% When refering to FROGS in any publication, please cite the paper "Fast 3D ground penetrating 
% radar simulations for glaciers" by J. Hunziker, E.C. Slob and J. Irving, Computers & Geosciences, 
% Vol. 173, 2023: https://doi.org/10.1016/j.cageo.2023.105320

c2m = ellength/2; % Distance from center point to middle of side

thetatemp = elempos(5)+pi/2;
phitemp = elempos(4);

% Calculate the position of the midpoint of all sides of the square
midlow(1) = c2m*sin(thetatemp)*cos(phitemp);
midlow(2) = c2m*sin(thetatemp)*sin(phitemp);
midlow(3) = c2m*cos(thetatemp);
midtop(1) = c2m*sin(thetatemp-pi)*cos(phitemp);
midtop(2) = c2m*sin(thetatemp-pi)*sin(phitemp);
midtop(3) = c2m*cos(thetatemp-pi);
% Before the shift, the midright and midleft points are in the
% xy-plane (i.e., z = 0);
midright(1) = c2m*cos(phitemp+pi/2);
midright(2) = c2m*sin(phitemp+pi/2);
midright(3) = 0.0;
midleft(1) = c2m*cos(phitemp-pi/2);
midleft(2) = c2m*sin(phitemp-pi/2);
midleft(3) = 0.0;

% Calculate the corners
botright = midlow+midright;
botleft = midlow+midleft;
topright = midtop+midright;
topleft = midtop+midleft;

plotx = elempos(1)+[botright(1);botleft(1);topleft(1);topright(1)];
ploty = elempos(2)+[botright(2);botleft(2);topleft(2);topright(2)];
plotz = elempos(3)+[botright(3);botleft(3);topleft(3);topright(3)];
