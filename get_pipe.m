function elempos = get_pipe(rad,len,azi,xshift,yshift,depth,ellength,thick,epsr_pipe,epsr_below,dohalfpipe)
% Create a pipe or a halfpipe constructed of square planar elements
%
% Input arguments: 
% - rad: Radius of pipe [m]
% - len: Length of pipe [m]
% - azi: Azimuthal rotation of pipe [rad]
% - xshift: Shift of pipe in positive x-direction [m]
% - yshift: Shift of pipe in positive y-direction [m]
% - depth: Depth of pipe in negative z-direction [m]
% - ellength: Length of scattering element [m]
% - thick: Thickness of plane [m] 
%          (If the thicness is zero, the standard reflection 
%          coefficient is chosen automatically, meaning that the 
%          layer is infinitely thick.)
% epsr_pipe: Relative electric permittivity of pipe [-]
% epsr_below: Relative electric permittivity of material below pipe [-]
%             If the pipe has no thickness, it becomes infinitely thick. 
%             In that case, there is no material below the pipe. This value 
%             is then disregarded, but still needs to be provided. 
% - dohalfpipe: Do halfpipe open upwards (1), open downwards (-1) or fullpipe (0)
%
% Output: 
% - elempos: 9 x nel matrix, where nel is the number of square scattering 
%            elements. The nine rows are the x-, y-, and z-coordinate of 
%            the center of the squares as well as the azimuth phi, the dip
%            theta, the surface area and the thickness of the scattering 
%            element. The last two entries are the relative electric
%            permittivity of the pipe and of the material below. 
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
% 2022 (submitted).

circ = 2*pi*rad;
num_of_el = round(circ/ellength);
dalpha = 2*pi/num_of_el; % angular interval to have an element [rad]

yel = len/ellength;
start = (len-ellength)/2;
yvec=linspace(-start,start,yel);

elempos = zeros(9,num_of_el*length(yvec));
for iy = 1:length(yvec)
    for iel=1:num_of_el
        elempos(1,iel+(iy-1)*num_of_el) = -rad*cos((iel-1)*dalpha+pi/2);
        elempos(2,iel+(iy-1)*num_of_el) = yvec(iy);
        elempos(3,iel+(iy-1)*num_of_el) = rad*sin((iel-1)*dalpha+pi/2);
        elempos(4,iel+(iy-1)*num_of_el) = azi;
        elempos(5,iel+(iy-1)*num_of_el) = (iel-1)*dalpha;
        elempos(6,iel+(iy-1)*num_of_el) = ellength^2;
        elempos(7,iel+(iy-1)*num_of_el) = thick;
        elempos(8,iel+(iy-1)*num_of_el) = epsr_pipe;
        elempos(9,iel+(iy-1)*num_of_el) = epsr_below;
    end
end

if dohalfpipe == 1
    icount = 1;
    for iel=1:size(elempos,2)
        if elempos(3,iel)<=0.0
            newelempos(:,icount) = elempos(:,iel);
            icount = icount + 1;
        end
    end
    clear elempos
    elempos = newelempos;
    clear newelempos
elseif dohalfpipe == -1
    icount = 1;
    for iel=1:size(elempos,2)
        if elempos(3,iel)>=0.0
            newelempos(:,icount) = elempos(:,iel);
            icount = icount + 1;
        end
    end
    clear elempos
    elempos = newelempos;
    clear newelempos
end

% Rotate the pipe around the z-axis
azirot = [ cos(azi),-sin(azi), 0;...
           sin(azi), cos(azi), 0;...
           0,        0,        1];
elempos(1:3,:) = azirot*elempos(1:3,:);

elempos(1,:) = elempos(1,:) + xshift;
elempos(2,:) = elempos(2,:) + yshift;
elempos(3,:) = elempos(3,:) - depth;
