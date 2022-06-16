function elempos = get_plane(xlen,ylen,azi,dip,xshift,yshift,depth,ellength,thick,epsr_plane,epsr_below)
% Create a plane constructed of square planar elements
%
% Input arguments: 
% - xlen: Length of plane in x-direction [m]
% - ylen: Length of plane in y-direction [m]
% - azi: Azimuthal rotation of plane [rad]
% - dip: Dip of plane [rad]
% - xshift: Shift of plane in positive x-direction [m]
% - yshift: Shift of plane in positive y-direction [m]
% - depth: Depth of plane in negative z-direction [m]
% - ellength: Length of scattering element [m]
% - thick: Thickness of plane [m] 
%          (If the thicness is zero, the standard reflection 
%          coefficient is chosen automatically, meaning that the 
%          layer is infinitely thick.)
% - epsr_plane: Relative electric permittivity of plane [-]
% - epsr_below: Relative electric permittivity of material below plane [-]
%               If the plane has no thickness, it becomes infinitely thick. 
%               In that case, there is no material below the plane. This 
%               value is then disregarded, but still needs to be provided. 
%
% Output: 
% - elempos: 9 x nel matrix, where nel is the number of square scattering 
%            elements. The nine rows are the x-, y-, and z-coordinate of 
%            the center of the squares as well as the azimuth phi, the dip
%            theta, the surface area and the thickness of the scattering 
%            element. The last two entries are the relative electric
%            permittivity of the plane and of the material below. 
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

if abs(xlen/ellength-round(xlen/ellength))>1e-6
    error('The length of the plane in x-direction is not an integer multiple of the element length.')
end
if abs(ylen/ellength-round(ylen/ellength))>1e-6
    error('The length of the plane in y-direction is not an integer multiple of the element length.')
end

xel = round(xlen/ellength);
yel = round(ylen/ellength);
startx = (xlen-ellength)/2;
starty = (ylen-ellength)/2;
xvec=linspace(-startx,startx,xel);
yvec=linspace(-starty,starty,yel);
elempos(1,:) = repmat(xvec,1,yel); % x-pos of center of element [m]
elempos(2,:) = reshape(repmat(yvec',1,xel)',1,xel*yel);% y-pos of center of element [m]
elempos(3,:) = zeros(1,length(elempos(1,:))); % z-pos of center of element [m]
elempos(4,:) =  azi; % azimuth of normal vector of element from x-axis [rad]
elempos(5,:) =  dip; % dip of normal vector of element from vertical [rad]
elempos(6,:) = ellength^2; % area of element [m^2]
elempos(7,:) = thick; % thickness of layer [m]
elempos(8,:) = epsr_plane; % Relative electric permittivity of plane [-]
elempos(9,:) = epsr_below; % Relative electric permittivity of material below plane [-]

% Rotate the plane around the y-axis (dip)
diprot = [ cos(dip), 0, sin(dip);...
           0,        1, 0       ;...
          -sin(dip), 0, cos(dip)];
elempos(1:3,:) = diprot*elempos(1:3,:);

% Rotate the plane around the z-axis (azimuth)
azirot = [ cos(azi),-sin(azi), 0;...
           sin(azi), cos(azi), 0;...
           0,        0,        1];
elempos(1:3,:) = azirot*elempos(1:3,:);

elempos(1,:) = elempos(1,:) + xshift;
elempos(2,:) = elempos(2,:) + yshift;
elempos(3,:) = elempos(3,:) - depth;
