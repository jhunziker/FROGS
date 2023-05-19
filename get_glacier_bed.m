function elempos = get_glacier_bed(xmin,xmax,ymin,ymax,points,ellength,thick,epsr_sed,epsr_below,doplot)
% Define the glacier bed with a few points, interpolate to a surface and
% approximate that surface with square elements. 
%
% Note: If your glacier bed is just an inclined plane, use get_plane to
%       generate it. 
%
% Input arguments: 
% - xmin, xmax, ymin, ymax: Maximum horizontal extent of glacier bed [m]
% - points: x-, y-, and z-coordinates of points building the glacier bed [m]
%   Note, that the four corners must be defined. 
% - ellength: Length of scattering element [m]
% - thick: Thickness of sediments [m] 
%          (If the thickness is zero, the standard reflection 
%          coefficient is chosen automatically, meaning that the 
%          layer is infinitely thick.)
% - epsr_sed: Relative electric permittivity of sediments [-]
% - epsr_below: Relative electric permittivity of material below sediments [-]
%             If the sediments have no thickness, they become infinitely thick. 
%             In that case, there is no material below the sediments. This value 
%             is then disregarded, but still needs to be provided. 
%
% Output: 
% - elempos: 9 x nel matrix, where nel is the number of scattering 
%            elements. The nine rows are the x-, y-, and z-coordinate of 
%            the center of the elements as well as the azimuth phi, the dip
%            theta, the surface area and the thickness of the scattering 
%            element. The last two entries are the relative electric
%            permittivity of the sediments and of the material below. 
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

xvec = linspace(xmin+ellength/2,xmax-ellength/2,round((xmax-xmin)/ellength));
yvec = linspace(ymin+ellength/2,ymax-ellength/2,round((ymax-ymin)/ellength));
[xgrid_center,ygrid_center] = ndgrid(xvec,yvec);
[xgrid_A,ygrid_A] = ndgrid(xvec-ellength/2,yvec-ellength/2);
[xgrid_B,ygrid_B] = ndgrid(xvec+ellength/2,yvec-ellength/2);
[xgrid_C,ygrid_C] = ndgrid(xvec+ellength/2,yvec+ellength/2);
[xgrid_D,ygrid_D] = ndgrid(xvec-ellength/2,yvec+ellength/2);
gbed_center = griddata(points(:,1),points(:,2),points(:,3),xgrid_center,ygrid_center,'cubic');
gbed_A = griddata(points(:,1),points(:,2),points(:,3),xgrid_A,ygrid_A,'cubic');
gbed_B = griddata(points(:,1),points(:,2),points(:,3),xgrid_B,ygrid_B,'cubic');
gbed_C = griddata(points(:,1),points(:,2),points(:,3),xgrid_C,ygrid_C,'cubic');
gbed_D = griddata(points(:,1),points(:,2),points(:,3),xgrid_D,ygrid_D,'cubic');

elempos = zeros(5,length(xvec)*length(yvec));
for iy = 1:length(yvec)
    for ix=1:length(xvec)
        % Calculate the x-component of the gradient
        if ix==length(xvec)
            gradx = (gbed_center(ix,iy)-gbed_center(ix-1,iy))/ellength;
        elseif ix==1
            gradx = (gbed_center(ix+1,iy)-gbed_center(ix,iy))/ellength;
        else
            gradx = (gbed_center(ix+1,iy)-gbed_center(ix-1,iy))/(2*ellength);
        end
        % Calculate the y-component of the gradient
        if iy==length(yvec)
            grady = (gbed_center(ix,iy)-gbed_center(ix,iy-1))/ellength;
        elseif iy==1
            grady = (gbed_center(ix,iy+1)-gbed_center(ix,iy))/ellength;
        else
            grady = (gbed_center(ix,iy+1)-gbed_center(ix,iy-1))/(2*ellength);
        end
        % Calculate the value of the gradient
        gradval = sqrt(gradx^2+grady^2);
        
        elempos(1,ix+(iy-1)*length(xvec)) = xvec(ix);
        elempos(2,ix+(iy-1)*length(xvec)) = yvec(iy);
        elempos(3,ix+(iy-1)*length(xvec)) = gbed_center(ix,iy);
        % The angle between two vectors is
        % acos(dotprod(a,b)/(length(a)*length(b)))
        % Here, a points to where the azimuth is 0: a = [-1,0], and b is
        % the gradient. 
        % Note, the dot production only measures the angle between the two
        % vectors, but not if one vector is to the left or to the right of
        % the first one. Therefore, the sign of the y-component of the
        % gradient needs to be added. 
        if gradval<1e-6
            elempos(4,ix+(iy-1)*length(xvec)) = 0.0;
            elempos(5,ix+(iy-1)*length(xvec)) = 0.0;
        else
            elempos(4,ix+(iy-1)*length(xvec)) = -sign(grady)*acos(-1.0*gradx/gradval);
            elempos(5,ix+(iy-1)*length(xvec)) = atan(gradval);
        end
        % Calculate the area of the quadrangle, assuming that the corners
        % are at the altitude of the glacier bed. 
        % https://de.wikipedia.org/wiki/Viereck#Formeln
        tempa = sqrt((xgrid_A(ix,iy)-xgrid_B(ix,iy))^2 + (ygrid_A(ix,iy)-ygrid_B(ix,iy))^2 + (gbed_A(ix,iy)-gbed_B(ix,iy))^2);
        tempb = sqrt((xgrid_B(ix,iy)-xgrid_C(ix,iy))^2 + (ygrid_B(ix,iy)-ygrid_C(ix,iy))^2 + (gbed_B(ix,iy)-gbed_C(ix,iy))^2);
        tempc = sqrt((xgrid_C(ix,iy)-xgrid_D(ix,iy))^2 + (ygrid_C(ix,iy)-ygrid_D(ix,iy))^2 + (gbed_C(ix,iy)-gbed_D(ix,iy))^2);
        tempd = sqrt((xgrid_D(ix,iy)-xgrid_A(ix,iy))^2 + (ygrid_D(ix,iy)-ygrid_A(ix,iy))^2 + (gbed_D(ix,iy)-gbed_A(ix,iy))^2);
        tempe = sqrt((xgrid_A(ix,iy)-xgrid_C(ix,iy))^2 + (ygrid_A(ix,iy)-ygrid_C(ix,iy))^2 + (gbed_A(ix,iy)-gbed_C(ix,iy))^2);
        tempf = sqrt((xgrid_B(ix,iy)-xgrid_D(ix,iy))^2 + (ygrid_B(ix,iy)-ygrid_D(ix,iy))^2 + (gbed_B(ix,iy)-gbed_D(ix,iy))^2);
        elempos(6,ix+(iy-1)*length(xvec)) = sqrt(4*tempe^2*tempf^2-(tempb^2+tempd^2-tempa^2-tempc^2)^2)/4;
        % Put thickness in the elempos-matrix
        elempos(7,ix+(iy-1)*length(xvec)) = thick;
        % Put the relative-permittivity values in the elempos-matrix
        elempos(8,ix+(iy-1)*length(xvec)) = epsr_sed;
        elempos(9,ix+(iy-1)*length(xvec)) = epsr_below;
    end
end


if doplot
    plotx = [xgrid_A(:),xgrid_B(:),xgrid_C(:),xgrid_D(:)].';
    ploty = [ygrid_A(:),ygrid_B(:),ygrid_C(:),ygrid_D(:)].';
    plotz = [gbed_A(:),gbed_B(:),gbed_C(:),gbed_D(:)].';
    
    figure(12);
    fill3(plotx,ploty,plotz,elempos(6,:))
    grid on
    axis equal
    xlabel('x [m]')
    ylabel('y [m]')
    zlabel('z [m]')
    
%     figure;
%     h=surf(yvec,xvec,gbed_center,reshape(squeeze(elempos(6,:)),length(xvec),length(yvec)));
%     colormap('jet');
%     colorbar
%     set(h,'EdgeColor','none');
%     axis image
%     xlabel('Distance [m]')
%     ylabel('Distance [m]')
%     title('Topography of glacier bed with size of scattering elements')
%     set(gca,'XDir','reverse');
end
