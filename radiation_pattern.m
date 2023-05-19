function [E1,E2,E3] = radiation_pattern(I,dz,eta,n,theta_c,k1,k2,theta,phi,r,doradpat)
% Calculate radiation pattern of antennas at an interface between two halfspaces. 
% References: 
% - Engheta, N., Papas, C.H., Elachi, C., 1982. Radiation patterns of interfacial dipole antennas. 
%   Radio Science 17, 1557-1566.
% - Arcone, S.A., 1995. Numerical studies of the radiation patterns of resistively loaded dipoles. 
%   Journal of Applied Geophysics 33, 39-52.
% 
% Input arguments: 
% - I: the current amplitude [A]
% - dz: the dipole element length [m]
% - eta: impedance of free space [Ohm]
% - n: Index of refraction below the interface [-]
% - theta_c: Critical angle [rad] 
% - k1: wavenumber in upper halfspace [m^(-1)]
% - k2: wavenumber in lower halfspace [m^(-1)]
% - theta: Dip angle between antenna and scattering element [rad] 
% - phi: Azimuth angle between antenna and scattering element [rad]
% - r: One-way radial distance between antenna and scattering element [m] 
% - doradpat: Switch to calculate the analytical radiation pattern (1) or to simulate
%             equal radiation in all directions (0).  
% 
% Output: 
% - E1, E2, E3: Electrical field components in spherical coordinates [V/m] 
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

% Calculate radiation pattern for given angles and distance
K1 = 1i*I*dz*k1.*eta.*exp(1i*(k1*r))./(2*pi*r); % upper halfspace
K2 = 1i*I*dz*k2.*eta.*exp(1i*(k2*r))./(2*pi*r); % lower halfspace

if doradpat==1 % Analytical radiation pattern in spherical coordinates
    if theta < pi/2 % The halfspace of air
        frac1 = (cos(theta))^2/(cos(theta)+sqrt(n^2-(sin(theta))^2));
        frac2 = (sin(theta))^2*cos(theta)*(cos(theta)-sqrt(n^2-(sin(theta))^2))/(n^2*cos(theta)+sqrt(n^2-(sin(theta))^2));
        E1 = K1*cos(phi)/n*(frac1-frac2);
        E2 = -K1/n*cos(theta)*sin(phi)/(cos(theta)+sqrt(n^2-(sin(theta))^2));
    else % The lower halfspace
        if theta >= pi-theta_c
            frac1 = (sqrt(1-n^2*(sin(theta))^2) + n*cos(theta))/(n*sqrt(1-n^2*(sin(theta))^2) - cos(theta));
            frac2 = ((cos(theta))^2)/(sqrt(1-n^2*(sin(theta))^2) - n*cos(theta));
            E1 = K2*cos(phi)*((sin(theta))^2*cos(theta)*frac1-frac2);
            E2 = K2*(cos(theta)*sin(phi)/(sqrt(1-n^2*(sin(theta))^2)-n*cos(theta)));
        else % theta <= pi-theta_c
            frac1 = (sqrt(n^2*(sin(theta))^2-1)-1i*n*cos(theta))/(n*sqrt(n^2*(sin(theta))^2-1)+1i*cos(theta));
            frac2 = (cos(theta))^2/(sqrt(n^2*(sin(theta))^2-1)+1i*n*cos(theta));
            E1 = K2*cos(phi)*((sin(theta))^2*cos(theta)*frac1+1i*frac2);
            E2 = -1i*K2*(cos(theta)*sin(phi)/(sqrt(n^2*(sin(theta))^2-1)+1i*n*cos(theta)));
        end
    end
else % Radially equal radiation pattern in spherical coordinates
    if theta < pi/2 % The halfspace of air
        E1 = K1;
        E2 = K1;
    else % The lower halfspace
        E1 = K2;
        E2 = K2;
    end
end
E3 = zeros(size(E1));
