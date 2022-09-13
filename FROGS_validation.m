% FROGS - Fast Radar-On-Glacier Simulations
% Simulate Ground Penetrating Radar data on glaciers in a fast and memory-efficient way.  
%
% Validation example of the paper "Fast 3D ground penetrating radar simulations for glaciers" 
% by J. Hunziker, E.C. Slob and J. Irving, Computers & Geosciences, 2022 (submitted). 
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

clear all; close all; clc; tic; 

% Switches
doradpat = 1; % (1) Take radiation pattern into account or not (0)

% Variables
I = 1; % Current Amplitude [A]
dz = 0.5; % Infinitesimal dipole element length [m]
epsr1 =  1.0; % Relative electric permittivity of the upper halfspace [-] (1.0 = air)
epsr2 =  3.2; % Relative electric permittivity of the lower halfspace [-] (3.2 = ice)
nf = 2*2048; % Number of frequencies [-]
maxfreq = 2e9; % Maximum frequency considered [Hz]
freq_c = 100*1e6; % Center frequency of antenna [Hz]
subfac = 2; % Factor for subsampling of frequency spectrum [-]
t0 = 0.0; % Time shift of wavelet
srcx = 0.0; % X-position of source antenna [m]
srcy = 0.0; % Y-position of source antenna [m]
srcz = 0.0; % Z-position of source antenna [m]
recx = 0.0; % X-position of receiver antenna [m]
recy = 0.0; % Y-position of receiver antenna [m]
recz = 0.0; % Z-position of receiver antenna [m]
phi_ant_src = 165.0/180.0*pi; % Rotation angle of source antenna [rad]
phi_ant_rec = 165.0/180.0*pi; % Rotation angle of receiver antenna [rad]
ellength = 0.5; % Length of scattering element [m]
critdist = 20.0; % Scattering elements need to be within this distance to be considered [m]
tapw = 10.0; % Width of slope of taper function
CoreNum = 32; % Amount of CPU's for parallel loop
fs = 9; % Fontsize
lw = 2; % Linewidth

% Open parallel pool
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool(CoreNum);
else
    disp('Parallel Computing for Multiple Cores is enabled.');
end

% Initialize the general elempos
elempos_gen = [];

% Plane of elements
planethick = 0.5; % Thickness of plane [m]
                  % (If the thickness is zero, the standard reflection
                  % coefficient is chosen automatically, meaning that the
                  % layer is infinitely thick.)
epsr_plane = 25.0; % Relative electric permittivity of plane [-] (81.0 = water)
epsr_below =  7.0; % Relative electric permittivity of material below plane [-] (3.2 = ice)
                   % If the plane has no thickness, it becomes infinitely
                   % thick. In that case, there is no material below the
                   % plane. This value is then disregarded.
planexshift = 0.0; % Shift the plane in positive x-direction [m]
planeyshift = 0.0; % Shift the plane in positive y-direction [m]
planedepth = 50.0; % Depth of plane in negative z-direction [m]
planexlen = 60.0; % Length of plane in x-direction [m] (Make this and the next value larger to avoid edge effects)
planeylen = 96.0; % Length of plane in y-direction [m]
planedip = 0.0/180.0*pi; % Dip rotation of plane [rad]
planeazi = 165.0/180.0*pi; % Azimuthal rotation of plane [rad]
elempos_gen = [elempos_gen,get_plane(planexlen,planeylen,planeazi,planedip,planexshift,planeyshift,planedepth,ellength,planethick,epsr_plane,epsr_below)];

% Parameters
c = 299792458; % Speed of light in vacuum [m/s]
mu0 = 4*pi*1e-7; % Magnetic permeability in vaccum [H/m]
eps0 = 1/(c^2*mu0); % Electric permittivity in vaccum [F/m]
eta = mu0*c; % impedance of free space [Ohm]
omega = 2*pi*linspace(0,maxfreq,nf); omega = omega(2:end); % Angular frequency [s^(-1)]
v = 1/sqrt(epsr2*eps0*mu0); % Wavespeed in the lower halfspace [m/s]
n = c/v; % Index of refraction below the interface [-]
theta_c = asin(1/n); % Critical angle [rad]

% Create wavelet
omega_c = 2*pi*freq_c; % Angular center frequency of antenna [s^(-1)]
wave = -(omega/omega_c).^2.*exp(-(omega/omega_c).^2).*exp(1i*omega*t0);

% Subsample omega for calculation
cutel = find(abs(wave)>1e-4*max(abs(wave)),1,'last');
omega_sub = linspace(1,omega(cutel),(2^nextpow2(cutel))/subfac);

% Calculate wavenumber on reduced frequency vector
k1 = omega_sub*sqrt(epsr1*eps0*mu0); % wavenumber in upper halfspace [m^(-1)]
k2 = omega_sub*sqrt(epsr2*eps0*mu0); % wavenumber in lower halfspace [m^(-1)]

% Plot geometry
plotx = zeros(4,size(elempos_gen,2));
ploty = zeros(size(plotx));
plotz = zeros(size(plotx));
for iel=1:size(elempos_gen,2) % Loop over the elements
    [plotx(:,iel),ploty(:,iel),plotz(:,iel)] = get_corners(elempos_gen(:,iel),ellength);
end

figure(1);
fill3(plotx,ploty,plotz,elempos_gen(6,:))
hold on
plotsrcx = [srcx+dz/2*cos(phi_ant_src),srcx-dz/2*cos(phi_ant_src)];
plotsrcy = [srcy+dz/2*sin(phi_ant_src),srcy-dz/2*sin(phi_ant_src)];
plot3(plotsrcx,plotsrcy,zeros(size(plotsrcx)),'b','Linewidth',lw)
plotrecx = [recx+dz/2*cos(phi_ant_rec),recx-dz/2*cos(phi_ant_rec)];
plotrecy = [recy+dz/2*sin(phi_ant_rec),recy-dz/2*sin(phi_ant_rec)];
plot3(plotrecx,plotrecy,zeros(size(plotrecx)),'--r','Linewidth',lw)
hold off
grid on
axis equal
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
set(gca,'XDir','reverse','YDir','reverse','Fontsize',fs)

isrc = 1; % DELETE THIS LINE WHEN INTRODUCING A LOOP OVER SOURCE POSITIONS

% Check that the critical distance is larger than the antenna separation
if sqrt((srcx(isrc)-recx(isrc))^2 + (srcy(isrc)-recy(isrc))^2)>critdist
    error('Antenna separation is larger than distance of scattering elements taken into account')
end

% Extract elements that are within a critical distance from the antennas. 
% Scattering elements further away are ignored. 
tempcounter = 0;
for iel=1:size(elempos_gen,2)
    for isrc=1:length(srcx)
        srcdist = sqrt((srcx(isrc)-elempos_gen(1,iel))^2+(srcy(isrc)-elempos_gen(2,iel))^2);
        recdist = sqrt((recx(isrc)-elempos_gen(1,iel))^2+(recy(isrc)-elempos_gen(2,iel))^2);
        if min(srcdist,recdist)<critdist
            tempcounter = tempcounter + 1;
            elempos(:,tempcounter) = elempos_gen(:,iel);
            break
        end
    end
end

% Create taper function for elements at the edge of the considered area
tapx = linspace(0.0,critdist,round(critdist/(0.6*ellength)));
tapy = ones(size(tapx));
slopstartel = find(tapx>=critdist-tapw,1,'first');
sloplen = length(tapx)-slopstartel+1;
slopx = linspace(0,pi,sloplen);
slopy = ((cos(slopx)+1.0)/2.0).^2;
tapy(slopstartel:end) = slopy;

% Source dependent parameters
r_src = sqrt((srcx(isrc)-elempos(1,:)).^2 + (srcy(isrc)-elempos(2,:)).^2 + (srcz(isrc)-elempos(3,:)).^2); % One-way radial distance from the source [m]
phi_src = pi-(atan2(elempos(2,:)-srcy(isrc),elempos(1,:)-srcx(isrc))); % Azimuth for radiation pattern [rad] (phi = 0 --> E-plane, phi = pi/2 --> H-plane)
if max(isnan(phi_src))==1
    error('Phi ist not a number.')
end
theta_src = acos((elempos(3,:)-srcz(isrc))./r_src); % Dip [rad] (theta = pi/2 --> horizontal direction)
if max(theta_src)>pi
    error('Theta is larger than pi.')
end

% Receiver dependent parameters
r_rec = sqrt((recx(isrc)-elempos(1,:)).^2 + (recy(isrc)-elempos(2,:)).^2 + (recz(isrc)-elempos(3,:)).^2); % One-way radial distance from the receiver [m]
phi_rec = pi-(atan2(elempos(2,:)-recy(isrc),elempos(1,:)-recx(isrc))); % Azimuth for radiation pattern [rad] (phi = 0 --> E-plane, phi = pi/2 --> H-plane)
if max(isnan(phi_rec))==1
    error('Phi ist not a number.')
end
theta_rec = acos((elempos(3,:)-recz(isrc))./r_rec); % Dip [rad] (theta = pi/2 --> horizontal direction)
if max(theta_rec)>pi
    error('Theta is larger than pi.')
end

% Calculate distance between source and receiver antennas
ant_dist = sqrt((srcx(isrc)-recx(isrc))^2 + (srcy(isrc)-recy(isrc))^2 + (srcz(isrc)-recz(isrc))^2);
ant_phi = abs(phi_ant_src-phi_ant_rec);

Erec = zeros(1,length(omega_sub));
% Parallel loop over scattering elements
parfor iel=1:size(elempos,2)
    % Radiation pattern at the source antenna
    [Etheta_src,Ephi_src,Er_src] = radiation_pattern(I,dz,eta,n,theta_c,k1,k2,theta_src(iel),phi_src(iel)+phi_ant_src,r_src(iel),doradpat);
    if (ant_dist>=1e-6 || ant_phi>=1e-6)
        [Etheta_rec,Ephi_rec,Er_rec] = radiation_pattern(I,dz,eta,n,theta_c,k1,k2,theta_rec(iel),phi_rec(iel)+phi_ant_rec,r_rec(iel),doradpat);
    else
        Etheta_rec = Etheta_src;
        Ephi_rec = Ephi_src;
        Er_rec = Er_src;
    end
    
    % Transform to Cartesian coordinates
    transy_src = [ cos(theta_src(iel)), 0, sin(theta_src(iel));...
        0,               1, 0;...
        -sin(theta_src(iel)), 0, cos(theta_src(iel))];
    transz_src = [ cos(-pi+phi_src(iel)),-sin(-pi+phi_src(iel)), 0;...
        sin(-pi+phi_src(iel)), cos(-pi+phi_src(iel)), 0;...
        0,                0,                1];    
    transy_rec = [ cos(theta_rec(iel)), 0, sin(theta_rec(iel));...
        0,               1, 0;...
        -sin(theta_rec(iel)), 0, cos(theta_rec(iel))];
    transz_rec = [ cos(-pi+phi_rec(iel)),-sin(-pi+phi_rec(iel)), 0;...
        sin(-pi+phi_rec(iel)), cos(-pi+phi_rec(iel)), 0;...
        0,                0,                1];    
    E_src = [-1,0,0;0,-1,0;0,0,1]*transz_src*transy_src*[Etheta_src;Ephi_src;Er_src];
    E_rec = [-1,0,0;0,-1,0;0,0,1]*transz_rec*transy_rec*[Etheta_rec;Ephi_rec;Er_rec];
    E = E_src.*E_rec;

    if elempos(7,iel)>-1e-6 % Reflection coefficient
        % Rotate the electric field into the coordinate system of the plane of the
        % reflector. As start, a horizontal plane with a normal pointing upward
        % parallel to the z-axis is considered.
        % First, rotate around the y-axis (dip of the normal vector from vertical)
        % Second, rotate around the z-axis (azimuth from the x-axis)
        % https://en.wikipedia.org/wiki/Rotation_matrix#Basic_rotations
        Ry = [ cos(elempos(5,iel)), 0, sin(elempos(5,iel));...
            0,                   1, 0;...
            -sin(elempos(5,iel)), 0, cos(elempos(5,iel))];
        Rz = [ cos(elempos(4,iel)),-sin(elempos(4,iel)), 0;...
            sin(elempos(4,iel)), cos(elempos(4,iel)), 0;...
            0,                   0,                   1];
        E = Rz*Ry*E;
        
        % Reflection coefficients for normal incidence according to Bradford and Deeds (2006)
        gamma2 = omega_sub*sqrt(epsr2*eps0*mu0); % vertical wavenumber in the ice
        gamma3 = omega_sub*sqrt(elempos(8,iel)*eps0*mu0); % vertical wavenumber in the scatterer
        gamma4 = omega_sub*sqrt(elempos(9,iel)*eps0*mu0); % vertical wavenumber below the scatterer
        
        if abs(elempos(7,iel))>1e-6
            % Bradford and Deeds (2006)
            RTE = (gamma2-gamma4-1i*(gamma2.*gamma4./gamma3-gamma3).*tan(gamma3*elempos(7,iel)))./(gamma2+gamma4-1i*(gamma2.*gamma4./gamma3+gamma3).*tan(gamma3*elempos(7,iel)));
            RTM = (gamma2*elempos(9,iel)*eps0-gamma4*epsr2*eps0-1i*(gamma2.*gamma4*elempos(8,iel)*eps0./gamma3-gamma3*epsr2*eps0*elempos(9,iel)/elempos(8,iel)).*tan(gamma3*elempos(7,iel)))./...;
                (gamma2*elempos(9,iel)*eps0+gamma4*epsr2*eps0-1i*(gamma2.*gamma4*elempos(8,iel)*eps0./gamma3+gamma3*epsr2*eps0*elempos(9,iel)/elempos(8,iel)).*tan(gamma3*elempos(7,iel)));
        else
            % Fresnel reflection coefficient
            RTE = (gamma2-gamma3)./(gamma2+gamma3);
            RTM = (gamma2-gamma3)./(gamma2+gamma3);
        end
        
        % Apply scattering through equation 9 of Shakas & Linde (2015)
        E = elempos(6,iel)*[RTE;RTE;RTM].*E;
        
        % Rotate the scattered electric field back to the original coordinate grid
        Ryinv = [ cos(-elempos(5,iel)), 0, sin(-elempos(5,iel));...
            0,                    1, 0;...
            -sin(-elempos(5,iel)), 0, cos(-elempos(5,iel))];
        Rzinv = [ cos(-elempos(4,iel)),-sin(-elempos(4,iel)), 0;...
            sin(-elempos(4,iel)), cos(-elempos(4,iel)), 0;...
            0,                    0,                    1];
        E = Ryinv*Rzinv*E;
        
    else % point scatterer
        epsr_pointscat = elempos(8,iel);
        scat_rad = elempos(7,iel); % Radius of point scatterer [m]
        vol = 4/3*pi*scat_rad^3; % Volume of a sphere [m^3]
        Aeps = log((epsr_pointscat)/epsr2)*vol; % Perturbation in the logarithmic permittivity
        % Unit vector pointing to the receiver antenna
        r_unit = [(recx(isrc)-elempos(1,iel))/r_rec(iel);(recy(isrc)-elempos(2,iel))/r_rec(iel);(0.0-elempos(3,iel))/r_rec(iel)];
        % Equation 19 of Saintenoy and Tarantola (2001)
        temp = zeros(size(E));
        for ifreq=1:length(omega_sub)
            temp(:,ifreq) = cross(Aeps*cross(r_unit,-(omega_sub(ifreq)^2)*E(:,ifreq)),r_unit);
        end
        E = temp/(4*pi*r_rec(iel)*v^2);
    end
   
    % Determine taper value
    srcdist = sqrt((srcx(isrc)-elempos(1,iel))^2+(srcy(isrc)-elempos(2,iel))^2);
    recdist = sqrt((recx(isrc)-elempos(1,iel))^2+(recy(isrc)-elempos(2,iel))^2);
    tapval = interp1(tapx,tapy,min(srcdist,recdist));
 
    % Sum up the fields from the various elements
    Erec = Erec + tapval*(E(1,:) + E(2,:) + E(3,:));
end

% Interpolate the data to the full omega-vector
Erec_fullfreq = zeros(size(omega));
[tempy,dtemp] = sincinterpol(omega_sub(2)-omega_sub(1),Erec,4,0);
tempx = linspace(1,length(tempy),length(tempy))*dtemp+(omega_sub(1)-dtemp);
Erec_fullfreq(1:cutel) = interp1(tempx,tempy,omega(1:cutel),'spline');

% Convolve the data with a wavelet
Erec_fullfreq = Erec_fullfreq.*wave;

% Complete the frequency vector.
% Remember the zero-frequency was not calculated and is here set to zero.
fErec = [zeros(1,1),Erec_fullfreq];
fErec = [zeros(1,1),conj(fliplr(fErec(2:end))),fErec];

% Transforming the data to the time domain
fvec = omega/(2*pi);
df = fvec(2)-fvec(1);
tErec(isrc,:)=length(fErec)*fftshift(ifft(fftshift(conj(fErec),2),[],2),2)*df;

% Set-up time vector
fvec = omega/(2*pi);
df = fvec(2)-fvec(1);
dt = 1/(2*nf*df);
tvec = linspace(-nf,nf-1,2*nf)*dt;

% Copy FROGS trace for comparison
trace_FROGS = squeeze(tErec(1,:));

% Load the trace from the semi-analytical code (Slob & Fokkema, 2002)
filenameES = 'trace_ES.txt';
ntES = 8192; % number of samples [-]
dtES = 0.3125; % dt [ns]
tvecES = linspace(0,ntES-1,ntES)*dtES;
fid = fopen(filenameES);
temp = fgetl(fid);
fclose(fid);
for it=1:ntES
    trace_ES(it) = str2num(temp((it-1)*25+1:it*25));
end

% Plot FROGS trace and semi-analytical trace
h=figure;
subplot(2,1,1)
plot(tvecES,trace_ES/max(trace_ES),'b','Linewidth',lw);
hold on
plot(tvec*1e9,trace_FROGS/max(trace_FROGS),'r--','Linewidth',lw);
hold off
xlim([550,675])
ylim([-1.2,1.2])
xlabel('Time [ns]')
legend('semi-analytical layered model code','FROGS')
set(gca,'Fontsize',fs)

toc
