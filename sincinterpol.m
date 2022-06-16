function [fxout dx2] = sincinterpol(dx,fxin,nfac,dotrunc)
% Applies Sinc-Interpolation also known as Whittaker-Shannon interpolation 
% to a vector. Sinc-Interpolation interpolates the data within its 
% bandwidth. Thus, aliased data stays aliased. The improvement is only
% visual.
% 
% USAGE: [fxout dx2] = sincinterpol(dx,fxin,nfac,dotrunc)
% 
% INPUT:
% dx:       the spacing of the coordinates
% fxin:     the vector to be interpolated
% nfac:     the factor that specifies how many more points there will be in
%           the interpolated function
% dotrunc:  optional argument
%           0: The amount of datapoints in the sinc-interpolated trace is
%              the amount of datapoints before interpolation times nfac
%           1: The amount of datapoints in the sinc-interpolated trace is
%              such that no datapoint lies behind the last real datapoint.
% 
% OUTPUT:
% fxout:    interpolated vector
% dx2:      new spacing of the coordinates
% 
% REQUIREMENT:
% The amount of points in the vector needs to be 2^n, where n is any
% positive integer.  
% 
% EXAMPLE:
% clear all; close all;
% 
% xsize = 4; % Datapoints on which the sine function is computed
% dx = pi/2; % Spacing for computation of sine function
% nfac = 27; % Factor by which the amount of datapoints is increased by sinc interpolation
% dotrunc = 0; % (1) if the interpolated trace should exceed the last datapoint, (0) otherwise
% fs = 18; % Fontsize
% 
% % xvec = (-xsize/2:xsize/2-1)*dx; 
% xvec = (0:xsize-1)*dx; 
% fx = sin(xvec);
% [fx2 dx2] = sincinterpol(dx,fx,nfac,dotrunc);
% xsize2 = length(fx2); 
% % xvec2 = (-xsize2/2:xsize2/2-1)*dx2;
% xvec2 = (0:xsize2-1)*dx2;
% 
% figure;
% plot(xvec,fx,'ok'); hold on
% plot(xvec,fx,'b'); 
% plot(xvec2,fx2,'r'); hold off
% grid on; 
% legend('Datapoints','Linear','Sinc','Location','NorthEast');
% title(['Interpolation of a Sine Function (',num2str(xsize),' datapoints)'],'Fontsize',fs)
% set(gca,'Fontsize',fs)
% 
% AUTHOR:
% Juerg Hunziker, 25th of October 2013
% 
% UPDATES:
% 2022: Use with FROGS. 
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

if nargin<4
    dotrunc = 0;
end

xsize = length(fxin);
if 2^nextpow2(xsize)-xsize ~=0
    error('The size of fxin has to be 2^n, where n is any positive integer.')
end
dx2 = dx/nfac;
dkx = 2*pi/(xsize*dx);
temp = (xsize)*ifft(fftshift(fxin))*dx;
temp2 = zeros(1,xsize*nfac);
temp2(1:length(temp)/2) = temp(1:length(temp)/2);
temp2(xsize*nfac-length(temp)/2+1:xsize*nfac) = temp(length(temp)/2+1:length(temp));
fxout = fftshift(fft(temp2))*dkx/(2*pi);
if isreal(fxin)
    fxout = real(fxout);
end
if dotrunc == 1
    delpoints = nfac-1;
    fxout(end-delpoints+1:end) = [];
end
