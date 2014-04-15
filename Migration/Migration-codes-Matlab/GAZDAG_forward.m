function D=GAZDAG_forward(Model,dx,z,dt,tmax,v,padt,padx)
% Author Nasser Kazemi
%   Copyright (C) 2014 University of Alberta
%   
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%   
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%   
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%Contact information:
%(electronic):
%Nasser Kazemi
%kazemino@ualberta.ca

%(paper):
%Nasser Kazemi
%Department of Physics
%CCIS,4-183
%University of Alberta
%Edmonton Alberta CANADA
%T6G 2E1

% This program implies post stack GAZDAG migration
% data is post stack data
% dt is time sampling interval
% h is a vector of CMP locations
% z0 is the first depth to migrate
% deltaz is the depth interval og migrated data
% zmax is the maximum depth to migrate
% v is velocity file v=v(z)
%fmin is minimum frequency to migrate
% fmax is maximum frequency to migrate
% padt   is padding with zero factor in time direction
% padx   is padding with zero factor in spatial direction


[m,n]=size(Model);
tt=floor(tmax/dt)+1;
NNx=padx*(2^nextpow2(n));% padd with zeros in horizontal wavenumber direction
NNt=padt*2^nextpow2(tt);
w=linspace(0, (pi-2*pi/NNt)/dt, NNt/2); 
Kx = linspace(-pi/dx, (pi-2*pi/NNx)/dx, NNx); % Horizontal wavenumber vector
%Kz=linspace(0, (pi-2*pi/NNz)/dz, NNz/2+1);            % Frequency vector (Just choose one half of the frequencies due to symmetry properties )
M=fftshift(fft(Model,NNx,2),2); %  apply 2Dfft with zero frequency in the center
P=M;               % Grab only one half of the data in f-k domain
i=sqrt(-1);
dz=z(2)-z(1);   
Dz=zeros(NNt,NNx);              % temporary file for each depth interval
for l=1:length(z)
for j=1:length(w)  
    for mm=1:length(Kx)         % loop over horizontal wavenumber
        kkz=(w(j)).^2/(v(end+1-l).^2)-Kx(mm).^2; % calculate depth shifting wavenumber
    Kz=-sqrt(kkz);
    if (kkz <0 || w(j)==0 || Kx(mm)==0) % don't do any thing for the non-physical part of Kz
        Dz(NNt/2+j,mm)=Dz(NNt/2+j,mm);
    else
        Dz(NNt/2+j,mm)=(Dz(NNt/2+j,mm)+P(end+1-l,mm)).*exp(i*Kz*dz); % apply phase shift to the physical part of wave field

    end
    end
end
end
        
    
Dz(2:NNt/2,:)=conj(Dz(end:-1:NNt/2+2,end:-1:1)); % complete temporary file based on symmetry properties for real valued signals
r=real(ifft2(ifftshift(Dz)));                          % apply inverse fft from Kx---> x 
D=r(1:tt,1:n);                          % truncate de- migrated image




    


