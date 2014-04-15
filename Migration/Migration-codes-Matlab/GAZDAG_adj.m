function M=GAZDAG_adj(data,dt,h,z0,deltaz,zmax,v,fmin,fmax,padt,padx)
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
% deltaz is the depth interval of migrated data
% zmax is the maximum depth to migrate
% v is velocity file v=v(z)
%fmin is minimum frequency to migrate
% fmax is maximum frequency to migrate
% padt   is padding with zero factor in time direction
% padx   is padding with zero factor in spatial direction


[m,n]=size(data);
NNt=padt*(2^nextpow2(m));% padd with zeros in time direction
NNx=padx*(2^nextpow2(n));% padd with zeros in horizontal wavenumber direction
deltax=h(2)-h(1);     % spatial sampling interval
Kx = linspace(-pi/deltax, (pi-2*pi/NNx)/deltax, NNx); % Horizontal wavenumber vector
w=linspace(0, (pi-2*pi/NNt)/dt, NNt/2);            % Frequency vector (Just choose one half of the frequencies due to symmetry properties )
D=fftshift(fft2(data,NNt,NNx)); %  apply 2Dfft with zero frequency in the center
P=D(NNt/2:end,:);               % Grab only one half of the data in f-k domain
i=sqrt(-1);
ifmin=floor(fmin*dt*NNt)+1;     % calculate the index for minimum frequency to migrate
ifmax=floor(fmax*dt*NNt);       % calculate the index for maximum frequency to migrate  
izmax=floor((zmax-z0)/deltaz)+1;                
Mz=zeros(size(D));              % temporary file for each depth interval
M=zeros(izmax,n);           % migrated model
for l=1:izmax               % loop over depth
for j=ifmin:ifmax               % loop over frequency
    for mm=1:length(Kx)         % loop over horizontal wavenumber
        kkz=(w(j)).^2/(v(l).^2)-Kx(mm).^2; % calculate depth shifting wavenumber
    Kz=-sqrt(kkz);
    if (kkz <0 || w(j)==0 || Kx(mm)==0) % dont do any thing for the non-physical part of Kz
        Mz(NNt/2+j,mm)=P(j,mm);
    else
        Mz(NNt/2+j,mm)=P(j,mm).*exp(-i*Kz*(z0+(l-1)*deltaz)); % apply phase shift to the physical part of wave field

    end
    end
end
        
    
Mz(2:NNt/2,:)=conj(Mz(end:-1:NNt/2+2,end:-1:1)); % complete temporary file based on symmetry properties for real valued signals
r=ifft(ifftshift(Mz),[],2);                          % apply inverse fft from Kx---> x 
M(l,:)=real(sum(r(:,1:n)));                          % apply imaging condition :sum over frequency
end



    


