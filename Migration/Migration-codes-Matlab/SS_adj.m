function M=SS_adj(data,dt,h,z0,deltaz,zmax,v,fmin,fmax,padt,padx)
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

% This is a simple program for post stack Split Step migration
% data   is post stack data
% dt     is time sampling interval
% h      is a vector of CMP locations
% z0     is the first depth to migrate
% deltaz is the depth interval of migrated data
% zmax   is the maximum depth to migrate
% v      is velocity file v=v(x,z)
% nrefv   is the number of velocities
%fmin    is minimum frequency to migrate
% fmax   is maximum frequency to migrate
% padt   is padding with zero factor in time direction
% padx   is padding with zero factor in spatial direction

[m,n]=size(data);
NNt=padt*(2^nextpow2(m)); % padd with zeros in time direction 
NNx=padx*(2^nextpow2(n)); % padd with zeros in horizontal wavenumber direction
deltax=h(2)-h(1);      % calculate spatial sampling interval from header information

pi = 4.*atan(1.);
%% Define Horizontal Wavenumber vector using symmetry properties and
%% spatial sampling information
%% 
Kx=zeros(NNx,1);
       dkx = 2*pi/(NNx*deltax);

       Kx(1) = 0.;
       for i=2:NNx/2+1
       Kx(i) = Kx(i-1)+dkx;
       end
       Kx(NNx/2+2) = -(NNx/2-1)*dkx;
       for i=NNx/2+3:NNx
       Kx(i) = Kx(i-1) + dkx;
       end
w=zeros(NNt/2+1,1);
       dw = 2*pi/(NNt*dt);
%%.......................................................................%%       
%% Define the first half vector of frequency: Note we will use only first
%% half of the frequency vecotr and apply conjugate symmetry properties of
%% real signals
       w(1) = 0.;
       for i=2:NNt/2+1;
       w(i) = w(i-1)+dw;
       end
%%.......................................................................%%
D=fft(data,NNt);                  % Apply Fourier transform to the data (t-x)---> (f-x)
P=D(1:NNt/2+1,:);                 % Grab only one half of the data in f-x domain
i=sqrt(-1);
ifmin=floor(fmin*dt*NNt)+1;       % Calculate the index for minimum frequency to migrate
ifmax=floor(fmax*dt*NNt);         % Calculate the index for maximum frequency to migrate  
izmax=floor((zmax-z0)/deltaz)+1;                
Mz=zeros(NNx,1);              % Temporary file for each depth interval
M=zeros(izmax,n);                 % Migrated model
v=v/2;                            % we are using exploding reflector
%% Start PSPI migration%%
for j=ifmin:ifmax                 % loop over frequency
    P_s=P(j,:);                   % Grap one slice at a time 
for l=1:izmax                     % loop over depth
    P0=fft(P_s,NNx);              % Move frequency slice from x to Kx
    V=v(l,:);
    V_ave=sum(V)/length(V);        % Define average velocities for GAZDAG phase shift step 
    for mm=1:length(Kx)                        % loop over horizontal wavenumber
        kkz=(w(j)).^2/(V_ave.^2)-Kx(mm).^2; % calculate depth shifting wavenumber
    Kz=sqrt(kkz);
    if (kkz <0)                                % Don't do any thing for the non-physical part 
        Mz(mm,1)=0;
    else
        Mz(mm,1)=P0(mm).*exp(i*Kz*(deltaz));  % Apply phase shift to the physical part of wave field
    end
    end
     Mzx=ifft(Mz);                         % Move phase shifted versions to spatial domain
     Mzx=Mzx(1:n);
     
     
     
 %% Apply correction based on dp to the GAZDAG type shifted wavefields    
        dp = 1./V-1/V_ave;        
 
            P_s=Mzx.'.*exp(i*w(j)*dp*deltaz);
     M(l,:)=M(l,:)+2*real(P_s); % Apply imaging condition ( Note: Imaging condition is in the frequency domain) 
end
 if ( mod(j,10)==0)
     fprintf (' freq %d      max_freq %d \n',floor(j/(NNt*dt)),floor(ifmax/(NNt*dt)));
 end
end
     
     