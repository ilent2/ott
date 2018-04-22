function [z,kz] = axial_equilibrium(tmatrix,beam,z)
%AXIAL_EQUILIBRIUM find equilibrium position and stiffness along beam axis
%
% [z,kz] = AXIAL_EQUILIBRIUM(T,beam) attempts to locate the equilibrium
% position for the T-matrix T in beam starting with an initial
% guess at z = 0.
%
% [z,kz] = axial_equilibrium(..., initial_guess) specifies an initial guess.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% This function is not directly concerned with force/torque calculation
ott.warning('ott:axialEquilibrium:move', ...
    'This function will move in a future release');

ott.warning('internal');

if nargin < 4
    z = 0;
end

% Ensure T-matrix and beam are the same size
Nmax = max(tmatrix.Nmax, beam.Nmax);
tmatrix.Nmax = Nmax;
beam.Nmax = Nmax(1);

% Normalise the beam power
beam.power = 1.0;

%start with three points each side over a 1/8 wavelength:
zd=.25/(beam.k_medium/2/pi);
dz=zd/5;
zs0=[z-zd/2:dz:z+zd/2];
zs=zs0;

fz=zeros(size(zs));

chkflg=false;
jj=0;

while ~chkflg
    jj=jj+1;
    for ii=1:length(zs)
        if fz(ii)==0
            tbeam = beam.translateZ(zs(ii));
            sbeam = tmatrix * tbeam;
            [~,~,fz(ii),~,~,~]=ott.forcetorque(tbeam, sbeam);
        end
    end
    
    dzf=diff(fz)./diff(zs);
    
    indf=find(fz<=0,1);
    indzf=find(dzf<=0,1);
    
    if length(indf)
        if indf==1
            zs=[zs0(1:end-1)-zd*jj,zs];
            fz=[zs0(1:end-1)*0,fz];
        elseif indf==length(fz)
            zs=[zs,zs0(2:end)+zd*jj];
            fz=[fz,zs0(2:end)*0];
        else
            chkflg=true;
        end
    else
        if indzf>0
            zs=[zs,zs0(2:end)+zd*jj];
            fz=[fz,zs0(2:end)*0];
        elseif indzf<0
            zs=[zs0(1:end-1)-zd*jj,zs];
            fz=[zs0(1:end-1)*0,fz];
        else
            zs=[zs0(1:end-1)-zd*jj,zs,zs0(2:end)+zd*jj];
            fz=[zs0(1:end-1)*0,fz,zs0(2:end)*0];
        end
    end
    if jj>5
        ott.warning('external');
        error('No stable equilibrium near z!')
    end
end

% Centre and scale for polyfit
zsc = zs(max([1,indf-3]):min([length(zs),indf+3]));
zsc_mid = mean(zsc);
zsc = (zsc - zsc_mid)*beam.k_medium;

pl=polyfit(zsc - zsc_mid,fz(max([1,indf-3]):min([length(fz),indf+3])),3);

% %testing
% plot(zs,fz,'.')
% hold on
% plot(zs,polyval(pl,zs),'r')
% hold off

rts=roots(pl);
rts=rts(imag(rts)==0);

dzf=polyval([3*pl(1),2*pl(2),pl(3)],rts);

% Invert the scaling
rts = rts/beam.k_medium + zsc_mid;
dzf = dzf*beam.k_medium;

rtsi=find(dzf<0);

if ~isempty(rtsi)
    [~,indz]=min(abs(rts(rtsi)-z));
    z=rts(rtsi(indz));
    kz=dzf(rtsi(indz));
else
    ott.warning('external');
    error('No stable equilibrium near z!')
end

ott.warning('external');
