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

import ott.*
ott_warning('internal');

if nargin < 4
    z = 0;
end

% Ensure T-matrix and beam are the same size
Nmax = max(tmatrix.Nmax, beam.Nmax);
tmatrix.Nmax = Nmax;
beam.Nmax = Nmax;

equiv_ka = nmax2ka(Nmax);

% Normalise the beam power
power=beam.power();
beam = beam / power;

%start with three points each side over a 1/8 wavelength:
zd=.25;
dz=zd/5;
zs0=[z-zd/2:dz:z+zd/2];
zs=zs0;

fz=zeros(size(zs));

chkflg=false;
jj=0;

p=zeros(length(a),1);
q=p;

while ~chkflg
    jj=jj+1;
    for ii=1:length(zs)
        if fz(ii)==0
            tbeam = beam.translateZ(zs(ii));
            sbeam = tmatrix * tbeam;
            [~,~,fz(ii),~,~,~]=forcetorque(tbeam, sbeam);
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
        ott_warning('external');
        error('No stable equilibrium near z!')
    end
end
pl=polyfit(zs(max([1,indf-3]):min([length(zs),indf+3])),fz(max([1,indf-3]):min([length(fz),indf+3])),3);

% %testing
% plot(zs,fz,'.')
% hold on
% plot(zs,polyval(pl,zs),'r')
% hold off

rts=roots(pl);
rts=rts(imag(rts)==0);

dzf=polyval([3*pl(1),2*pl(2),pl(3)],rts);

rtsi=find(dzf<0);

if length(rtsi)
    [d,indz]=min(abs(rts(rtsi)-z));
    z=rts(rtsi(indz));
    kz=dzf(rtsi(indz));
else
    ott_warning('external');
    error('No stable equilibrium near z!')
end

ott_warning('external');

return
