function [t,t2,a,b] = tmatrix_mie_layered(nmax,k_medium,k_particle,radius)
% TMATRIX_MIE_LAYERED mie scattering and internal coefficients for a
% uniform or layered sphere arranged as a sparse t-matrix.
%
% [t,t2] = TMATRIX_MIE(nmax,k_medium,k_particle,radius)
% k_medium is the wavenumber in the surrounding medium.
% k_particle is the wavenumber in the layers starting at the core going to
%   the outermost layer.
% radius is the distance from the center of the core for each boundary.
% t is the t-matrix of scattered modes. t2 is the t-matrix of internal
% modes.
%
% Note: Will not compute for particles of 500 layers or more. Not all of
% the recursions in the reference are implemented because the author was
% lazy. Nmax of 100 is numerically stable.
%
% Reference:
% "Improved recursive algorithm for light scattering by a multilayered
% sphere", Wen Yang, Applied Optics 42(9), 2003
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

ott_warning('internal');

k_layer=[k_particle,k_medium]; %medium is on the outside...
n_layer=k_layer/2/pi;
radius=[radius,radius(end)]; %radius for the medium layer is the "same" as final layer.

n=[1:nmax]; %n for nmax

lastElement=length(k_layer);

%generate all special functions first:
if length(k_particle)>1
    [jN,jNd] = sbesselj(n,[k_layer.*radius,k_layer(2:end).*radius(1:end-1)]);
    [hN,hNd] = sbesselh1(n,[k_layer.*radius,k_layer(2:end).*radius(1:end-1)]);
else
    [jN,jNd] = sbesselj(n,k_layer(:).*radius(:));
    [hN,hNd] = sbesselh1(n,k_layer(:).*radius(:)); 
end
jN=jN.';
hN=hN.';
jNd=jNd.';
hNd=hNd.';

d1_N=jNd./jN;
d3_N=hNd./hN;
r_N=jN./hN;

ha_0=1/n_layer(1);
hb_0=1;

for ii=1:length(k_particle)
    ha_n=ha_0;
    hb_n=hb_0;
    
    m_0 = n_layer(ii);
    
    d1_0=d1_N(:,ii);
    d3_0=d3_N(:,ii);
    r_0=r_N(:,ii);
    
    if ii>1
        m_n = n_layer(ii-1);
        
        d1_n=d1_N(:,lastElement+ii-1);
        d3_n=d3_N(:,lastElement+ii-1);
        r_n=r_N(:,lastElement+ii-1);
    else
        m_n=1;
        
        d1_n=0;
        d3_n=0;
        r_n=0;
    end
    
    g0=m_0*ha_n-m_n*d1_n;
    g1=m_0*ha_n-m_n*d3_n;
    g0b=m_n*hb_n-m_0*d1_n;
    g1b=m_n*hb_n-m_0*d3_n;
    q_0=r_n./r_0;
    ha_0=(g1.*d1_0-q_0.*g0.*d3_0)./(g1-q_0.*g0);
    hb_0=(g1b.*d1_0-q_0.*g0b.*d3_0)./(g1b-q_0.*g0b);

end
m_0 = n_layer(lastElement-1);
m_1 = n_layer(lastElement);

m=m_0/m_1;

d1_1=d1_N(:,lastElement);
d3_1=d3_N(:,lastElement);
r_1=r_N(:,lastElement);

% %TM/TE coeffs...
al_1=r_1.*(ha_0-m*d1_1)./(ha_0-m*d3_1);
bl_1=r_1.*(m*hb_0-d1_1)./(m*hb_0-d3_1);

%swap the modes to TE/TM... and negatize.
b = -al_1;
a = -bl_1;

%t-matrix indices...
indexing=combined_index(1:nmax^2+2*nmax)';

t=sparse([1:2*(nmax^2+2*nmax)],[1:2*(nmax^2+2*nmax)],[a(indexing);b(indexing)]);

ott_warning('external');

if nargout>1
    
    r_0=(jN(:,lastElement)./jN(:,lastElement-1));
    
    d = r_0.*(d3_1 - d1_1 )  ./ (ha_0-m*d3_1);
    c = r_0.*(d3_1 - d1_1 )  ./ (m*hb_0 - d3_1);
    
    ott_warning('ott:tmatrix_mie_layered:internalcoefficientwarning', ...
        ['The internal coefficients are for the outermost layer only...' ...
         ' the real ones are only defined for each layer.']);
    t2=sparse([1:2*(nmax^2+2*nmax)],[1:2*(nmax^2+2*nmax)],[c(indexing);d(indexing)]);
    
end
