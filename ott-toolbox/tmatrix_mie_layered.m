function [t,t2,a,b] = tmatrix_mie_layered(nmax,k_medium,k_particle,radius)
% tmatrix_mie.m : mie scattering and internal coefficients for a uniform or
%                 layered sphere arranged as a sparse t-matrix.
%
% usage:
% [t,t2] = tmatrix_mie(nmax,k_medium,k_particle,radius)
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
% This file is part of the package Optical tweezers toolbox 1.3
% Copyright 2006-2013 The University of Queensland.
% See README.txt or README.m for license and details.
%
% http://www.physics.uq.edu.au/people/nieminen/software.html

n_layer=[k_particle,k_medium]/2/pi; %medium is on the outside...
radius=radius; %radius for the medium layer is the "same" as final layer.

n=[1:nmax]; %n for nmax

[ha_0,hb_0]=layerrecursion(n,n_layer,radius,1,1/n_layer(1),1);

m_0 = n_layer(end-1);
m_1 = n_layer(end);

r1 = 2*pi*m_1 * radius(end);

j1 = (sbesselj(n,r1)).';

h1 = (sbesselh1(n,r1)).';

j1d = (sbesselj(n-1,r1) - n.*sbesselj(n,r1)/r1).';

h1d = (sbesselh1(n-1,r1) - n.*sbesselh1(n,r1)/r1).';

d1_1=j1d./j1;
d3_1=h1d./h1;
r_1=j1./h1;

% %TM/TE coeffs...
al_1=r_1.*(m_1*ha_0-m_0*d1_1)./(m_1*ha_0-m_0*d3_1);
bl_1=r_1.*(m_0*hb_0-m_1*d1_1)./(m_0*hb_0-m_1*d3_1);
% aL_1=((ha_0/m_0+n.'/r1).*j1-sbesselj(n-1,r1).')./((ha_0/m_0+n.'/r1).*h1-sbesselh1(n-1,r1).');
% bL_1=((m_0*hb_0+n.'/r1).*j1-sbesselj(n-1,r1).')./((m_0*hb_0+n.'/r1).*h1-sbesselh1(n-1,r1).');

%swap the modes to TE/TM... and negatize.
b = -al_1;
a = -bl_1;

%t-matrix indices...
indexing=combined_index(1:nmax^2+2*nmax)';

t=sparse([1:2*(nmax^2+2*nmax)],[1:2*(nmax^2+2*nmax)],[a(indexing);b(indexing)]);

if nargout>1
    
    warning('This release does not calculate internal coefficients.') %because I have to use different ratio functions to find c/d
    t2=sparse([1:2*(nmax^2+2*nmax)],[1:2*(nmax^2+2*nmax)],[a(indexing);b(indexing)]);
    
end

end

function [ha_0,hb_0]=layerrecursion(n,n_layer,radius,recursionnumber,ha_n,hb_n)
%Internal function which should never be bare. This function has a lot of
%overhead but takes very little time to complete so I'm happy with it as
%is. I made the recursion start from 1 because it should...

m_0 = n_layer(recursionnumber);
r0 = 2*pi*m_0 * radius(recursionnumber);
j0 = (sbesselj(n,r0)).';
h0 = (sbesselh1(n,r0)).';
j0d = (sbesselj(n-1,r0) - n.*sbesselj(n,r0)/r0).';
h0d = (sbesselh1(n-1,r0) - n.*sbesselh1(n,r0)/r0).';

d1_0=j0d./j0;
d3_0=h0d./h0;
r_0=j0./h0;

if recursionnumber>1
    m_n = n_layer(recursionnumber-1);
    rn = 2*pi*m_0 * radius(recursionnumber-1);
    jn = (sbesselj(n,rn)).';
    hn = (sbesselh1(n,rn)).';
    jnd = (sbesselj(n-1,rn) - n.*sbesselj(n,rn)/rn).';
    hnd = (sbesselh1(n-1,rn) - n.*sbesselh1(n,rn)/rn).';
    
    d1_n=jnd./jn;
    d3_n=hnd./hn;
    r_n=jn./hn;
else
    m_n=1;
    rn=0;
    jn=0;
    hn=0;
    jnd=0;
    hnd=0;
    
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

if recursionnumber<length(radius) %doesn't calculate last layer with this recursion...
    recursionnumber=recursionnumber+1;
    [ha_0,hb_0]=layerrecursion(n,n_layer,radius,recursionnumber,ha_0,hb_0);
end

end