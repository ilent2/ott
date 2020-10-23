function stRtfunc = sphMakeGeometry(nNbTheta, a, c, theta)
%% sphMakeGeometry
% Evaluates the functions defining the geometry r(theta) for spheroids
% 
% sphMakeGeometry(nNbTheta, a, c, theta) calculates the geometry for
% spheroids and the nodes and weights of the Gaussian quadrature.
% By reflection symmetry, only 0<theta<pi/2 are used; quadrature weights are doubled
% to compensate. If theta is specified, the quadrature properties are not computed.
%
% Input:
%   nNbTheta:   the number of angles to use (in half of the shape)
%	a:          the semi-axis length for axes along x, y
%	c:          the semi-axis length along z (axis of rotation)
%   theta:      (optional) theta values to calculate the geometry for; else uses
%                   lgwt to generate quadrature
%
% Output:
%		stRtfunc	A structure containing fields
%           - theta: [T x 1] angles where things are computed (if not
%                    specified in the parameters, these will be the quadrature nodes
%		    - r: [T x 1] r(theta) defining the geometry
%           - drdt: [T x 1] Derivative of r(theta)
%           - wTheta: [T x 1] Weights from quadrature
%           - nNbTheta: [1 x 1] Number of theta
%           - a [1 x 1] Input parameter
%           - c [1 x 1] Input parameter
%			- h: [1 x 1] The aspect ratio, max(r)/min(r)
%           - r0: [1 x 1] Equivalent-volume-sphere radius
%           - sInt: String specifying 'Gauss' for quadrature or 'Pts'
%             if theta are specified
%
% Dependency: 
% auxPrepareIntegrals

import ott.tmatrix.smarties.*;

sInt='Gauss';

if (nargin == 4) % points only
    if (size(theta,2)~=1)
        disp 'makeGeometry: theta should be a column vector... transposing';
        theta=transpose(theta);
    end
    stRtfunc.wTheta=0*theta;
    stRtfunc.theta=theta;
    stRtfunc.nNbTheta=length(theta);
    sInt='Pts';
else
    stRtfunc=auxPrepareIntegrals(2*nNbTheta,sInt);
    stRtfunc.theta=stRtfunc.theta(1:nNbTheta); % [T x 1]
	theta = stRtfunc.theta;
	stRtfunc.wTheta = stRtfunc.wTheta(1:nNbTheta)*2;
	stRtfunc.nNbTheta = nNbTheta;

end

stRtfunc.a=a;
stRtfunc.c=c;
stRtfunc.sInt=sInt;


% Defines geometry
sint=sin(theta);
cost=cos(theta);

% r for spheroid with different semi major and semi minor axes at the origin
stRtfunc.r = a*c./realsqrt(c^2*sint.^2+a^2*cost.^2); %[T x 1]
stRtfunc.drdt = (a^2-c^2)/(a*c)^2*sint.*cost.*(stRtfunc.r).^3;   %[T x 1]

stRtfunc.h = max([a,c])/min([a,c]); % aspect ratio, >= 1

stRtfunc.r0 = (a^2*c)^(1/3); % the radius of volume-equivalent sphere

end
