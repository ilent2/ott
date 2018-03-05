function [alpha,beta,gamma]=axisangletoeuler(nxyz,w);
% axisangletoeuler.m --- converts the axis angle into alpha, beta, gamma of 
% the wigner-D functions.

if nargin<2
    w=norm(nxyz);
end

if norm(nxyz)==0|w==0
	alpha=0;
	beta=0;
	gamma=0;
	if nargout==1
		alpha=[alpha,beta,gamma];
	end
	return
end

nxyz=nxyz/norm(nxyz);

%this will work stably
phi=atan2(nxyz(2),nxyz(1));
theta=atan2(sqrt(sum(nxyz(1:2).^2)),nxyz(3));

%definitions in Quantum theory of Angular momentum pp. 26.
aplusgon2=atan2(cos(theta)*sin(w/2),cos(w/2));
gamma=aplusgon2-phi+pi/2;
alpha=2*phi-pi+gamma;

beta=2*atan2(sin(aplusgon2)*sin(theta),cos(theta));

if nargout==1
    alpha=[alpha,beta,gamma];
end
