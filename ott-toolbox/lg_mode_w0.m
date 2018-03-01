function w0 = lg_mode_w0( mode, angle )
% lg_mode_w0.m - gives paraxial waist parameter w0 to match a given 1/e^2 angle
%
% Usage:
% w0 = lg_mode_w0( mode, angle );
% where:
% w0 is the waist parameter in units of wavelength
% mode = [ p l ] [BUT CURRENTLY p MUST BE ZERO]
%  or
% mode = l (and p = 0 is assumed)
% angle = 1/e^2 angle in degrees
% 
% Accuracy not guaranteed for very large l (ie > 200)
%
% PACKAGE INFO

% Precalculated results for small l
psi1 = [ 2.25262074788770;
   3.14619322065846;
   3.95382908471947;
   4.71535334785585;
   5.44688108968330;
   6.15686752035254;
   6.85040872419631;
   7.53085901035229;
   8.20056601716351;
   8.86124914810352;
   9.51421256077447;
   10.16047339493525;
   10.80084315888962;
   11.43598163902789;
   12.06643388775774 ];

% Approximation formula for large l
psi2 = [ -4.277910129402566e-017;
   4.133454181733631e-014;
   -1.705821602342175e-011;
   3.938104134621960e-009;
   -5.612672442622125e-007;
   5.185761533022050e-005;
   -0.00328272527981;
   0.69510313966730;
   2.23523934080996 ];

psi3 = [ -4.859972437674526e-006;
   0.52767748587717;
   9.47541396001228 ];

if length(mode) > 1
   l = abs(mode(2));
else
   l = mode;
end

if l == 0
   % Gaussian
   psi = 1;
elseif l < 16
   psi = psi1(l);
elseif l < 201
   psi = polyval(psi2,l);
else
   psi = polyval(psi3,l);
end

% OK, now we have psi!
%
% psi = k^2 w_0^2 tan(theta)^2 /4
%
% We want w_0 in wavelengths, so k = 2*pi

w0 = sqrt(psi)/(pi*tan(angle*pi/180));

return
