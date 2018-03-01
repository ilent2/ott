function [r,n,rotsym,xyz_vec] = shapesurface(theta,phi,shape,parameters)
% shapesurface.m
% Generates radii and normals to the surface for a range of
% shapes for which the surface is a function of angle (theta,phi).
%
% Usage:
% [r,n,rotsym] = shapesurface(theta,phi,shape,parameters);
%
% where 
% n is the unit-length normal, stored as
%   a points-by-3 matrix with columns [ nr ntheta nphi ]
% rotsym = 3 -> the shape is a cube
% rotsym = 2 -> the shape is a sphere
% rotsym = 1 -> the shape is rotationally symmetric about the z-axis
% rotsym = 0 -> no axisymmetry
% and
% shape = 0: ellipsoid, parameters = [ a b c ]
% shape = 1: cylinder, z-axis, parameters = [ r h ]
% shape = -1: cylinder, x-axis, parameters = [ r h ]
% shape = 2: superellipsoid, paramteres = [ a b c e n ]
% shape = 3: cone-tipped cylinder, parameters = [ r h d ], d = cone height
% shape = 4: cube, parameters = [d], d = width of cube
%
% PACKAGE INFO

% For a surface defined by
% r = r(theta,phi)
% the surface area element and normal are
% dS n = r^2 sin(theta) sigma(theta,phi) dtheta dphi
% where
% sigma(theta,phi) = rhat - 1/r dr/dtheta * thetahat
%                    - 1/(r sin(theta) * dr/dphi * phihat
%
% See comments for each particular individual case

% Extract parameters and check for symmetry
switch shape
   
case  0
   % Ellipsoid
   a = parameters(1);
   b = parameters(2);
   c = parameters(3);
   if a == b
      rotsym = 1;
      if a == c
         rotsym = 2;
      end
   else
      rotsym = 0;
   end
   
case 1
    % Cylinder along the z axis
    radius = parameters(1);
    height = parameters(2);
    rotsym = 1;
   
case 2
   % Superellipsoid
   a = parameters(1);
   b = parameters(2);
   c = parameters(3);
   ew = parameters(4);
   ns = parameters(5);
   if a == b & ew == 1
      rotsym = 1;
      if a == c & ns == 1
         rotsym = 2;
      end
   else
      rotsym = 0;
   end
      
case 3
    % Cone-tipped cylinder along the z axis
    radius = parameters(1);
    height = parameters(2);
    coneheight = parameters(3);
    rotsym = 1;
    
case 4 
    % Cube
    rotsym = 0;
    d = parameters(1);
   
otherwise
   error('Unknown shape');
   
end
   
% Do we just want to ask about rotational symmetry only?
if isempty(theta)
   r = [];
   n = [];
   return;
end


theta = theta(:);
phi = phi(:);
[theta,phi] = matchsize(theta,phi);

switch shape
   
case 0
   %%%%%%%%%%%%%
   % Ellipsoid %
   %%%%%%%%%%%%%
   % r = 1/sqrt( cos(phi)^2 sin(theta)^2 / a^2 +
   %     sin(phi)^2 sin(theta)^2 / b2 + cos(theta)^2 / c^2 )
   % dr/dtheta = -r^3 ( cos(phi)^2 / a^2 + sin(phi)^2 / b^2
   %             - 1/c^2 ) sin(theta) cos(theta)
   % dr/dphi = -r^3 sin(phi) cos(phi) (1/b^2 - 1/a^2) sin(theta)^2
   % sigma = rhat + thetahat r^2 sin(theta) cos(theta) *
   %     ( cos(phi)^2/a^2 + sin(phi)^2/b^2 - 1/c^2 )
   %    + phihat r^2 sin(theta) sin(phi) cos(phi) (1/b^2 - 1//a^2)
   r = 1 ./ sqrt( (cos(phi).*sin(theta)/a).^2 + ...
      (sin(phi).*sin(theta)/b).^2 + (cos(theta)/c).^2 );
   sigma_r = ones(size(r));
   sigma_theta = r.^2 .* sin(theta) .* cos(theta) .* ...
      ( (cos(phi)/a).^2 + (sin(phi)/b).^2 - 1/c^2 );
   sigma_phi = r.^2 .* sin(theta) .* sin(phi) .* cos(phi) .* ...
      (1/b^2 - 1/a^2);
   sigma_mag = sqrt( sigma_r.^2 + sigma_theta.^2 + sigma_phi.^2 );
   n = [ sigma_r./sigma_mag sigma_theta./sigma_mag ...
         sigma_phi./sigma_mag ];
   
case 1
    %%%%%%%%%%%%
    % Cylinder %
    %%%%%%%%%%%%

    % Find angle of edges from xy plane
    edge_angle = atan(height/(2*radius));

    r = zeros(size(theta));
    n = [r r r];
    sides = find( ( theta >= pi/2 - edge_angle ) & ( theta <= pi/2 + edge_angle ) );
    top = find( theta < pi/2 - edge_angle );
    bottom = find( theta > pi/2 + edge_angle );
    ends = [ top; bottom ];

    r(sides) = radius ./ sin(theta(sides));
    r(ends) = (height/2) ./ abs(cos(theta(ends)));
    
    n(sides,1) = sin(theta(sides));
    n(sides,2) = cos(theta(sides));
    n(ends,1) = abs(cos(theta(ends)));
    n(ends,2) = - sin(theta(ends)) .* sign(cos(theta(ends)));
    
   
 
case 2
   %%%%%%%%%%%%%%%%%%
   % Superellipsoid %
   %%%%%%%%%%%%%%%%%%
   % Defined by (cartesian)
   % { (x/a)^(2/e) + (y/b)^(2/e) }^(e/n) + (z/c)^(2/n) = 1
   % where n = NS roundedness, e = EW roundedness
   % n = e = 1 reduces to ellipsoid
   % COnverting to spherical coordinates,
   % r = [ sin(theta)^(2/n) { (cos(phi)/a)^(2/e) + (sin(phi)/b)^(2/e) }^(e/n)
   %     + (cos(theta)/c)^(2/n) ]^(-n/2)
   % dr/dtheta = -r^((2+n)/n) sin(theta) cos(theta)
   %   x [ sin(theta)^((1-n)/n) { (cos(phi)/a)^(2/e) + (sin(phi)/b)^(2/e) }^(e/n)
   %             - cos(theta)^((1-n)/n) * 1/c^(2/n) ) 
   % dr/dphi = -r^((2+n)/n) sin(phi) cos(phi)
   %   x { (cos(phi)/a)^(2/e) + (sin(phi)/b)^(2/e) }^((e-n)/n)
   %   x ( sin(phi)^((1-e)/e)/b^(2/e) - cos(phi)^((1-e)/e)/a^(2/e)) sin(theta)^(2/n)
   % sigma(theta,phi) = rhat - 1/r dr/dtheta * thetahat
   %                    - 1/(r sin(theta)) * dr/dphi * phihat
   acp = abs(cos(phi));
   asp = abs(sin(phi));
   act = abs(cos(theta));
   ast = abs(sin(theta));
   cpsp = (acp/a).^(2/ew) + (asp/b).^(2/ew);
   r = ( ast.^(2/ns) .* cpsp.^(ew/ns) + (act/c).^(2/ns) ).^(-ns/2);
   sigma_r = ones(size(r));
   sigma_theta = r.^(2/ns) .* ast .* cos(theta) .* ...
      ( ast.^((1-ns)/ns) .* cpsp.^(ew/ns) - act.^((1-ns)/ns)/c^(2/ns) );
   sigma_phi = r.^(2/ns) .* ast.^((2-ns)/ns) .* sin(phi) .* cos(phi) .* ...
      cpsp.^((ew-ns)/ns) .* ( asp.^((1-ew)/ew)/b^(2/ew) - acp.^((1-ew)/ew)/a^(2/ew) );
   sigma_mag = sqrt( sigma_r.^2 + sigma_theta.^2 + sigma_phi.^2 );
   n = [ sigma_r./sigma_mag sigma_theta./sigma_mag ...
         sigma_phi./sigma_mag ];
   
case 3
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Cone-tipped Cylinder %
    %%%%%%%%%%%%%%%%%%%%%%%%

    % Find angle of edges from xy plane
    edge_angle = atan(height/(2*radius));
    cone_angle = atan(coneheight/radius);
    
    r = zeros(size(theta));
    n = [r r r];
    sides = find( ( theta >= pi/2 - edge_angle ) & ( theta <= pi/2 + edge_angle ) );
    top = find( theta < pi/2 - edge_angle );
    bottom = find( theta > pi/2 + edge_angle );
    ends = [ top; bottom ];

    r(sides) = radius ./ sin(theta(sides));
    r(top) = (height/2 + coneheight) * cos(cone_angle) ./ cos( theta(top) - cone_angle );
    r(bottom) = (height/2 + coneheight) * cos(cone_angle) ./ ...
        abs(cos( theta(bottom) - pi + cone_angle ));
    
    n(sides,1) = sin(theta(sides));
    n(sides,2) = cos(theta(sides));
    n(top,1) = cos( theta(top) - cone_angle );
    n(top,2) = sin( cone_angle - theta(top) );
    n(bottom,1) = abs(cos( theta(bottom) - pi + cone_angle ));
    n(bottom,2) = sin( pi - cone_angle - theta(bottom) );
 
%    [theta(top) r(top) n(top,1) n(top,2) ]
%    [theta(bottom) r(bottom) n(bottom,1) n(bottom,2) ]
    
case 4
    %%%%%%%%
    % CUBE %
    %%%%%%%%
    
    % First lets pretend our cube is a sphere - then we will collapse the
    % sphere down to a cube
    r_prov = sqrt(3*(d/2)^2);
% We want this to also work for vector inputs of theta and phi

theta_vec = theta(:);
phi_vec = phi(:);
[theta_vec,phi_vec] = matchsize(theta_vec,phi_vec);
r_cube_vec = zeros(length(theta_vec),1);
n_vec = zeros(length(theta_vec),3);
xyz_vec = zeros(length(theta_vec),3);

for j = 1:length(theta_vec)
    theta = theta_vec(j);
    phi = phi_vec(j);

% The cartesian coordiantes will be
[x y z] = rtp2xyz(r_prov,theta,phi);

% Now to collapse the sphere down to a cube

if abs(x) == d/2 & abs(y) == d/2 & abs(z) == d/2
a=1;
x2 = a*x;
y2 = a*y;
z2 = a*z;
n = [1 0 0];
elseif abs(x) > d/2 & abs(x) > abs(y) & abs(x) > abs(z)
a = (d/2)/abs(x);
x2 = a*x;
y2 = a*y;
z2 = a*z;
if x > d/2
    n = [sin(theta)*cos(phi) cos(theta)*cos(phi) -sin(phi)];
    norm = sqrt(n(1)^2 + n(2)^2 + n(3)^2);
    n = n./norm;
elseif x < -d/2
    n = [-sin(theta)*cos(phi) -cos(theta)*cos(phi) sin(phi)];
    norm = sqrt(n(1)^2 + n(2)^2 + n(3)^2);
    n = n./norm;
end
elseif abs(y) > d/2 & abs(y) > abs(x) & abs(y) > abs(z)
a = (d/2)/abs(y);
x2 = a*x;
y2 = a*y;
z2 = a*z;
if y > d/2
    n = [sin(theta)*sin(phi) cos(theta)*sin(phi) cos(phi)];
    norm = sqrt(n(1)^2 + n(2)^2 + n(3)^2);
    n = n./norm;
elseif y < -d/2
    n = [-sin(theta)*sin(phi) -cos(theta)*sin(phi) -cos(phi)];
    norm = sqrt(n(1)^2 + n(2)^2 + n(3)^2);
    n = n./norm;
end
elseif abs(z) > d/2 & abs(z) > abs(x) & abs(z) > abs(y)
a = (d/2)/abs(z);
x2 = a*x;
y2 = a*y;
z2 = a*z;
if z > d/2
    n = [cos(theta) -sin(theta) 0];
    norm = sqrt(n(1)^2 + n(2)^2 + n(3)^2);
    n = n./norm;
elseif z < -d/2
    n = [-cos(theta) sin(theta) 0];
    norm = sqrt(n(1)^2 + n(2)^2 + n(3)^2);
    n = n./norm;
end

end

% The new polar coordinates are
[r_cube theta_cube phi_cube] = xyz2rtp(x2, y2, z2);

r_cube_vec(j) = r_cube;
n_vec(j,:) = n;

end

r = r_cube_vec;
n = n_vec;

otherwise
   error('Unknown shape');
end


return
