function tests = chebyshev
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testSimple2d(testCase)

  radius = 1.0;
  sscale = 0.1;
  order = 10;
  shape = ott.shapes.Chebyshev2d(radius, sscale, order);

  testCase.verifyEqual(shape.radius, radius, ...
    'shape radius not set');
  testCase.verifyEqual(shape.sscale, sscale, ...
    'shape perterbation not set');
  testCase.verifyEqual(shape.order, order, ...
    'shape order not set');
end

function testSimple3d(testCase)

  radius = 1.0;
  sscale = 0.1;
  order = 10;
  shape = ott.shapes.Chebyshev3d(radius, sscale, order);

  testCase.verifyEqual(shape.radius, radius, ...
    'shape radius not set');
  testCase.verifyEqual(shape.sscale, sscale, ...
    'shape perterbation not set');
  testCase.verifyEqual(shape.order, order, ...
    'shape order not set');
end

% function testSurfaceNorms(testCase)
% 
%   radius = 1.0;
%   sscale = 0.1;
%   order = 10;
% %   shape = ott.shapes.Chebyshev2d(radius, sscale, order);
%   shape = ott.shapes.Chebyshev3d(radius, sscale, order);
%   
%   figure();
%   shape.surf();
%   hold on;
%   phi = pi/2.0 * ones(1, 20);
%   theta = linspace(0, pi, 20);
%   n = shape.normals(theta, phi);
%   r = shape.radii(theta, phi);
%   [nxyz, xyz] = ott.utils.rtpv2xyzv(n, [r(:), theta(:), phi(:)]);
%   s = 0.1;
%   plot3([xyz(:, 1), xyz(:, 1)+s*nxyz(:, 1)].', ...
%     [xyz(:, 2), xyz(:, 2)+s*nxyz(:, 2)].', ...
%     [xyz(:, 3), xyz(:, 3)+s*nxyz(:, 3)].');
% end
