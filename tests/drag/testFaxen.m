function tests = testFaxen
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstruction(testCase)

  radius = 1.0;
  viscosity = 1.4;
  separation = 1.2;
  a = ott.drag.FaxenSphere(radius, separation, viscosity);

  testCase.verifyEqual(a.separation, separation, 'Separation not set');
  testCase.verifyEqual(a.radius, radius, 'Separation not set');
  testCase.verifyEqual(a.viscosity, viscosity, 'Separation not set');

  testCase.verifyTrue(isdiag(a.forward), 'Drag tensor should be diagonal');

  testCase.verifyNotEqual(a.forward(1, 1), a.forward(3, 3), ...
    'A(1,1) should not equal A(3,3)');
  testCase.verifyEqual(a.forward(1, 1), a.forward(2, 2), ...
    'A(1,1) should equal A(2,2)');
  testCase.verifyEqual(a.forward(4, 4), a.forward(5, 5), ...
    'A(4,4) should equal A(5,5)');
  testCase.verifyNotEqual(a.forward(5, 5), a.forward(6, 6), ...
    'A(5,5) should equal A(6,6)');

end

function testFarLimit(testCase)

  radius = 1.0;
  viscosity = 1.0;
  separation = 10000.0;

  a = ott.drag.FaxenSphere(radius, separation, viscosity);
  b = ott.drag.StokesSphere(radius, viscosity);

  testCase.verifyEqual(a.forward, b.forward, ...
      'RelTol', 1.0e-3, ...
      'Far-limit doesn''t match forward');
  testCase.verifyEqual(a.inverse, b.inverse, ...
      'RelTol', 1.0e-3, ...
      'Far-limit doesn''t match inverse');

end

% function testPlot(testCase)
%   % Attempt to reproduce figure 3 from J. Leach, et al.
%   % https://doi.org/10.1103/PhysRevE.79.026301
% 
%   radius = 2.0e-6;    % 2 um
%   separation = logspace(0, 1, 100)*radius;
%   viscosity = 0.001;
% 
%   values = zeros(4, length(separation));
%   
%   sphere = ott.drag.StokesSphere(radius, viscosity);
%   sphere = [sphere(1, 1); sphere(3, 3); sphere(4, 4); sphere(6, 6)];
% 
%   for ii = 1:length(separation)
%     drag = ott.drag.FaxenSphere(radius, separation(ii), viscosity);
%     values(:, ii) = [drag(1, 1); drag(3, 3); drag(4, 4); drag(6, 6)];
%   end
% 
%   figure();
%   loglog(separation./radius, 1.0 - sphere./values);
%   axis([1e0, 1e1, 1e-2, 1e0]);
% 
% end

