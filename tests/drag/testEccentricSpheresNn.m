function tests = testEccentricSpheresNn
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstruction(testCase)

  innerRadius = 1.0;
  outerRadius = 1.2;
  separation = 0.1;
  viscosity = 1.6;
  a = ott.drag.EccentricSpheresNn(innerRadius, outerRadius, ...
    separation, viscosity);

  testCase.verifyEqual(a.innerRadius, innerRadius, 'Inner radius not set');
  testCase.verifyEqual(a.outerRadius, outerRadius, 'Outer radius not set');
  testCase.verifyEqual(a.separation, separation, 'Separation not set');
  testCase.verifyEqual(a.viscosity, viscosity, 'Viscosity not set');

end

function testFaxenLimit(testCase)
  % Compare against faxen's correction (valid near plane surface,
  % at least 1 radius spacing between wall and surface).
  %
  % This is not the best comparison, since moving far from the wall
  % reduces Faxen's correction to a Stokes sphere, but too close to the
  % wall and Faxen's correction breaks.  While for the eccentric sphere
  % case, if we are too far from the wall, the surface no longer looks
  % like a infinite plane (i.e. no longer Faxen's).

  radius = 1.0;
  outerRadius = 1e5;
  viscosity = 3.0;
  separation = 2.21;  % Faxen's should work to about 1 radius (i.e. sep=2)

  a = ott.drag.EccentricSpheresNn(radius, outerRadius, separation-radius, viscosity);
  b = ott.drag.FaxenSphere(radius, separation, viscosity);
  
  testCase.verifyEqual(a.forward, b.forward, ...
      'RelTol', 0.12);
  testCase.verifyEqual(a.inverse, b.inverse, ...
      'RelTol', 0.12);

end

function testFarLimit(testCase)

  radius = 1.0;
  outerRadius = 100000.0;
  viscosity = 1.0;
  separation = 10000.0;

  a = ott.drag.EccentricSpheresNn(radius, outerRadius, separation, viscosity);
  b = ott.drag.StokesSphere(radius, viscosity);

  testCase.verifyEqual(a.forward, b.forward, ...
      'RelTol', 2.0e-3, 'AbsTol', viscosity*3e-2, ...
      'Far-limit doesn''t match forward');
  testCase.verifyEqual(a.inverse, b.inverse, ...
      'RelTol', 2.0e-3, 'AbsTol', 1e-4./viscosity, ...
      'Far-limit doesn''t match inverse');

end

function testViscosity(testCase)

  radius = 1.0;
  outerRadius = 100.0;
  separation = 10.0;
  
  a = ott.drag.EccentricSpheresNn(radius, outerRadius, separation, 1.0);
  b = ott.drag.EccentricSpheresNn(radius, outerRadius, separation, 2.0);
  
  testCase.verifyEqual(a.forward, 0.5*b.forward, ...
    'RelTol', 1.0e-16, 'Viscosity has incorrect effect');
end


