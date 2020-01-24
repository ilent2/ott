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

  radius = 1.0;
  outerRadius = 1000.0;
  viscosity = 1.0;
  separation = 900.0;

  a = ott.drag.EccentricSpheresNn(radius, outerRadius, separation, viscosity);
  b = ott.drag.FaxenSphere(radius, separation, viscosity);
  
  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  import matlab.unittest.constraints.RelativeTolerance;
  testCase.verifyThat(a.forward, IsEqualTo(b.forward, ...
      'Within', AbsoluteTolerance(3e-2) | RelativeTolerance(1e-2)));
  testCase.verifyThat(a.inverse, IsEqualTo(b.inverse, ...
      'Within', AbsoluteTolerance(3e-2) | RelativeTolerance(1e-2)));

end

function testFarLimit(testCase)

  radius = 1.0;
  outerRadius = 100000.0;
  viscosity = 1.0;
  separation = 10000.0;

  a = ott.drag.EccentricSpheresNn(radius, outerRadius, separation, viscosity);
  b = ott.drag.StokesSphere(radius, viscosity);

  testCase.verifyEqual(a.forward, b.forward, ...
      'AbsTol', 3.0e-2, ...
      'Far-limit doesn''t match forward');
  testCase.verifyEqual(a.inverse, b.inverse, ...
      'AbsTol', 3.0e-2, ...
      'Far-limit doesn''t match inverse');

end

