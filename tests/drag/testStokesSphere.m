function tests = testSphere
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstruction(testCase)

  radius = 1.0;
  viscosity = 1.4;
  a = ott.drag.StokesSphere(radius, viscosity);

  testCase.verifyEqual(a.radius, radius, 'Radius not set');
  testCase.verifyEqual(a.viscosity, viscosity, 'Viscosity not set');

  testCase.verifyTrue(isdiag(a.forward), 'Drag tensor should be diagonal');

  testCase.verifyEqual(a.forward(1, 1), a.forward(2, 2), ...
    'A(1, 1) should equal A(2, 2)');
  testCase.verifyEqual(a.forward(1, 1), a.forward(3, 3), ...
    'A(1, 1) should equal A(3, 3)');
  testCase.verifyEqual(a.forward(4, 4), a.forward(5, 5), ...
    'A(4, 4) should equal A(5, 5)');
  testCase.verifyEqual(a.forward(4, 4), a.forward(6, 6), ...
    'A(4, 4) should equal A(6, 6)');
end

function testSimple(testCase)
  error('not yet implemented');
end

