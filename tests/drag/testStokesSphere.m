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

  testCase.verifyEqual(a(1, 1), a(2, 2), 'A(1, 1) should equal A(2, 2)');
  testCase.verifyEqual(a(1, 1), a(3, 3), 'A(1, 1) should equal A(3, 3)');
  testCase.verifyEqual(a(4, 4), a(5, 5), 'A(4, 4) should equal A(5, 5)');
  testCase.verifyEqual(a(4, 4), a(6, 6), 'A(4, 4) should equal A(6, 6)');
end

