function tests = testSphere
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
  
  warning('on', 'ott:drag:StokesCylinder:outside_range');
end

function testConstruction(testCase)

  radius = 1.0;
  height = 4.0;
  viscosity = 1.4;
  a = ott.drag.StokesCylinder(radius, height, viscosity);

  testCase.verifyEqual(a.radius, radius, 'Radius not set');
  testCase.verifyEqual(a.height, height, 'height not set');
  testCase.verifyEqual(a.viscosity, viscosity, 'Viscosity not set');

  testCase.verifyTrue(isdiag(a.forward), 'Drag tensor should be diagonal');
end

function testAspectWarning(testCase)

  radius = 1.0;
  height = 0.5;
  drag = ott.drag.StokesCylinder(radius, height);
  
  testCase.verifyWarning(@() drag.forward, ...
      'ott:drag:StokesCylinder:outside_range');

end


