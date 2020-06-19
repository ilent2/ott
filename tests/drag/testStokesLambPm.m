function tests = testStokesLambPm()
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstruct(testCase)
  
  shape = ott.shapes.Sphere();
  viscosity = 2.0;
  drag = ott.drag.StokesLambPm(shape, viscosity);
  
  testCase.verifyEqual(drag.shape, shape, 'shape');
  testCase.verifyEqual(drag.viscosity, viscosity, 'eta');
end

function testSphere(testCase)
  % Compare to StokesSphere
  
  shape = ott.shapes.Sphere();
  viscosity = 2.0;
  
  drag = ott.drag.StokesLambPm(shape, viscosity);
  target = ott.drag.StokesSphere(shape.radius, viscosity);
  
  testCase.verifyEqual(drag.forward, target.forward, ...
    'AbsTol', 1.0e-13, 'forward');
end

function testCylinder(testCase)

  % Want cylinder with L/(2*R) = 1/0.5 (extent of cylinder method)
  shape = ott.shapes.Cylinder(1, 4);
  viscosity = 2.0;
  
  drag = ott.drag.StokesLambPm(shape, viscosity);
  drag.seriesOrder = 10;   % Need a larger order
  target = ott.drag.StokesCylinder(shape.radius, shape.height, viscosity);
  
  testCase.verifyEqual(drag.forward, target.forward, ...
    'RelTol', 0.2, 'AbsTol', 1.0e-13, 'forward');

end
