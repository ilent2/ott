function tests = testStokesLambNn()
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstruct(testCase)
  
  shape = ott.shapes.Sphere();
  viscosity = 2.0;
  drag = ott.drag.StokesLambNn(shape, viscosity);
  
  testCase.verifyEqual(drag.shape, shape, 'shape');
  testCase.verifyEqual(drag.viscosity, viscosity, 'eta');
end

function testSphere(testCase)

  radius = 1.0;
  shape = ott.shapes.Sphere(radius);
  viscosity = 2.0;
  
  drag = ott.drag.StokesLambNn(shape, viscosity);
  target = ott.drag.StokesSphere(shape.radius, viscosity);
  
  testCase.verifyEqual(drag.forward, target.forward, ...
    'AbsTol', 2.6e-2, 'RelTol', 1.0e-2, 'forward');

end

function testCylinder(testCase)

  % Want cylinder with L/(2*R) = 1/0.5 (extent of cylinder method)
  shape = ott.shapes.Cylinder(1, 4);
  viscosity = 2.0;
  
  drag = ott.drag.StokesLambNn(shape, viscosity);
  target = ott.drag.StokesCylinder(shape.radius, shape.height, viscosity);
  
  testCase.verifyEqual(drag.forward, target.forward, ...
    'RelTol', 0.2, 'AbsTol', 3, 'forward');

end


