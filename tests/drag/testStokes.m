function tests = testStokes
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstruct(testCase)

  % Can't construct abstract class, construct StokesData instead
  data = rand(6, 6);
  drag = ott.drag.StokesData(data);

  testCase.verifyEqual(drag.gamma, data(1:3, 1:3), ...
    'Translation part not correct');
  testCase.verifyEqual(drag.delta, data(4:6, 4:6), ...
    'Rotation part not correct');
  testCase.verifyEqual(drag.crossterms, data(1:3, 4:6), ...
    'crossterms part not correct');

end

function testArray(testCase)

  % Test the getDefaultScalarElement method
  data(5) = ott.drag.StokesData(randn(6, 6));
  testCase.verifyEqual(data(1).gamma, eye(3), 'empty');

end

function testCastData(testCase)

  drag = ott.drag.StokesSphere(1.0);
  cast = ott.drag.StokesData(drag);

  testCase.verifyEqual(cast.forward, drag.forward, 'forward');
  testCase.verifyEqual(cast.inverse, drag.inverse, 'inverse');

end

function testFromShape(testCase)

  shape = ott.shape.Sphere(1);
  target = ott.drag.StokesSphere(shape.radius);
  
  % Single shape
  drag = ott.drag.Stokes.FromShape(shape);
  testCase.verifyEqual(drag, target, 'sphere');
  
  % Two distant shapes
  shapeArray = [shape, shape];
  shapeArray(2).position = [100;0;0];
  drag = ott.drag.Stokes.FromShape(shapeArray);
  testCase.verifyWarningFree(@() ott.drag.Stokes.FromShape(shapeArray), ...
    'far no warn');
  testCase.verifyEqual(drag(1), target, 'distant 1');
  testCase.verifyEqual(drag(2), target, 'distant 2');
  
  % Two nearby shapes
  shapeArray(2).position = [3;0;0];
  testCase.verifyWarning(@() ott.drag.Stokes.FromShape(shapeArray), ...
    'ott:drag:Stokes:no_interaction_terms', 'near warn');
  
  % Concentric spheres
  shapeArray(2).position = [0;0;0];
  shapeArray(2).radius = 0.1;
  drag = ott.drag.Stokes.FromShape(shapeArray);
  testCase.verifySize(drag, [1, 1], 'concentric');
  testCase.verifyInstanceOf(drag, 'ott.drag.EccentricSpheresNn', 'concen');
  
  % Cylinder
  shape = ott.shape.Cylinder(1, 1);
  drag = ott.drag.Stokes.FromShape(shape);
  testCase.verifyInstanceOf(drag, 'ott.drag.StokesStarShaped', 'cylinder');
  
  % Sphere-wall
  shapeArray(2) = ott.shape.Plane();
  shapeArray(1).position = [0;0;10];
  drag = ott.drag.Stokes.FromShape(shapeArray);
  shapeArray2 = ott.shape.Shape(drag);
  testCase.verifyEqual(shapeArray2(1).radius, shapeArray(1).radius, 'sphere-wall');

end

function testInverse(testCase)

  data = rand(6, 6);
  a = ott.drag.StokesData('forward', data);
  testCase.verifyEqual(inv(a), inv(data), 'Forward inverse not correct');

  b = ott.drag.StokesData('inverse', data);
  testCase.verifyEqual(b.forward, inv(data), 'Inverse inverse not correct');

end

function testMtimes(testCase)

  data = rand(6, 6);
  a = ott.drag.StokesData('inverse', data);
  testCase.verifyEqual(a*eye(6), inv(data), 'Mtimes doesn''t work');

end

function testRotation(testCase)

  import ott.utils.*;

  a = ott.drag.StokesData('forward', diag([1, 2, 3, 4, 5, 6]));
  a = a.rotateZ(pi/2);
  testCase.verifyEqual(double(a), diag([2, 1, 3, 5, 4, 6]), ...
    'AbsTol', 1.0e-6, 'Forward rotation not correct');
  testCase.verifyEqual(inv(a), diag(1./[2, 1, 3, 5, 4, 6]), ...
    'AbsTol', 1.0e-6, 'Inverse rotation not correct');
  testCase.verifyEqual(a.rotation, rotz(90)*eye(3), ...
    'Rotation not set');

  % Change the rotation back
  a.rotation = eye(3);
  testCase.verifyEqual(double(a), diag([1, 2, 3, 4, 5, 6]), ...
    'AbsTol', 1.0e-6, 'Rotation not cleared with clear');

end

function testCoverage(testCase)

  data = diag([1, 2, 3, 4, 5, 6]);
  a = ott.drag.StokesData('forward', data);
  
  testCase.verifyEqual(diag(a), diag(data), 'diag');
  testCase.verifyEqual(vecnorm(a), sum(data, 1).', 'vecnorm');
  
  testCase.verifyEqual(a.igamma, inv(diag([1,2,3])), 'igamma');
  testCase.verifyEqual(a.idelta, inv(diag([4,5,6])), 'idelta');
  testCase.verifyEqual(a.icrossterms, zeros(3, 3), 'icross');
  
end