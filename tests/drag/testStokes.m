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

function testCastData(testCase)

  drag = ott.drag.StokesSphere(1.0);
  cast = ott.drag.StokesData(drag);

  testCase.verifyEqual(cast.forward, drag.forward, 'forward');
  testCase.verifyEqual(cast.inverse, drag.inverse, 'inverse');

end

function testFromShape(testCase)

  shape = ott.shape.Sphere();
  drag = ott.drag.Stokes.FromShape(shape);
  target = ott.drag.StokesSphere(shape.radius);

  testCase.verifyEqual(drag, target, 'sphere');

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
  a.rotation = eye(3);;
  testCase.verifyEqual(double(a), diag([1, 2, 3, 4, 5, 6]), ...
    'AbsTol', 1.0e-6, 'Rotation not cleared with clear');

end
