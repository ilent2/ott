function tests = testFaxen
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstruct(testCase)

  a = ott.drag.Stokes();
  testCase.verifyTrue(isempty(a.forward), 'Forward not empty');
  testCase.verifyTrue(isempty(a.inverse), 'Inverse not empty');

  data = rand(6, 6);
  b = ott.drag.Stokes('forward', data);
  testCase.verifyEqual(b.forward, data, 'Forward not set correctly');
  testCase.verifyEqual(b.translation, data(1:3, 1:3), ...
    'Translation part not correct');
  testCase.verifyEqual(b.rotation, data(4:6, 4:6), ...
    'Rotation part not correct');

end

function testInverse(testCase)

  data = rand(6, 6);
  a = ott.drag.Stokes('forward', data, 'finalize', true);
  testCase.verifyEqual(inv(a), inv(data), 'Forward inverse not correct');

  b = ott.drag.Stokes('inverse', data, 'finalize', true);
  testCase.verifyEqual(b.forward, inv(data), 'Inverse inverse not correct');

end

function testMtimes(testCase)

  data = rand(6, 6);
  a = ott.drag.Stokes('inverse', data, 'finalize', true);
  testCase.verifyEqual(a*eye(6), inv(data), 'Mtimes doesn''t work');

end

function testRotation(testCase)

  import ott.utils.*;

  a = ott.drag.Stokes('forward', diag([1, 2, 3, 4, 5, 6]));
  a = a.rotateZ(pi/2);
  testCase.verifyEqual(double(a), diag([2, 1, 3, 5, 4, 6]), ...
    'AbsTol', 1.0e-6, 'Forward rotation not correct');
  testCase.verifyEqual(inv(a), diag(1./[2, 1, 3, 5, 4, 6]), ...
    'AbsTol', 1.0e-6, 'Inverse rotation not correct');
  testCase.verifyEqual(a.orientation, rotz(90)*eye(3), ...
    'Rotation not set');

  % Change the rotation back
  a = a.clearRotation();
  testCase.verifyEqual(double(a), diag([1, 2, 3, 4, 5, 6]), ...
    'AbsTol', 1.0e-6, 'Rotation not cleared with clear');

  % Apply rotation and set orientation
  a = a.rotateZ(pi/2);
  a = a.setOrientation(eye(3));
  testCase.verifyEqual(double(a), diag([1, 2, 3, 4, 5, 6]), ...
    'AbsTol', 1.0e-6, 'Rotation not cleared with orientation');

end
