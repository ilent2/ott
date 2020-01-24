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

