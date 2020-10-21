function tests = testSpharmData
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstruct(testCase)

  data = ott.utils.SpharmData();
  testCase.verifyEmpty(data.uphi, 'uphi');
  testCase.verifyEmpty(data.utheta, 'utheta');
  testCase.verifyEmpty(data.um, 'um');
  testCase.verifyEmpty(data.uci, 'uci');

  testCase.verifyEmpty(data.expimphi, 'expimphi');
  testCase.verifyEmpty(data.Y, 'Y');
  testCase.verifyEmpty(data.Ytheta, 'dtY');
  testCase.verifyEmpty(data.Yphi, 'dpY');

end

function testEvaluate(testCase)

  data = ott.utils.SpharmData();
  ci = 1:8;
  theta = linspace(0, 5, 5);
  phi = linspace(0, 10, 5);
  [Y, Ydt, Ydp, data2] = data.evaluate(ci, theta, phi);

  testCase.verifySize(Y, [8, 5], 'y size');
  testCase.verifySize(Ydt, [8, 5], 'ydt size');
  testCase.verifySize(Ydp, [8, 5], 'ydp size');

  testCase.verifyEqual(data2.uci, 1:8, 'data2.uci');
  testCase.verifyEqual(data2.um, -2:2, 'data2.um');

end

function testValues(testCase)

  data = ott.utils.SpharmData();
  ci = [1, 2]+3;
  theta = [0.1, 0.3];
  phi = [0.2, 0.1];
  [n, m] = ott.utils.combined_index(ci);
  [Y, Ydt, Ydp, data] = data.evaluate(ci, theta, phi);
  [Ye, Ydte, Ydpe] = ott.utils.spharm(unique(n), m, theta, phi);

  testCase.verifyEqual(Y, Ye.', 'AbsTol', 1e-15, 'Y');
  testCase.verifyEqual(Ydt, Ydte.', 'AbsTol', 1e-15, 'Yde');
  testCase.verifyEqual(Ydp, Ydpe.', 'AbsTol', 1e-15, 'Ydp');

  % Test second retrieval
  theta = [theta, 0.4];
  phi = [phi, 0.6];
  ci = [1, 2, 3, 4]+3;
  [n, m] = ott.utils.combined_index(ci);
  [Y, Ydt, Ydp, ~] = data.evaluate(ci, theta, phi);
  [Ye, Ydte, Ydpe] = ott.utils.spharm(unique(n), m, theta, phi);

  testCase.verifyEqual(Y, Ye.', 'AbsTol', 1e-15, 'Y 2');
  testCase.verifyEqual(Ydt, Ydte.', 'AbsTol', 1e-15, 'Yde 2');
  testCase.verifyEqual(Ydp, Ydpe.', 'AbsTol', 1e-15, 'Ydp 2');

end

