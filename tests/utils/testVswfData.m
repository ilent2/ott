function tests = testVswfData
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstruct(testCase)

  data = ott.utils.VswfData();

end

function testEvaluate(testCase)

  data = ott.utils.VswfData();
  ci = 1:8;
  theta = linspace(0, 5, 5);
  phi = linspace(0, 10, 5);
  [Y, Ydt, Ydp, data2] = data.evaluateYtp(ci, theta, phi);

  n = 1:5;
  kr = linspace(0, 5, 5);
  [hn, dhn, data2] = data.evaluateBessel(n, kr, 'incoming');
  [hn, dhn, data2] = data.evaluateBessel(n, kr, 'outgoing');
  [hn, dhn, data2] = data.evaluateBessel(n, kr, 'regular');

end

function testEmptyInputs(testCase)

  data = ott.utils.VswfData();
  ci = zeros(1, 0);
  theta = linspace(0, 5, 5);
  phi = linspace(0, 10, 5);
  [Y, Ydt, Ydp, data2] = data.evaluateYtp(ci, theta, phi);

end

