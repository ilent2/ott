function tests = testPowerSpectrum
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testHandles(testCase)

  h = figure();
  testCase.addTeardown(@() close(h));

  xy = randn(2, 100);
  t = linspace(0, 1, 100);
  
  h = ott.tools.powerSpectrum(t, xy);
  testCase.verifySize(h, [2, 1], 'num handles');
  testCase.verifyClass(h, 'matlab.graphics.chart.primitive.Line', 'class');

end

