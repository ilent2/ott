function tests = testTools
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testForcePosition(testCase)
  gui = ott.ui.tools.ForcePosition();
  testCase.addTeardown(@delete,gui);
end

function testPowerSpectrum(testCase)
  gui = ott.ui.tools.PowerSpectrum();
  testCase.addTeardown(@delete,gui);
end