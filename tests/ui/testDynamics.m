function tests = testDynamics
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testIsolated(testCase)
  gui = ott.ui.dynamics.Isolated();
  testCase.addTeardown(@delete,gui);
end