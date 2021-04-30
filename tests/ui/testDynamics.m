function tests = testLauncher
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testIsolated(testCase)
  error('Not yet implemented');
end