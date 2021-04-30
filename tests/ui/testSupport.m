function tests = testLauncher
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testPlaceholder(testCase)
  error('Not yet implemented');
end
