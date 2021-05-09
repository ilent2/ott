function tests = testDrag
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testSimple(testCase)
  gui = ott.ui.drag.Simple();
  testCase.addTeardown(@delete,gui);
end