function tests = testShape
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testCadFile(testCase)
  gui = ott.ui.shape.CadFileLoader();
  testCase.addTeardown(@delete,gui);
end

function testSimple(testCase)
  gui = ott.ui.shape.Simple();
  testCase.addTeardown(@delete,gui);
end

function testVisualise(testCase)
  gui = ott.ui.shape.Visualise();
  testCase.addTeardown(@delete,gui);
end
