function tests = testTmatrix
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testDda(testCase)
  gui = ott.ui.tmatrix.Dda();
  testCase.addTeardown(@delete,gui);
end

function testMie(testCase)
  gui = ott.ui.tmatrix.Mie();
  testCase.addTeardown(@delete,gui);
end

function testPointmatch(testCase)
  gui = ott.ui.tmatrix.Pointmatch();
  testCase.addTeardown(@delete,gui);
end

function testSimple(testCase)
  gui = ott.ui.tmatrix.Simple();
  testCase.addTeardown(@delete,gui);
end
