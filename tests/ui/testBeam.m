function tests = testBeam
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testAnnular(testCase)
  gui = ott.ui.beam.Annular();
  testCase.addTeardown(@delete,gui);
end

function testGaussian(testCase)
  gui = ott.ui.beam.Gaussian();
  testCase.addTeardown(@delete,gui);
end

function testPlaneWave(testCase)
  gui = ott.ui.beam.PlaneWave();
  testCase.addTeardown(@delete,gui);
end

function testPmParaxial(testCase)
  gui = ott.ui.beam.PmParaxial();
  testCase.addTeardown(@delete,gui);
end

function testScattered(testCase)
  gui = ott.ui.beam.Scattered();
  testCase.addTeardown(@delete,gui);
end

function testVisualise(testCase)
  gui = ott.ui.beam.Visualise();
  testCase.addTeardown(@delete,gui);
end
