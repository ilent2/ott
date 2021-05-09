function tests = testParticle
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testFixed(testCase)
  gui = ott.ui.particle.Fixed();
  testCase.addTeardown(@delete,gui);
end

function testFromShape(testCase)
  gui = ott.ui.particle.FromShape();
  testCase.addTeardown(@delete,gui);
end
