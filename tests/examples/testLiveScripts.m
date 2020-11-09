function tests = testLiveScripts
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../examples/liveScripts');
end

function testBeams(testCase)
  beams();
end

function testForce(testCase)
  force();
end

function testParticles(testCase)
  particles();
end

