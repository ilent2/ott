function tests = testLiveScripts
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../examples/liveScripts');

  testCase.TestData.oldFigures = get(groot, 'Children');
end

function teardownOnce(testCase)

  % Clean up figures
  currentFigures = get(groot, 'Children');
  old = ismember(currentFigures, testCase.TestData.oldFigures);
  currentFigures(old) = [];
  close(currentFigures);

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

