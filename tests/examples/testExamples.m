function tests = testExamples
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../examples/');

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
  ottBeams;
end

function testDdaVaterite(testCase)
  ottDdaVaterite;
end

function testDynamics(testCase)
  ottDynamics;
end

function testForce(testCase)
  ottForce();
end

function testLandscape(testCase)
  ottLandscape();
end

function testNeuralNetwork(testCase)
  ottNeuralNetwok
end

function testNonSpherical(testCase)
  ottNonSpherical
end

function testParticles(testCase)
  ottParticles
end

function testWallEffects(testCase)
  ottWallEffect
end

