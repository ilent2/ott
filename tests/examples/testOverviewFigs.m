function tests = testOverviewFigs
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../examples/packageOverview');

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

function testBscs(testCase)
  bscs();
end

function testDrag(testCase)
  drag();
end

function testParticles(testCase)
  particles();
end

function testShapes(testCase)
  shapes();
end

function testTmatrix(testCase)
  tmatrix();
end

