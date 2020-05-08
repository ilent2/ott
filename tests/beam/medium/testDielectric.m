function tests = testVacuum
  % Test for beam combination functionality
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testWater(testCase)

  vacuum = ott.beam.medium.Dielectric.Water;
  testCase.verifyEqual(vacuum.index, 1.33);
  testCase.verifyEqual(vacuum.relative_index, 1.33);
  
end
