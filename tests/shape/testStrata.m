function tests = testStrata
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  shape = ott.shape.Strata();

  testCase.verifyEqual(shape.normal, [0; 0; 1], 'normal');

end

