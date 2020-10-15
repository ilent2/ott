function tests = testInverse
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  a = ott.shape.Cube();

  shape = ott.shape.Inverse(a);

  testCase.verifyEqual(shape.internal, a, 'shape');

end

