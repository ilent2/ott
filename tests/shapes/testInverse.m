function tests = testInverse
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  a = ott.shapes.Cube();

  shape = ott.shapes.Inverse(a);

  testCase.verifyEqual(shape.internal, a, 'shape');

end

