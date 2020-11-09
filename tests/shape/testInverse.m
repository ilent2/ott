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
  testCase.verifyEqual(shape.volume, Inf, 'vol');
  testCase.verifyEqual(shape.maxRadius, Inf, 'R');
  testCase.verifyEqual(shape.starShaped, a.starShaped, 'star');
  testCase.verifyEqual(shape.zRotSymmetry, a.zRotSymmetry', 'z');
  testCase.verifyEqual(shape.xySymmetry, a.xySymmetry', 'xy');

end

function testInverseInverse(testCase)

  a = ott.shape.Cube();
  b = ~a;
  c = ~b;

  testCase.verifyEqual(c, a, 'cube');
end

