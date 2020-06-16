function tests = testCube
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  width = 1.0;

  shape = ott.shapes.Cube(width);

  testCase.verifyEqual(shape.width, width, 'width');
  testCase.verifyEqual(shape.starShaped, true, 'star');
  testCase.verifyEqual(shape.volume, width.^3, 'volume');
  testCase.verifyEqual(shape.xySymmetry, true, 'xysym');
  testCase.verifyEqual(shape.zRotSymmetry, 4, 'zsym');
  testCase.verifyEqual(shape.maxRadius, sqrt(3*width.^2./4), ...
      'AbsTol', 1.0e-10, 'maxRadius');

end

