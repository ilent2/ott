function tests = testCylinder
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  radius = 1.0;
  height = 2.0;

  shape = ott.shape.Cylinder(radius, height);

  testCase.verifyEqual(shape.radius, radius, 'rad');
  testCase.verifyEqual(shape.height, height, 'height');
  testCase.verifyEqual(shape.volume, pi*radius.^2*height, 'volume');
  testCase.verifyEqual(shape.zRotSymmetry, 0, 'zrot');
  testCase.verifyEqual(shape.xySymmetry, true, 'xysym');
  testCase.verifyEqual(shape.boundingBox, [-1, 1; -1, 1; -1, 1]);
  testCase.verifyEqual(shape.perimeter, 4*radius+2*height, 'perimeter');

end

