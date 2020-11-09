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
  testCase.verifyEqual(shape.maxRadius, ...
      sqrt(radius^2 + (height/2)^2), 'maxRadius');
  testCase.verifyEqual(shape.volume, pi*radius.^2*height, 'volume');
  testCase.verifyEqual(shape.zRotSymmetry, 0, 'zrot');
  testCase.verifyEqual(shape.xySymmetry, true, 'xysym');
  testCase.verifyEqual(shape.boundingBox, [-1, 1; -1, 1; -1, 1]);
  testCase.verifyEqual(shape.perimeter, 4*radius+2*height, 'perimeter');

  % Scale
  shape = shape / 2;
  testCase.verifyEqual(shape.radius, radius/2, 'scale radius');
  testCase.verifyEqual(shape.height, height/2, 'scale height');

end

function testSurf(testCase)

  f = figure();
  testCase.addTeardown(@() close(f));
  shape = ott.shape.Cylinder();
  a = shape.surf('visualise', true, 'showNormals', true);

end

