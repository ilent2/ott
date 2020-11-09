function tests = testEllipsoid
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  radii = [1; 2; 3];
  shape = ott.shape.Ellipsoid(radii);

  testCase.verifyEqual(shape.radii, radii, 'radii');

  % test volume where sphere
  shape = ott.shape.Ellipsoid([1, 1, 1]);
  testCase.verifyEqual(shape.volume, 4/3*pi, 'AbsTol', 1.0e-14, 'volume');
  testCase.verifyEqual(shape.isSphere, true, 'isSphere');
  testCase.verifyEqual(shape.maxRadius, 1, 'maxRadius');
  testCase.verifyEqual(shape.xySymmetry, true, 'xy');
  testCase.verifyEqual(shape.zRotSymmetry, 0, 'xy');
  testCase.verifyEqual(shape.boundingBox, [-1,1;-1,1;-1,1], 'bb');

end

function testCasts(testCase)

  sph = ott.shape.Sphere();

  elp = ott.shape.Ellipsoid(sph);

  testCase.verifyEqual(elp.isSphere, true, 'isSphere');

  sph = ott.shape.Sphere(elp);

end

