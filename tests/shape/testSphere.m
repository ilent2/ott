function tests = testSphere
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  radius = 1.5;
  shape = ott.shape.Sphere(radius);

  testCase.verifyEqual(shape.radius, radius, 'radius');
  testCase.verifyEqual(shape.volume, 4/3*pi*radius.^3, ...
      'AbsTol', 1.0e-14, 'vol');
  testCase.verifyEqual(shape.perimeter, 2*pi*radius, 'perimeter');
  testCase.verifyEqual(shape.isSphere, true, 'isSphere');
  testCase.verifyEqual(shape.xySymmetry, true, 'xysym');
  testCase.verifyEqual(shape.maxRadius, radius, 'maxRadius');

  % Scale
  shape = shape ./ 2;
  testCase.verifyEqual(shape.radius, radius./2, 'scale');

end

function testNormals(testCase)

  rtp = randn(3, 1);
  target = [1;0;0];
  shape = ott.shape.Sphere();
  testCase.verifyEqual(shape.normalsRtp(rtp), target, 'rtp');

  testCase.verifyEqual(shape.normalsXyz(rtp), rtp./vecnorm(rtp), 'xyz');

end

function testCasts(testCase)

  shape = ott.shape.Sphere();

  a = ott.shape.Ellipsoid(shape);
  testCase.verifyInstanceOf(a, 'ott.shape.Ellipsoid');

  b = ott.shape.Superellipsoid(shape);
  testCase.verifyInstanceOf(b, 'ott.shape.Superellipsoid');

end

