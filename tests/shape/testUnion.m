function tests = testUnion
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  a = ott.shape.Cube();
  b = ott.shape.Cylinder();

  shape = ott.shape.Union([a, b]);

  testCase.verifyEqual(shape.shapes, [a, b], 'shapes');

  c = ott.shape.Sphere();
  shape = shape | c;
  testCase.verifyEqual(shape.shapes, [a, b, c], 'shapes 3');

end

function testInsideXyz(testCase)

  radius = 1.0;
  offset1 = [0;0;-2];
  offset2 = [0;0;2];

  shape1 = ott.shape.Sphere(radius, 'position', offset1);
  shape2 = ott.shape.Sphere(radius, 'position', offset2);
  union = shape1 | shape2;

  xyz = [0;0;1] .* linspace(-5, 5, 20);
  inside1 = shape1.insideXyz(xyz, 'origin', 'global');
  inside2 = shape2.insideXyz(xyz, 'origin', 'global');
  insideU = union.insideXyz(xyz);

  testCase.verifyEqual(double(insideU), double(inside1 | inside2), ...
      'Union doesn''t match');
end

function testSurf(testCase)

  radius = 1.0;
  offset1 = [0;0;-2];
  offset2 = [0;0;2];

  shape1 = ott.shape.Sphere(radius, 'position', offset1);
  shape2 = ott.shape.Sphere(radius, 'position', offset2);
  union = shape1 | shape2;

  % Test functions from IsosurfSurf
  union.surfPoints();
  union.surf('visualise', false);
  union.normalsRtp([1;0;0]);
  union.normalsXyz([1;0;0]);

end

function testVolume(testCase)

  radius = 1.0;
  offset1 = [0;0;-2];
  offset2 = [0;0;2];

  shape1 = ott.shape.Sphere(radius, 'position', offset1);
  shape2 = ott.shape.Sphere(radius, 'position', offset2);
  union = shape1 | shape2;

  testCase.verifyEqual(union.volume, 2*shape1.volume, ...
      'RelTol', 1e-2, 'volume');

  union.xySymmetry = true;
  testCase.verifyEqual(union.volume, 2*shape1.volume, ...
      'RelTol', 1e-2, 'mirror');

  union.zRotSymmetry = 0;
  testCase.verifyEqual(union.volume, 2*shape1.volume, ...
      'RelTol', 1e-2, 'z0');

  union.zRotSymmetry = 2;
  testCase.verifyEqual(union.volume, 2*shape1.volume, ...
      'RelTol', 1e-2, 'z2');

  union.zRotSymmetry = 4;
  testCase.verifyEqual(union.volume, 2*shape1.volume, ...
      'RelTol', 1e-2, 'z4');

end

function testScale(testCase)

  a = ott.shape.Cube(1);
  b = ott.shape.Cube(1);

  a.position = [2;0;0];

  u = a | b;
  testCase.verifyEqual(u.volume, 2*a.volume, ...
      'RelTol', 1e-2, 'original volume');

  u = u * 0.5;
  testCase.verifyEqual(u.volume, a.volume/8 + b.volume/8, ...
      'RelTol', 1e-2, 'scaled');
  testCase.verifyEqual(u.shapes(1).position, [1;0;0], 'position');

end

