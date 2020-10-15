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

