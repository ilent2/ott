function tests = testShape
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  rotation = randn(3, 3);
  position = randn(3, 1);

  % Can't call constructor directly, create an empty
  shape = ott.shape.Empty('position', position, 'rotation', rotation);

  testCase.verifyEqual(shape.position, position, 'position');
  testCase.verifyEqual(shape.rotation, rotation, 'rotation');

  % Constructor also creates arrays

  rotation = randn(3, 6);
  position = randn(3, 2);

  % Can't call constructor directly, create an empty
  shape = ott.shape.Empty('position', position, 'rotation', rotation);

  testCase.verifyEqual(numel(shape), 2, 'size');
  testCase.verifyEqual([shape.position], position, 'arr position');
  testCase.verifyEqual([shape.rotation], rotation, 'arr rotation');

end

function testHetrogeneousArray(testCase)

  plane = ott.shape.Plane([0;0;1]);
  sphere = ott.shape.Sphere(1.0);

  array = [plane, sphere, sphere, plane];

  testCase.verifyEqual(size(array), [1, 4]);
  testCase.verifyClass(array(2:3), 'ott.shape.Sphere');

end

function testEmpty(testCase)

  sz = [0, 5];
  shape = ott.shape.Shape.empty(sz);
  testCase.verifyEqual(size(shape), sz);
end

function testInsideXyzHelper(testCase)

  radius = 1.0;
  shape = ott.shape.Sphere(radius);

  oxyz = [0;0;1];
  shape.position = [0;0;0];

  xyz = shape.insideXyzParseArgs(oxyz);
  testCase.verifyEqual(xyz, oxyz, 'Values changed (1 input)');

  shape.position = -[0;0;1];
  oxyz = [1;0;0];
  expected = [1;0;1];

  rtp = shape.insideXyzParseArgs(oxyz, 'origin', 'global');
  testCase.verifyEqual(rtp, expected, 'Values changed (1 input + shift)');

end

function testInsideRtpHelper(testCase)

  radius = 1.0;
  shape = ott.shape.Sphere(radius);

  ortp = [1;0;1];
  shape.position = [0;0;0];

  rtp = shape.insideRtpParseArgs(ortp);
  testCase.verifyEqual(rtp, ortp, 'Values changed (1 input)');

  shape.position = -[0;0;1];
  ortp = [1;0;0];
  expected = [2;0;0];

  rtp = shape.insideRtpParseArgs(ortp, 'origin', 'global');
  testCase.verifyEqual(rtp, expected, 'Values changed (1 input + shift)');

end

function testInsideXyz(testCase)

  radius = 1.0;
  shape = ott.shape.Sphere(radius);

  % Choose three points inside the shape and one outside
  b = [true, true, true, false];
  x = [0.5*radius.*rand(1, 3), 4.0];
  y = [0.5*radius.*rand(1, 3), 4.0];
  z = [0.5*radius.*rand(1, 3), 4.0];

  xyz = [x(:), y(:), z(:)].';

  testCase.verifyEqual(shape.insideXyz(xyz), b, ...
    'insideXyz with 1 argument failed');

  testCase.verifyEqual(shape.insideXyz(xyz, 'origin', 'global'), b, ...
    'insideXyz with 1 argument failed and optional arg');
end

function testIntersectBoundingBox(testCase)

  shape = ott.shape.Sphere(1.0);

  x0 = [0;0;-5];
  x1 = x0 + [0;0;1];
  ints = shape.intersectBoundingBox(x0, x1);
  testCase.verifyEqual(ints, [[0;0;-1],[0;0;1]], 'z ray on axis');

  x0 = [10;10;-5];
  x1 = x0 + [0;0;1];
  ints = shape.intersectBoundingBox(x0, x1);
  testCase.verifyEqual(ints, nan(3, 2), 'far away ray');

  x0 = [10;10;-5];
  x1 = x0 + [0;0;-1];
  ints = shape.intersectBoundingBox(x0, x1);
  testCase.verifyEqual(ints, nan(3, 2), 'leaving ray');

  x0 = [0;0;0];
  x1 = x0 + [0;0;-1];
  ints = shape.intersectBoundingBox(x0, x1);
  testCase.verifyEqual(ints, [[0;0;-1], nan(3, 1)], 'inside ray');
end

function testIntersect(testCase)

  radius = 1.0;
  x0 = zeros(3, 500);
  x1 = x0 + randn(3, 500);
  shape = ott.shape.Sphere(radius);
  ints = shape.intersect(x0, x1);

%   figure();
%   plot3(ints(1, :), ints(2, :), ints(3, :), '*');
%   axis equal;

  testCase.verifyEqual(vecnorm(ints), radius.*ones(1, size(ints, 2)), ...
    'AbsTol', 1.0e-4, 'Radius check');

end

function testGetBoundingBoxShape(testCase)

  shape = ott.shape.Cube();
  shape.position = randn(3, 1);

  bbshape = shape.getBoundingBoxShape();

  testCase.verifyEqual(bbshape.verts, shape.verts);

end


