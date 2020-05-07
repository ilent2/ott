function tests = shape
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testSimple(testCase)

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.IsOfClass;
  import matlab.unittest.constraints.RelativeTolerance;
  tol = 1.0e-3;

  r = 1.0;
  shape = ott.shapes.Shape.simple('sphere', r);
  testCase.verifyThat(shape, IsOfClass(?ott.shapes.Sphere), ...
      'Incorrect shape type');
  testCase.verifyThat(shape.perimiter, IsEqualTo(2*pi*r, ...
      'Within', RelativeTolerance(tol)), ...
      'Incorrect sphere radius');

end

function testArray(testCase)

  num_planes = 5;

  normals = randn(3, num_planes);
  planes = ott.shapes.Plane(normals);
  
  testCase.verifyEqual(size(planes), [1, num_planes], 'size');
  testCase.verifyEqual([planes.normal], normals, 'normals');
  
end

function testHetrogeneousArray(testCase)

  plane = ott.shapes.Plane([0;0;1]);
  sphere = ott.shapes.Sphere(1.0);
  
  array = [plane, sphere, sphere, plane];
  
  testCase.verifyEqual(size(array), [1, 4]);
  testCase.verifyClass(array(2:3), 'ott.shapes.Sphere');

end

function testEmpty(testCase)

  sz = [0, 5];
  shape = ott.shapes.Shape.empty(sz);
  testCase.verifyEqual(size(shape), sz);
end

function testInsideXyzHelper(testCase)

  radius = 1.0;
  shape = ott.shapes.Shape.simple('sphere', radius);
  
  oxyz = [0;0;1];
  shape.position = [0;0;0];
  
  xyz = shape.insideXyzParseArgs(oxyz);
  testCase.verifyEqual(xyz, oxyz, 'Values changed (1 input)');
  
  shape.position = -[0;0;1];
  oxyz = [1;0;0];
  expected = [1;0;1];
  
  rtp = shape.insideXyzParseArgs(oxyz, 'origin', 'world');
  testCase.verifyEqual(rtp, expected, 'Values changed (1 input + shift)');

end

function testInsideRtpHelper(testCase)

  radius = 1.0;
  shape = ott.shapes.Shape.simple('sphere', radius);
  
  ortp = [1;0;1];
  shape.position = [0;0;0];
  
  rtp = shape.insideRtpParseArgs(ortp);
  testCase.verifyEqual(rtp, ortp, 'Values changed (1 input)');
  
  shape.position = -[0;0;1];
  ortp = [1;0;0];
  expected = [2;0;0];
  
  rtp = shape.insideRtpParseArgs(ortp, 'origin', 'world');
  testCase.verifyEqual(rtp, expected, 'Values changed (1 input + shift)');

end

function testInsideXyz(testCase)

  radius = 1.0;
  shape = ott.shapes.Shape.simple('sphere', radius);
  
  % Choose three points inside the shape and one outside
  b = [true, true, true, false].';
  x = [0.5*radius.*rand(1, 3), 4.0];
  y = [0.5*radius.*rand(1, 3), 4.0];
  z = [0.5*radius.*rand(1, 3), 4.0];
  
  xyz = [x(:), y(:), z(:)].';
  
  testCase.verifyEqual(shape.insideXyz(xyz), b, ...
    'insideXyz with 1 argument failed');
  
  testCase.verifyEqual(shape.insideXyz(xyz, 'origin', 'world'), b, ...
    'insideXyz with 1 argument failed and optional arg');
end

function testIntersectBoundingBox(testCase)

  shape = ott.shapes.Sphere(1.0);
  
  vec = ott.utils.Vector([0;0;-5], [0;0;1]);
  ints = shape.intersectBoundingBox(vec);
  testCase.verifyEqual(ints, [0;0;-1;0;0;1], 'z ray on axis');
  
  vec = ott.utils.Vector([10;10;-5], [0;0;1]);
  ints = shape.intersectBoundingBox(vec);
  testCase.verifyEqual(ints, nan(6, 1), 'far away ray');

  vec = ott.utils.Vector([10;10;-5], [0;0;-1]);
  ints = shape.intersectBoundingBox(vec);
  testCase.verifyEqual(ints, nan(6, 1), 'leaving ray');
  
  vec = ott.utils.Vector([0;0;0], [0;0;-1]);
  ints = shape.intersectBoundingBox(vec);
  testCase.verifyEqual(ints, [nan(3, 1); 0;0;-1], 'inside ray');
end

function testIntersect(testCase)

  radius = 1.0;
  vecs = ott.utils.Vector([0;0;0], randn(3, 500));
  shape = ott.shapes.Sphere(radius);
  ints = shape.intersect(vecs);
  
%   figure();
%   plot3(ints(1, :), ints(2, :), ints(3, :), '*');
%   axis equal;
  
  testCase.verifyEqual(vecnorm(ints), radius.*ones(1, size(ints, 2)), ...
    'AbsTol', 1.0e-4, 'Radius check');

end
