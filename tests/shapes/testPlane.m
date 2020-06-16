function tests = testPlane
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  normal = [0;0;1];
  offset = 0.0;

  shape = ott.shapes.Plane(normal, 'offset', offset);

  testCase.verifyEqual(shape.normal, normal, 'normal');
  testCase.verifyEqual(shape.offset, offset, 'normal');
  testCase.verifyEqual(shape.xySymmetry, false, 'xysym');
  testCase.verifyEqual(shape.zRotSymmetry, 0, 'zsym');
  testCase.verifyEqual(shape.starShaped, false, 'star');
  testCase.verifyEqual(shape.maxRadius, Inf, 'maxR');
  testCase.verifyEqual(shape.volume, Inf, 'volume');

end

function testInsideXyz(testCase)

  plane = ott.shapes.Plane();

  xyz = [0, 0; 0, 0; 1, -1];
  b = plane.insideXyz(xyz);
  testCase.verifyEqual(b, [true, false]);

end

function testSurf(testCase)

  normal = rand(3, 1);
  offset = 0.0;

  plane = ott.shapes.Plane(normal, 'offset', offset);

  h = figure();

  p = plane.surf('scale', 2.0);
  testCase.verifyClass(p, 'matlab.graphics.primitive.Patch');

  close(h);

end

function testIntersect(testCase)

  normal = [0;0;1];
  offset = 0.0;

  plane = ott.shapes.Plane(normal, 'offset', offset);

  ray = ott.utils.Vector([0;0;-1], [1;0;0]);
  locs = plane.intersect(ray);
  testCase.verifyEqual(locs, nan(3, 1), 'parallel ray');

  ray = ott.utils.Vector([0;0;-1], [0;0;1]);
  locs = plane.intersect(ray);
  testCase.verifyEqual(locs, [0;0;0], 'nice ray');

  ray = ott.utils.Vector([0;0;-1], [0;0;-1]);
  locs = plane.intersect(ray);
  testCase.verifyEqual(locs, nan(3, 1), 'negative ray');

end


