function tests = testSlab
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  normal = [0;0;1];
  depth = 0.5;

  slab = ott.shape.Slab(depth, normal);

  % Check properties
  testCase.verifyEqual(slab.normal, normal);
  testCase.verifyEqual(slab.depth, depth);
  testCase.verifyEqual(slab.maxRadius, Inf);
  testCase.verifyEqual(slab.volume, Inf);
  testCase.verifyEqual(slab.zRotSymmetry, 0, 'zrot');
  testCase.verifyEqual(slab.xySymmetry, true, 'xySym');
  testCase.verifyEqual(slab.starShaped, true, 'starShaped');
  testCase.verifyEqual(slab.boundingBox, [-Inf, Inf; -Inf, Inf; -0.25, 0.25], 'bb');
  
  % Scale
  slab = slab/2;
  testCase.verifyEqual(slab.normal, normal, 'snormal');
  testCase.verifyEqual(slab.depth, depth/2, 'sdepth');

end

function testIntersect(testCase)

  slab = ott.shape.Slab(1);
  
  testCase.verifyEqual(slab.intersect([0;0;0], [10;0;0]), ...
      nan(3, 1), 'nan');
  testCase.verifyEqual(slab.intersectAll([0;0;-10], [0;0;10]), ...
      [0,0;0,0;-0.5,0.5], 'all');

end

function testInsideXyz(testCase)

  slab = ott.shape.Slab();
  
  xyz = [0, 0, 0; 0, 0, 0; 1, 0.25, -1];
  b = slab.insideXyz(xyz);
  testCase.verifyEqual(b, [false, true, false]);

end

function testNormalsXyz(testCase)

  slab = ott.shape.Slab();

  testCase.verifyEqual(slab.normalsXyz([0;0;1]), [0;0;1], 'pos');
  testCase.verifyEqual(slab.normalsXyz([0;0;-1]), [0;0;-1], 'neg');

end

function testSurf(testCase)

  plane = ott.shape.Slab();

  h = figure();
  testCase.addTeardown(@close, h);

  p = plane.surf('scale', 2.0);
  testCase.verifyClass(p, 'matlab.graphics.primitive.Patch');

end

