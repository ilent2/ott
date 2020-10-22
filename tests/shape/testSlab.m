function tests = testSlab
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  normal = [0;0;1];
  depth = 0.5;

  slab = ott.shape.Slab(normal, depth);

  % Check properties
  testCase.verifyEqual(slab.normal, normal);
  testCase.verifyEqual(slab.depth, depth);
  testCase.verifyEqual(slab.maxRadius, Inf);
  testCase.verifyEqual(slab.volume, Inf);

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

