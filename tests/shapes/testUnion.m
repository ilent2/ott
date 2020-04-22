function tests = testUnion
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testVoxels(testCase)

  radius = 1.0;
  offset1 = [0;0;-2];
  offset2 = [0;0;2];
  
  shape1 = ott.shapes.Sphere(radius, offset1);
  shape2 = ott.shapes.Sphere(radius, offset2);
  union = ott.shapes.Union([shape1, shape2]);
  
  xyz = union.voxels();
end

function testInsideXyz(testCase)

  import matlab.unittest.constraints.IsEqualTo;
  
  radius = 1.0;
  offset1 = [0;0;-2];
  offset2 = [0;0;2];
  
  shape1 = ott.shapes.Sphere(radius, offset1);
  shape2 = ott.shapes.Sphere(radius, offset2);
  union = ott.shapes.Union([shape1, shape2]);
  
  xyz = [0;0;1] .* linspace(-5, 5, 20);
  inside1 = shape1.insideXyz(xyz, 'origin', 'world');
  inside2 = shape2.insideXyz(xyz, 'origin', 'world');
  insideU = union.insideXyz(xyz);
  
  testCase.verifyEqual(double(insideU), double(inside1 | inside2), ...
    'Union doesn''t match');
  
end
  