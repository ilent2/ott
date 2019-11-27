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
  
  xyz = [0;0;1] .* linspace(-10, 10, 100);
  inside1 = shape1.insideXyz(xyz(1, :).', xyz(2, :).', xyz(3, :).', 'origin', 'world');
  inside2 = shape2.insideXyz(xyz(1, :).', xyz(2, :).', xyz(3, :).', 'origin', 'world');
  insideU = union.insideXyz(xyz(1, :).', xyz(2, :).', xyz(3, :).');
  
  testCase.verifyThat(insideU, IsEqualTo(inside1 | inside2), ...
    'Union doesn''t match');
  
end
  