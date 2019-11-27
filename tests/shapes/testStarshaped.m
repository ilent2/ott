function tests = axisym
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testVolume(testCase)

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.RelativeTolerance;
  tol = 1.0e-3;

  r = 1.0;
  shape = ott.shapes.Sphere(r);
  testCase.verifyThat(shape.volume, IsEqualTo(4/3*pi, ...
      'Within', RelativeTolerance(tol)), ...
      'Incorrect sphere volume');

  r = 1.0;
  shape = ott.shapes.Ellipsoid([r, r, r]);
  testCase.verifyThat(shape.volume, IsEqualTo(4/3*pi, ...
      'Within', RelativeTolerance(tol)), ...
      'Incorrect ellipsoid (sphere) volume');

  r = 1.0;
  shape = ott.shapes.Cylinder([r/2, r]);
  testCase.verifyThat(shape.volume, IsEqualTo(r*pi*(r/2).^2, ...
      'Within', RelativeTolerance(tol)), ...
      'Incorrect cylinder volume');

  r = 1.0;
  shape = ott.shapes.Cube([r]);
  testCase.verifyThat(shape.volume, IsEqualTo(1.0, ...
      'Within', RelativeTolerance(tol)), ...
      'Incorrect cube volume');

end

function testVoxelPositionsWithOffset(testCase)

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  tol = 1.0e-6;

  offset = [1;0;0];
  radius = 1.0;
  shape1 = ott.shapes.Sphere(radius);
  shape2 = ott.shapes.Sphere(radius, offset);

  spacing = 0.1;
  voxels1 = shape1.voxels(spacing, 'origin', 'world') + offset;
  voxels2 = shape2.voxels(spacing, 'origin', 'world');

  testCase.verifyThat(voxels2, IsEqualTo(voxels1, ...
      'Within', AbsoluteTolerance(tol)), ...
      'Incorrect voxel positions');
end

function testInsideXyz(testCase)

  import matlab.unittest.constraints.IsEqualTo;
  
  radius = 1.0;
  offset1 = [0;0;2];
  
  shape1 = ott.shapes.Sphere(radius, offset1);

  z = linspace(-10, 10, 100);
  mask = shape1.insideXyz(0, 0, z.', 'origin', 'world');
  
  testCase.verifyThat(size(mask), IsEqualTo(size(z.')), ...
    'Incorrect size of mask');

end
