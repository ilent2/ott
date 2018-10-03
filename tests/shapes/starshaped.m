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

