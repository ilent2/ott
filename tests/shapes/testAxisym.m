function tests = axisym
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testPerimiter(testCase)

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  import matlab.unittest.constraints.RelativeTolerance;
  tol = 1.0e-3;

  r = 1.0;
  shape = ott.shapes.Sphere(r);
  testCase.verifyThat(shape.perimiter, IsEqualTo(2*pi*r, ...
      'Within', RelativeTolerance(tol)), ...
      'Incorrect sphere radius');
    
  r = 1.0;
  shape = ott.shapes.Ellipsoid([r, r, r]);
  testCase.verifyThat(shape.perimiter, IsEqualTo(2*pi*r, ...
      'Within', RelativeTolerance(tol)), ...
      'Incorrect ellipsoid (sphere) radius');
    
  r = 1.0;
  shape = ott.shapes.Cylinder([r/2, r]);
  testCase.verifyThat(shape.perimiter, IsEqualTo(4*r, ...
      'Within', RelativeTolerance(tol)), ...
      'Incorrect cube radius');

end

function testBoundarypoints(testCase)

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  tol = 1.0e-3;

  shape = ott.shapes.Sphere([1.0]);

  numr = 100;
  [rtp, n, ds] = shape.boundarypoints('npts', numr);

  testCase.verifyThat(length(rtp(:, 1)), IsEqualTo(numr), ...
      'Incorrect number of values returned');

  testCase.verifyThat(rtp(:, 1), IsEqualTo(ones(size(rtp(:, 1))), ...
      'Within', AbsoluteTolerance(tol)), ...
      'Incorrect values for r');

  testCase.verifyThat(n(:, 1), IsEqualTo(ones(size(n(:, 1))), ...
      'Within', AbsoluteTolerance(tol)), ...
      'Incorrect values for r');
    
  shape = ott.shapes.Ellipsoid([1.0, 1.0, 1.0]);

  numr = 100;
  [rtp, n, ds] = shape.boundarypoints('npts', numr);

  testCase.verifyThat(length(rtp(:, 1)), IsEqualTo(numr), ...
      'Ellipsoid: Incorrect number of values returned');

  testCase.verifyThat(rtp(:, 1), IsEqualTo(ones(size(rtp(:, 1))), ...
      'Within', AbsoluteTolerance(tol)), ...
      'Ellipsoid: Incorrect values for r');

  testCase.verifyThat(n(:, 1), IsEqualTo(ones(size(n(:, 1))), ...
      'Within', AbsoluteTolerance(tol)), ...
      'Ellipsoid: Incorrect values for r');

end
