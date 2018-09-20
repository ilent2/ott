function tests = axisym
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
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

end
