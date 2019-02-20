function tests = testTmatrixDda
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testSphere(testCase)

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  
  % Large error, unless we want a smoother sphere/slower runtime
  tol = 0.5e-2;

  shape = ott.shapes.Sphere(0.1);

  nrel = 1.2;
  
  Tdda = ott.TmatrixDda.simple(shape, 'index_relative', nrel, ...
    'spacing', 1/50);
  Tmie = ott.TmatrixMie.simple(shape, 'index_relative', nrel);

  testCase.verifyThat(Tdda.Nmax, IsEqualTo(Tmie.Nmax), ...
      'Nmax does not match Mie T-matrix');

  testCase.verifyThat(full(Tdda.data), IsEqualTo(full(Tmie.data), ...
      'Within', AbsoluteTolerance(tol)), ...
      'T-matrix does not match Mie T-matrix within tolerance');

end

function testCube(testCase)
  % Test DDA against a cube calculated using point matching

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;

  % Large error, unless we want a smoother sphere/slower runtime
  tol = 0.5e-2;

  shape = ott.shapes.Shape.simple('cube', 0.1);

  nrel = 1.2;

  Tdda = ott.TmatrixDda.simple(shape, 'index_relative', nrel, ...
    'spacing', 1/50);
  Tpm = ott.TmatrixPm.simple(shape, 'index_relative', nrel);

  testCase.verifyThat(full(Tdda.data), IsEqualTo(full(Tpm.data), ...
      'Within', AbsoluteTolerance(tol)), ...
      'T-matrix does not match point matching T-matrix within tolerance');

end

