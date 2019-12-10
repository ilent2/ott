function tests = testTmatrixDda
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../');
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

function testTooManyDipoles(testCase)

  wavelength_0 = 1;

  %Height of the cylinder in radius lengths - this is our data to compare

  c_height = 10.*wavelength_0;
  radius = 1*wavelength_0;

  %Refractive index of water is 1.3
  n_medium = 1.3;

  %Refractive index of glass is 1.5
  n_particle = 1.5;

  shape = ott.shapes.Shape.simple('cylinder', [radius, c_height]);
    
  testCase.verifyError(@() ott.TmatrixDda.simple(shape, ...
    'index_medium', n_medium,...
    'index_particle', n_particle, 'wavelength0', wavelength_0), ...
    'OTT:TmatrixDda:too_many_dipoles');

end

