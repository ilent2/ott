function tests = testTmatrixDda
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../');
end

function testSphere(testCase)
  
  % Large error, unless we want a smoother sphere/slower runtime
  abs_tol = 2e-3;
  rel_tol = 22e-2;

  shape = ott.shapes.Sphere(0.1);

  nrel = 1.2;
  
  Tdda = ott.TmatrixDda.simple(shape, 'index_relative', nrel, ...
    'spacing', 1/35);
  Tmie = ott.TmatrixMie.simple(shape, 'index_relative', nrel);

  testCase.verifyEqual(Tdda.Nmax, Tmie.Nmax, ...
      'Nmax does not match Mie T-matrix');

  diag_dda = diag(Tdda.data);
  ndiag_dda = Tdda.data - diag(diag_dda);
  diag_mie = diag(Tmie.data);
  ndiag_mie = Tmie.data - diag(diag_mie);

  testCase.verifyEqual(full(ndiag_dda), full(ndiag_mie), ...
      'AbsTol', abs_tol, ...
      'T-matrix abstol failed');

  testCase.verifyEqual(full(diag_dda), full(diag_mie), ...
      'RelTol', rel_tol, ...
      'T-matrix retol failed');

end

function testDipoleSphere(testCase)

  xyz = [0;0;0];
  nrel = 1.2;
  
  radius = 0.01;
  d = (4*pi/3).^(1/3) .* radius;
  
  Tdda = ott.TmatrixDda(xyz, 'index_relative', nrel, ...
    'spacing', d, 'Nmax', 1, 'polarizability', 'CM');
  
  shape = ott.shapes.Sphere(radius);
  Tmie = ott.TmatrixMie.simple(shape, 'index_relative', nrel, ...
    'Nmax', 1);
  
  Tdda_diag = diag(Tdda.data);
  Tmie_diag = full(diag(Tmie.data));
  testCase.verifyEqual(Tdda_diag(4:end), Tmie_diag(4:end), ...
    'RelTol', 0.001, 'Upper Mie coefficeints not equal');
  testCase.verifyEqual(Tdda.data, full(Tmie.data), ...
    'AbsTol', 1e-6, 'Upper Mie coefficeints not equal');

end

function testSingleDipole(testCase)

  Nmax = 2;

  xyz = [0;0;0];
  nrel = 1.2;
  Tdda1 = ott.TmatrixDda(xyz, 'index_relative', nrel, ...
    'spacing', 1, 'Nmax', Nmax);
  
  nrel = [1.2; 1.2; 1.2];
  Tdda2 = ott.TmatrixDda(xyz, 'index_relative', nrel, ...
    'spacing', 1, 'Nmax', Nmax);
  testCase.verifyEqual(Tdda1.data, Tdda2.data);
  
  % Relative refractive index in Cartesian coordinates
  nrel = [1.2; 1.4; 1.5];
  Tdda3 = ott.TmatrixDda(xyz, 'index_relative', nrel, ...
    'spacing', 1, 'Nmax', Nmax);
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

function testHomogeneousMaterialSizes(testCase)
  % Test for the following part of the documentation
  %
  % The method supports homogenous and inhomogenous particles.
  % For homogeneous particles, specify the material as a scalar,
  % 3x1 vector or 3x3 polarizability matrix.

  shape = ott.shapes.Sphere(0.1);

  % Test homogeneous scalar
  nrel = 1.2;
  Tdda = ott.TmatrixDda.simple(shape, 'index_relative', nrel, ...
    'spacing', 1/20);
  
  nrel = [1.2; 1.3; 1.4];
  Tdda = ott.TmatrixDda.simple(shape, 'index_relative', nrel, ...
    'spacing', 1/20);
  
  nrel = [1; 1; 1];
  pol = [1.2, 0, 0; 0, 1.3, 0; 0, 0, 1.4];
  Tdda = ott.TmatrixDda.simple(shape, 'polarizability', pol, ...
    'spacing', 1/20, 'index_relative', nrel);
  
end
