function tests = testTmatrixDda
  % Unit tests for TmatrixDDA.  A lot of these tests are for hidden
  % functions that are likely to change in the next release (they
  % will probably move to other locations).  We had a bit of trouble
  % getting the signs correct.

  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../');
end

function testSphere(testCase)
  
  % Large error, unless we want a smoother sphere/slower runtime
  abs_tol = 1.0e-4;
  rel_tol = 3.0e-2;
  
  wavelength0 = 1.0e-6;

  shape = ott.shapes.Sphere(0.2*wavelength0);

  nrel = 1.2;
  
  % Small changes to spacing spacing seems to have a very large
  % effect on error, even/odd spacing for voxels also change things
%   spacing = 1/35;
  spacing = 1/35 .* wavelength0;
  
  Tdda = ott.TmatrixDda.simple(shape, 'index_relative', nrel, ...
    'spacing', spacing, 'index_medium', 1.0, 'wavelength0', wavelength0);
  Tmie = ott.TmatrixMie.simple(shape, 'index_relative', nrel, ...
    'index_medium', 1.0, 'wavelength0', wavelength0);
  
%   z = [0;0;1] .* linspace(-3, 3, 100);
%   beam = ott.BscPmGauss();
%   fz_mie = ott.forcetorque(beam, Tmie, 'position', z);
%   fz_dda = ott.forcetorque(beam, Tdda, 'position', z);
%   figure(), plot(z(3, :), [fz_mie(3, :); fz_dda(3, :)]);
%   legend({'Mie', 'DDA'});

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
      'RelTol', rel_tol, 'AbsTol', 1.0e-3, ...
      'T-matrix retol failed');

end

function testDipoleSphere(testCase)

  xyz = [0;0;0];
  nrel = 1.2;
  
  radius = 0.01;
  d = (4*pi/3).^(1/3) .* radius;
  
  Tdda = ott.TmatrixDda(xyz, 'index_relative', nrel, ...
    'spacing', d, 'Nmax', 1, 'polarizability', 'CM', ...
    'z_rotational_symmetry', 1, 'z_mirror_symmetry', false);
  
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
  spacing = 1/20;
  Tdda1 = ott.TmatrixDda.simple(shape, 'index_relative', nrel, ...
    'spacing', spacing);
  
  nrel = [1.2; 1.3; 1.4];
  pol = ott.utils.polarizability.LDR(spacing, nrel);
  pol = diag(pol);
  Tdda2 = ott.TmatrixDda.simple(shape, 'polarizability', pol, ...
    'spacing', spacing, 'index_relative', nrel);
  
  testCase.verifyEqual(Tdda2.data, Tdda1.data, ...
    'AbsTol', 1e-16, ...
    'Polarizability and index_relative not equivilant');
  
end

function testInhomogeneousPolarisability(testCase)

  shape = ott.shapes.Sphere(0.1);
  
  n_medium = 1.1;
  
  nrel = [1.2; 1.3; 1.4];
  spacing = 1/20;
  Tdda1 = ott.TmatrixDda.simple(shape, 'index_relative', nrel, ...
    'spacing', spacing, 'index_medium', n_medium, ...
    'wavelength0', 1.0);
  
  Ndipoles = 32;
  
  nrel = [1.2; 1.3; 1.4];
  pol = ott.utils.polarizability.LDR(spacing.*n_medium, nrel);
  pol = repmat(diag(pol), [1, Ndipoles]);
  nrel = repmat(diag(nrel), [1, Ndipoles]);
  Tdda2 = ott.TmatrixDda.simple(shape, 'polarizability', pol, ...
    'spacing', spacing, 'index_relative', nrel, ...
    'index_medium', n_medium, ...
    'wavelength0', 1.0);
  
  testCase.verifyEqual(Tdda2.data, Tdda1.data, ...
    'AbsTol', 1e-15, ...
    'Polarizability and index_relative not equivilant');

end

function testSphereRotSym(testCase)
  % Create a sphere T-matrix with z-axis rotational symmetry
  
  % Large error, unless we want a smoother sphere/slower runtime
  abs_tol = 1e-5;
  rel_tol = 22.2e-2;

  nrel = 1.2;
  shape = ott.shapes.Sphere(0.1);

  Tmie = ott.TmatrixMie.simple(shape, 'index_relative', nrel);
  Nmax = Tmie.Nmax(1);
  
  z_rotational_symmetry = 4;
  z_mirror_symmetry = false;
  
  spacing = 1/40;
  voxels = shape.voxels(spacing, 'even_range', true);

%   voxels = ott.utils.roty(35) * voxels;
  
  Tdda = ott.TmatrixDda(voxels, ...
      'index_relative', nrel, ...
      'Nmax', Nmax, ...
      'spacing', spacing, ...
      'z_mirror_symmetry', z_mirror_symmetry, ...
      'z_rotational_symmetry', z_rotational_symmetry);

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

function testSphereMirSym(testCase)
  % Create a sphere T-matrix with z-axis rotational symmetry
  
  % Large error, unless we want a smoother sphere/slower runtime
  abs_tol = 2e-3;
  rel_tol = 22.2e-2;

  nrel = 1.2;
  shape = ott.shapes.Sphere(0.1);

  Tmie = ott.TmatrixMie.simple(shape, 'index_relative', nrel);
  Nmax = Tmie.Nmax(1);
  
  z_rotational_symmetry = 1;
  z_mirror_symmetry = true;
  
  spacing = 1/40;
  voxels = shape.voxels(spacing, 'even_range', true);

  Tdda = ott.TmatrixDda(voxels, ...
      'index_relative', nrel, ...
      'Nmax', Nmax, ...
      'spacing', spacing, ...
      'z_mirror_symmetry', z_mirror_symmetry, ...
      'z_rotational_symmetry', z_rotational_symmetry);

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

function testSphereRotMirSym(testCase)
  % Create a sphere T-matrix with z-axis rotational symmetry
  
  % Large error, unless we want a smoother sphere/slower runtime
  abs_tol = 2e-3;
  rel_tol = 22.2e-2;

  nrel = 1.2;
  shape = ott.shapes.Sphere(0.1);

  Tmie = ott.TmatrixMie.simple(shape, 'index_relative', nrel);
  Nmax = Tmie.Nmax(1);
  
  z_rotational_symmetry = 4;
  z_mirror_symmetry = true;
  
  spacing = 1/40;
  voxels = shape.voxels(spacing, 'even_range', true);

  Tdda = ott.TmatrixDda(voxels, ...
      'index_relative', nrel, ...
      'Nmax', Nmax, ...
      'spacing', spacing, ...
      'z_mirror_symmetry', z_mirror_symmetry, ...
      'z_rotational_symmetry', z_rotational_symmetry);

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

function testSphere3Modes(testCase)
  
  % Large error, unless we want a smoother sphere/slower runtime
  abs_tol = 2e-3;
  rel_tol = 22.2e-2;
  
  wavelength0 = 1.0e-6;

  shape = ott.shapes.Sphere(0.1*wavelength0);

  nrel = 1.2;
  
  % Small changes to spacing spacing seems to have a very large
  % effect on error, even/odd spacing for voxels also change things
%   spacing = 1/35;
  spacing = 1/40 .* wavelength0;
  
  Tdda = ott.TmatrixDda.simple(shape, 'index_relative', nrel, ...
    'spacing', spacing, 'index_medium', 1.0, 'wavelength0', wavelength0, ...
    'modes', (1:3).');
  Tmie = ott.TmatrixMie.simple(shape, 'index_relative', nrel, ...
    'index_medium', 1.0, 'wavelength0', wavelength0);

  testCase.verifyEqual(Tdda.Nmax, Tmie.Nmax, ...
      'Nmax does not match Mie T-matrix');

  diag_dda = diag(Tdda.data);
  ndiag_dda = Tdda.data - diag(diag_dda);
  diag_mie = diag(Tmie.data);
  ndiag_mie = Tmie.data - diag(diag_mie);

  testCase.verifyEqual(full(ndiag_dda), full(ndiag_mie), ...
      'AbsTol', abs_tol, ...
      'T-matrix abstol failed');

  testCase.verifyEqual(full(diag_dda([1:3, 25:27])), full(diag_mie([1:3, 25:27])), ...
      'RelTol', rel_tol, ...
      'T-matrix retol failed');

  testCase.verifyEqual(full(diag_dda([4:24, 28:end])), 0.*full(diag_mie([4:24, 28:end])), ...
      'Uncalculated T-matrix modes incorrect');

end

function testDipoleFarfield(testCase)

  distance = 1e6;
  dipole_xyz = [0;1;0];
  targets_xyz = ott.utils.rtp2xyz([distance,0,0]).';
  M_dipole = eye(3);
  k = 2*pi;
  
  % Phase factor for near-field
  P = exp(1i*k*distance)./distance;
  
  farfield = ott.TmatrixDda.dipole_farfield(dipole_xyz, targets_xyz, M_dipole, k);
  nearfield = ott.TmatrixDda.dipole_nearfield(dipole_xyz, targets_xyz, M_dipole, k);

  testCase.verifyEqual(farfield, nearfield./P, 'AbsTol', 1.0e-3);
  
end

function testFieldModes(testCase)

  % Check near and far-field match in far-field limit
  distance = 1.0e6;
  k = 1;
  kr = k*distance;
  theta = pi/2;
  phi = pi/2;
  Nmax = 1;
  
  % Phase factor for near-field
  P = exp(1i*k*distance)./kr;
  
  farfield = ott.TmatrixDda.calculate_farfield_modes(theta, phi, Nmax);
  nearfield = ott.TmatrixDda.calculate_nearfield_modes(kr, theta, phi, Nmax);
  
  testCase.verifyEqual(farfield, nearfield./P, 'AbsTol', 1.0e-3, ...
    'near/far limit');
  
  % Test against some actual values
  total_modes = ott.utils.combined_index(Nmax, Nmax);
  ab = randn(2*total_modes, 1);
  beam = ott.Bsc(ab(1:end/2), ab(end/2+1:end), ...
    'outgoing', 'incident', 'k_medium', k);
  
  % Generate far-field grid of points for field evaluation
  [X, Y, Z] = sphere();
  rtp = ott.utils.xyz2rtp([X(:), Y(:), Z(:)]);

  Etmatrix = beam.farfield(rtp(:, 2), rtp(:, 3));
  Etmatrix = Etmatrix ./ max(abs(Etmatrix(:)));
  
  farfield = ott.TmatrixDda.calculate_farfield_modes(rtp(:, 2), rtp(:, 3), Nmax);
  farfield = farfield * ab;
  farfield = ott.utils.xyzv2rtpv(reshape(farfield, 3, []).', [X(:), Y(:), Z(:)]).';
  farfield = farfield ./ max(abs(farfield(:)));
  
  testCase.verifyEqual(farfield, Etmatrix, 'AbsTol', 1.0e-15, 'farfield');

  Etmatrix = beam.emFieldXyz([X(:), Y(:), Z(:)].');
  Etmatrix = Etmatrix ./ max(abs(Etmatrix(:)));
  
  nearfield = ott.TmatrixDda.calculate_nearfield_modes(k*rtp(:, 1), rtp(:, 2), rtp(:, 3), Nmax);
  nearfield = nearfield * ab;
  nearfield = reshape(nearfield, 3, []);
  nearfield = nearfield ./ max(abs(nearfield(:)));
  
  testCase.verifyEqual(nearfield, Etmatrix, 'AbsTol', 1.0e-15, 'nearfield');
end

function testCart2Sph(testCase)

  theta = 0.0;
  phi = 0.0;

  M = ott.TmatrixDda.cart2sph_mat(theta, phi);

  Md = ott.TmatrixDda.sph2cart_mat(theta, phi);
  
  testCase.verifyEqual(M * Md, eye(3), 'AbsTol', 1.0e-15);
  
end

function testFarfield(testCase)

  wavelength0 = 1.0e-6;
  shape = ott.shapes.Sphere(0.3*wavelength0);
  nrel = 1.2;
  
  % Small changes to spacing spacing seems to have a very large
  % effect on error, even/odd spacing for voxels also change things
  spacing = 1/20 .* wavelength0;
  
  TddaFar = ott.TmatrixDda.simple(shape, 'index_relative', nrel, ...
    'spacing', spacing, 'index_medium', 1.0, 'wavelength0', wavelength0, ...
    'use_nearfield', false);
  TddaNear = ott.TmatrixDda.simple(shape, 'index_relative', nrel, ...
    'spacing', spacing, 'index_medium', 1.0, 'wavelength0', wavelength0, ...
    'use_nearfield', true);
  
  testCase.verifyEqual(TddaFar.data, TddaNear.data, 'AbsTol', 1.0e-6);

end

function testDipoleNearfield(testCase)
  % Test fields from dipole_nearfield are outgoing
  
  offset = 0.5;
  distance1 = linspace(10, 11, 100).';
  distance2 = distance1 + offset;
  dipole_xyz = [0;0;0];
  tp = [pi/2, pi/2] .* ones(100, 2);
  rtp1 = [distance1,tp];
  rtp2 = [distance2,tp];
  targets_xyz1 = ott.utils.rtp2xyz(rtp1).';
  targets_xyz2 = ott.utils.rtp2xyz(rtp2).';
  M_dipole = eye(3);
  k = 2*pi;
  
  % Different sign convention for modes
  nearfield1 = ott.TmatrixDda.dipole_nearfield(dipole_xyz, targets_xyz1, M_dipole, k);
  nearfield1 = nearfield1(1:3:end, 1) .* distance1;
  
  nearfield2 = ott.TmatrixDda.dipole_nearfield(dipole_xyz, targets_xyz2, M_dipole, k);
  nearfield2 = nearfield2(1:3:end, 1) .* distance2;
  
  testCase.verifyEqual(real(nearfield1.*exp(-1i*offset*k)), real(nearfield2), ...
      'RelTol', 1.0e-2, 'AbsTol', 3.0e-2);

end

function testOneDipoleDirect(testCase)

  D = [0;0;0];  rtp = [0,0,0];
  k = 2*pi;
  spacing = 0.01;
  n_relative = 1.2;
  alpha = ott.utils.polarizability.LDR(spacing, n_relative);
  inv_alpha = ott.TmatrixDda.alpha_to_full_inv_alpha(alpha, size(rtp, 1));
  
  % Generate far-field grid of points for field evaluation
  [X, Y, Z] = sphere();
  rtpFar = ott.utils.xyz2rtp([X(:), Y(:), Z(:)]);
  
  % Generate a dipole with Mie
  Tmatrix = ott.TmatrixMie(0.005, 'index_relative', n_relative);
  Tmatrix.Nmax = 1;
  
  % Generate a plane wave beam
  beam = ott.BscPlane(0, 0, 'Nmax', 10);
  
  % Calculate Tmatrix far-field
  sbeam = Tmatrix * beam;
  sbeam = sbeam.scatteredField(beam);
  Etmatrix = sbeam.farfield(rtpFar(:, 2), rtpFar(:, 3));
  
  % Calculate beam fields
  Ei = beam.emFieldXyz([0;0;0]);
  
  % Calculate interatction matrix
  A = ott.TmatrixDda.interaction_A_total(k, rtp, inv_alpha, false, 1);
  
  %  Calcualte far-field matrix
  M_dipole = eye(3);
  F = ott.TmatrixDda.dipole_farfield(D, [X(:), Y(:), Z(:)].', M_dipole, k);

  % Calculat dipole fields
  Edipole = ott.TmatrixDda.solve_and_evaluate(A, F, Ei(:), false);
  Edipole = reshape(Edipole, 3, []);
  
  % Change Edipole to spherical coordinates
  Edipole = ott.utils.xyzv2rtpv(Edipole.', [X(:), Y(:), Z(:)]).';
  
  % Likely diffeernce in magnitude (due to different size/polarisability)
  % But both should be dipole shaped and same direction
  maxTmatrix = max(abs(Etmatrix(:)));
  maxDipole = max(abs(Edipole(:)));
  
  testCase.verifyEqual(Edipole ./ maxDipole, Etmatrix ./ maxTmatrix, ...
      'AbsTol', 1.0e-3, 'farfield');
    
  % Repeat calculation with near-field
  F = ott.TmatrixDda.dipole_nearfield(D, [X(:), Y(:), Z(:)].', M_dipole, k);
  
  Etmatrix = sbeam.emFieldXyz([X(:), Y(:), Z(:)].');
  
  % Calculat dipole fields
  Edipole = ott.TmatrixDda.solve_and_evaluate(A, F, Ei(:), false);
  Edipole = reshape(Edipole, 3, []);
  
  % Likely diffeernce in magnitude (due to different size/polarisability)
  % But both should be dipole shaped and same direction
  maxTmatrix = max(abs(Etmatrix(:)));
  maxDipole = max(abs(Edipole(:)));
  
  Edipole = Edipole ./ maxDipole;
  Etmatrix = Etmatrix ./ maxTmatrix;
  
  testCase.verifyEqual(Edipole, Etmatrix, ...
      'AbsTol', 1.0e-3, 'nearfield');
  
%   idx = 1;
%   
%   figure();
%   subplot(1, 2, 1);
%   surf(X, Y, Z, reshape(imag(Etmatrix(idx, :)), size(X)));
%   colorbar();
%   subplot(1, 2, 2);
%   surf(X, Y, Z, reshape(imag(Edipole(idx, :)), size(X)));
%   colorbar();
end

function testFullInvAlpha(testCase)

  alpha = 2.0;
  Ndipoles = 3;
  target = repmat(diag([0.5, 0.5, 0.5]), 1, Ndipoles);

  inv_alpha = ott.TmatrixDda.alpha_to_full_inv_alpha(alpha, Ndipoles);
  testCase.verifyEqual(inv_alpha, target);
  
end

function testRotsymMatrix(testCase)

  F_total = randn(100, 100);
  m = 1;
  z_rotation_safe = 1;
  [Feven, Fodd] = ott.TmatrixDda.combine_rotsym_matrix(F_total, ...
            m, z_rotation_safe);
          
  testCase.verifyEqual(Feven, F_total);
  testCase.verifyEqual(Fodd, F_total);

end

function testNearfieldMatrixTotal(testCase)

  Ndipoles = 30;
  xyz = randn(Ndipoles, 3);
  rtp = ott.utils.xyz2rtp(xyz);
  k = 2*pi;
  r_near = 10;
  
  Npts = 20;
  tp = randn(2, Npts);

  F_total = ott.TmatrixDda.nearfield_matrix_total(...
            tp(1, :), tp(2, :), r_near, k, xyz, rtp, false, 1, ...
            @ott.TmatrixDda.dipole_nearfield);
          
  testCase.verifySize(F_total, [3*Npts, Ndipoles*3], 'sz');
  
  weights = randn(3*Ndipoles, 1);
  fields = F_total * weights;
  
  fieldXyz = ott.utils.rtp2xyz([r_near.*ones(size(tp, 2), 1), tp.']);
  M_dipole = eye(3);
  field_target = zeros(3*Npts, 1);
  for ii = 1:Ndipoles
    D = xyz(ii, :).';
    F = ott.TmatrixDda.dipole_nearfield(D, fieldXyz.', M_dipole, k);
    field_target = field_target + F * weights((1:3) + (ii-1)*3);
  end
  
  testCase.verifyEqual(fields, field_target, 'AbsTol', 1.0e-12, 'fields');

end

function testTwoDipolesDirect(testCase)
  % Test fields from multiple dipoles are outgoing
  
  D1 = [0;0;0];
  D2 = [0;0;1];
  k = 2*pi;
  rtp = ott.utils.xyz2rtp([D1, D2].');
  spacing = 0.1;
  n_relative = 1.2;
  alpha = ott.utils.polarizability.LDR(...
                spacing, n_relative);
  inv_alpha = ott.TmatrixDda.alpha_to_full_inv_alpha(alpha, size(rtp, 1));
  
  % Calculate interatction matrix
  A = ott.TmatrixDda.interaction_A_total(k, rtp, inv_alpha, false, 1);
  
  % Calculate fields going out
  
  offset = 0.5;
  distance1 = linspace(10, 11, 100).';
  distance2 = distance1 + offset;
  tp = [pi/2, pi/2] .* ones(100, 2);
  rtp1 = [distance1,tp];
  rtp2 = [distance2,tp];
  targets_xyz1 = ott.utils.rtp2xyz(rtp1).';
  targets_xyz2 = ott.utils.rtp2xyz(rtp2).';
  M_dipole = eye(3);
  k = 2*pi;
  
  F1 = ott.TmatrixDda.dipole_nearfield(D1, targets_xyz1, M_dipole, k);
  F1 = [F1, ott.TmatrixDda.dipole_nearfield(D2, targets_xyz1, M_dipole, k)];
  F2 = ott.TmatrixDda.dipole_nearfield(D1, targets_xyz2, M_dipole, k);
  F2 = [F2, ott.TmatrixDda.dipole_nearfield(D2, targets_xyz2, M_dipole, k)];
  
  % Calculate incident beam
  Ei = ott.utils.vswfcart(2, -2, rtp(:, 1)*k, rtp(:, 2), rtp(:, 3), 'regular');
  Ei = Ei.';
  
  % Solve DDA problem
  E1 = ott.TmatrixDda.solve_and_evaluate(A, F1, Ei(:), false);
  E2 = ott.TmatrixDda.solve_and_evaluate(A, F2, Ei(:), false);
  
  E1 = E1(1:3:end);
  E2 = E2(1:3:end);
  
%   figure();
%   plot(distance1, [real(E1), real(E1*exp(-1i*pi/8))]);
  
  testCase.verifyEqual(real(E1*exp(-1i*offset*k)), real(E2), ...
      'RelTol', 0.1, 'AbsTol', 1.0e-3);
  
end

function testInteractionA(testCase)

  Ndipoles = 5;
  inv_alpha = randn(3, 3*Ndipoles);
  xyz = randn(Ndipoles, 3);
  k = 2*pi;

  A = ott.TmatrixDda.interaction_A_total(k, xyz, inv_alpha, false, 1);
  
  testCase.verifySize(A, [3, 3].*Ndipoles, 'sz');
  testCase.verifyEqual(A(1:3, 1:3), inv_alpha(1:3, 1:3), 'first diag');
  
  Az = ott.TmatrixDda.interaction_A_total(k, xyz, 0*inv_alpha, false, 1);
  testCase.verifyEqual(Az, Az.', 'AbsTol', 1.0e-12, 'symmetric');
  
end