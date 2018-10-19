function tests = tmatrixpm
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)

  % Ensure the ott package is in our path
  addpath('../');

  % Target T-matrix
  testCase.TestData.Tmie = ott.TmatrixMie(1.0, 'index_relative', 1.2);

  % Tolerance for comparisons
  testCase.TestData.tol = 1.0e-6;

end

function testSimpleSphere(testCase)
  % This tests with both mirror and axial symmetry optimisations

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  Tmie = testCase.TestData.Tmie;
  tol = testCase.TestData.tol;

  Tpm = ott.TmatrixPm.simple('sphere', 1.0, 'index_relative', 1.2);

  testCase.verifyThat(Tpm.Nmax, IsEqualTo(Tmie.Nmax), ...
      'Nmax does not match Mie T-matrix');

  testCase.verifyThat(Tpm.data, IsEqualTo(Tmie.data, ...
      'Within', AbsoluteTolerance(tol)), ...
      'T-matrix does not match Mie T-matrix within tolerance');

end

function testSphereNoSym(testCase)

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  Tmie = testCase.TestData.Tmie;
  tol = testCase.TestData.tol;

  % Create a sphere T-matrix without rotational or mirror symmetry
  shape = ott.shapes.Shape.simple('sphere', 1.0);
  Nmax = Tmie.Nmax;
  z_rotational_symmetry = 1;
  z_mirror_symmetry = false;

  % Get the coordinates of the shape
  rtp = shape.angulargrid(max(Nmax), 'full', true);
  normals = shape.normals(rtp(:, 2), rtp(:, 3));

  Tpm = ott.TmatrixPm(rtp, normals, ...
      'index_relative', 1.2, ...
      'Nmax', Nmax, ...
      'z_mirror_symmetry', z_mirror_symmetry, ...
      'z_rotational_symmetry', z_rotational_symmetry);

  testCase.verifyThat(Tpm.Nmax, IsEqualTo(Tmie.Nmax), ...
      'Nmax does not match Mie T-matrix');

  testCase.verifyThat(Tpm.data, IsEqualTo(full(Tmie.data), ...
      'Within', AbsoluteTolerance(tol)), ...
      'T-matrix does not match Mie T-matrix within tolerance');

end

function testSphereMirrorSym(testCase)

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  Tmie = testCase.TestData.Tmie;
  tol = testCase.TestData.tol;

  % Create a sphere T-matrix without rotational or mirror symmetry
  shape = ott.shapes.Shape.simple('sphere', 1.0);
  Nmax = Tmie.Nmax;
  z_rotational_symmetry = 1;
  z_mirror_symmetry = true;

  % Get the coordinates of the shape
  rtp = shape.angulargrid(max(Nmax), 'full', true);
  rtp = rtp(rtp(:, 2) < pi/2, :);
  normals = shape.normals(rtp(:, 2), rtp(:, 3));

  Tpm = ott.TmatrixPm(rtp, normals, ...
      'index_relative', 1.2, ...
      'Nmax', Nmax, ...
      'z_mirror_symmetry', z_mirror_symmetry, ...
      'z_rotational_symmetry', z_rotational_symmetry);

  testCase.verifyThat(Tpm.Nmax, IsEqualTo(Tmie.Nmax), ...
      'Nmax does not match Mie T-matrix');

  testCase.verifyThat(Tpm.data, IsEqualTo(Tmie.data, ...
      'Within', AbsoluteTolerance(tol)), ...
      'T-matrix does not match Mie T-matrix within tolerance');

end

function testSphereRot4Sym(testCase)

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  Tmie = testCase.TestData.Tmie;
  tol = testCase.TestData.tol;

  % Create a sphere T-matrix without rotational or mirror symmetry
  shape = ott.shapes.Shape.simple('sphere', 1.0);
  Nmax = Tmie.Nmax;
  z_rotational_symmetry = 4;
  z_mirror_symmetry = false;

  % Get the coordinates of the shape
  rtp = shape.angulargrid(max(Nmax), 'full', true);
  rtp = rtp(rtp(:, 3) < pi/2, :);
  normals = shape.normals(rtp(:, 2), rtp(:, 3));

  Tpm = ott.TmatrixPm(rtp, normals, ...
      'index_relative', 1.2, ...
      'Nmax', Nmax, ...
      'z_mirror_symmetry', z_mirror_symmetry, ...
      'z_rotational_symmetry', z_rotational_symmetry);

  testCase.verifyThat(Tpm.Nmax, IsEqualTo(Tmie.Nmax), ...
      'Nmax does not match Mie T-matrix');

  testCase.verifyThat(Tpm.data, IsEqualTo(Tmie.data, ...
      'Within', AbsoluteTolerance(tol)), ...
      'T-matrix does not match Mie T-matrix within tolerance');

  %% Mirror symmetry
  z_mirror_symmetry = true;

  rtp = rtp(rtp(:, 2) < pi/2, :);
  normals = shape.normals(rtp(:, 2), rtp(:, 3));

  Tpm = ott.TmatrixPm(rtp, normals, ...
      'index_relative', 1.2, ...
      'Nmax', Nmax, ...
      'z_mirror_symmetry', z_mirror_symmetry, ...
      'z_rotational_symmetry', z_rotational_symmetry);

  testCase.verifyThat(Tpm.Nmax, IsEqualTo(Tmie.Nmax), ...
      'Nmax does not match Mie T-matrix (mirror)');

  testCase.verifyThat(Tpm.data, IsEqualTo(Tmie.data, ...
      'Within', AbsoluteTolerance(tol)), ...
      'T-matrix does not match Mie T-matrix within tolerance (mirror)');

end

function testCube(testCase)

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;

  % Only seems to work to about 1% agreement (good enough)
  tol = 1.0e-1;

  % Turn off warnings
  warning('off', 'MATLAB:rankDeficientMatrix');

  shape = ott.shapes.Shape.simple('cube', 1.0);
  Nmax = 15;
  aNmax = 20;

  % 4th order rotational symmetry and no mirror
  z_rotational_symmetry = 4;
  z_mirror_symmetry = false;

  rtp = shape.angulargrid(aNmax, 'full', true);
  rtp = rtp(rtp(:, 3) < pi/2, :);
  normals = shape.normals(rtp(:, 2), rtp(:, 3));

  T1 = ott.TmatrixPm(rtp, normals, ...
      'index_relative', 1.2, ...
      'Nmax', Nmax, ...
      'z_mirror_symmetry', z_mirror_symmetry, ...
      'z_rotational_symmetry', z_rotational_symmetry);

  % 1st order rotational symmetry and no mirror
  z_rotational_symmetry = 1;
  z_mirror_symmetry = false;

  rtp = shape.angulargrid(aNmax, 'full', true);
  normals = shape.normals(rtp(:, 2), rtp(:, 3));

  T2 = ott.TmatrixPm(rtp, normals, ...
      'index_relative', 1.2, ...
      'Nmax', Nmax, ...
      'z_mirror_symmetry', z_mirror_symmetry, ...
      'z_rotational_symmetry', z_rotational_symmetry);

  % 4th order rotational symmetry and mirror
  z_rotational_symmetry = 4;
  z_mirror_symmetry = true;

  rtp = shape.angulargrid(aNmax, 'full', true);
  rtp = rtp(rtp(:, 2) < pi/2, :);
  rtp = rtp(rtp(:, 3) < pi/2, :);
  normals = shape.normals(rtp(:, 2), rtp(:, 3));

  T3 = ott.TmatrixPm(rtp, normals, ...
      'index_relative', 1.2, ...
      'Nmax', Nmax, ...
      'z_mirror_symmetry', z_mirror_symmetry, ...
      'z_rotational_symmetry', z_rotational_symmetry);

  % 1st order rotational symmetry and mirror
  z_rotational_symmetry = 1;
  z_mirror_symmetry = true;

  rtp = shape.angulargrid(aNmax, 'full', true);
  rtp = rtp(rtp(:, 2) < pi/2, :);
  normals = shape.normals(rtp(:, 2), rtp(:, 3));

  T4 = ott.TmatrixPm(rtp, normals, ...
      'index_relative', 1.2, ...
      'Nmax', Nmax, ...
      'z_mirror_symmetry', z_mirror_symmetry, ...
      'z_rotational_symmetry', z_rotational_symmetry);

  %% Verifications

  testCase.verifyThat(T1.Nmax, IsEqualTo([Nmax, Nmax]), ...
      'Nmax does not match (1)');
  testCase.verifyThat(T2.Nmax, IsEqualTo([Nmax, Nmax]), ...
      'Nmax does not match (2)');
  testCase.verifyThat(T3.Nmax, IsEqualTo([Nmax, Nmax]), ...
      'Nmax does not match (3)');
  testCase.verifyThat(T4.Nmax, IsEqualTo([Nmax, Nmax]), ...
      'Nmax does not match (4)');

  testCase.verifyThat(T2.data, IsEqualTo(full(T1.data), ...
      'Within', AbsoluteTolerance(tol)), ...
      'T-matrix does not match (2)');
  testCase.verifyThat(T3.data, IsEqualTo(T1.data, ...
      'Within', AbsoluteTolerance(tol)), ...
      'T-matrix does not match (3)');
  testCase.verifyThat(T4.data, IsEqualTo(T1.data, ...
      'Within', AbsoluteTolerance(tol)), ...
      'T-matrix does not match (4)');

end

function testInternal(testCase)
  % Test construction of internal T-matrix
  
  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;

  Tpm = ott.TmatrixPm.simple('sphere', 1.0, ...
      'index_relative', 1.2, 'internal', true);

  Tmie = ott.TmatrixMie(1.0, 'index_relative', 1.2, 'internal', true);

  testCase.verifyThat(Tpm.Nmax, IsEqualTo(Tmie.Nmax), ...
      'Nmax does not match Mie T-matrix');

  tol = 1.0e-2;
  Tebcm_data = Tpm.data ./ max(abs(Tpm.data(:)));
  Tmie_data = Tmie.data ./ max(abs(Tmie.data(:)));
  testCase.verifyThat(Tebcm_data, IsEqualTo(Tmie_data, ...
      'Within', AbsoluteTolerance(tol)), ...
      'T-matrix does not match Mie T-matrix within tolerance');
    
end
