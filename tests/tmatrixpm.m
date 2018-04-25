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
