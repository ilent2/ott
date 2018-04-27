function tests = shapesurface
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testValues(testCase)

  import ott.shapes.Shape;

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  tol = 1.0e-4;

  % Create a grid of test points
  [theta, phi] = ott.utils.angulargrid(3, 3);
  sz = numel(theta);

  % Test a sphere with radius 1
  shape = Shape.simple('sphere', 1.0);
  r = shape.radii(theta, phi);
  n = shape.normals(theta, phi);
  [~, ~, rotsym] = shape.axialSymmetry();
  testCase.verifyThat(r, IsEqualTo(repmat(1.0, sz, 1), ...
      'Within', AbsoluteTolerance(tol)), ...
      'Sphere radius should be 1.0');
  testCase.verifyThat(n, IsEqualTo(repmat([1.0, 0, 0], sz, 1), ...
      'Within', AbsoluteTolerance(tol)), ...
      'Sphere normals should be [1 0 0]');
  testCase.verifyThat(rotsym, IsEqualTo(0), ...
      'Sphere rotational symmetry incorrect');

  % Test a ellipsoid
  shape = Shape.simple('ellipsoid', [ 1, 2, 3 ]);
  r = shape.radii(theta, phi);
  n = shape.normals(theta, phi);
  [~, ~, rotsym] = shape.axialSymmetry();
  rTarget = [ 1.7321    1.0000    1.7321    2.2780    1.5119, ...
      2.2780    2.2780    1.5119    2.2780 ].';
  nTarget = [
    0.6547    0.7559         0
    1.0000    0.0000         0
    0.6547   -0.7559         0
    0.6670    0.4892    0.5620
    0.8030    0.0000    0.5960
    0.6670   -0.4892    0.5620
    0.6670    0.4892   -0.5620
    0.8030    0.0000   -0.5960
    0.6670   -0.4892   -0.5620 ];
  testCase.verifyThat(r, IsEqualTo(rTarget, ...
      'Within', AbsoluteTolerance(tol)), ...
      'Ellipsoid radii incorrect');
  testCase.verifyThat(n, IsEqualTo(nTarget, ...
      'Within', AbsoluteTolerance(tol)), ...
      'Ellipsoid normals incorrect');
  testCase.verifyThat(rotsym, IsEqualTo(2), ...
      'Ellipsoid rotational symmetry incorrect');

  % Test a cylinder
  shape = Shape.simple('cylinder', [ 1, 1 ]);
  r = shape.radii(theta, phi);
  n = shape.normals(theta, phi);
  [~, ~, rotsym] = shape.axialSymmetry();
  rTarget = [ 0.5774    1.0000    0.5774    0.5774    1.0000, ...
      0.5774    0.5774    1.0000    0.5774 ].';
  nTarget = [
    0.8660   -0.5000         0
    1.0000    0.0000         0
    0.8660    0.5000         0
    0.8660   -0.5000         0
    1.0000    0.0000         0
    0.8660    0.5000         0
    0.8660   -0.5000         0
    1.0000    0.0000         0
    0.8660    0.5000         0 ];
  testCase.verifyThat(r, IsEqualTo(rTarget, ...
      'Within', AbsoluteTolerance(tol)), ...
      'Cylinder radii incorrect');
  testCase.verifyThat(n, IsEqualTo(nTarget, ...
      'Within', AbsoluteTolerance(tol)), ...
      'Cylinder normals incorrect');
  testCase.verifyThat(rotsym, IsEqualTo(0), ...
      'Cylinder rotational symmetry incorrect');

  % Test a cylinder (x-axis aligned)
  % Doesn't work, not actually supported
  %[r, n, rotsym] = shapesurface(theta, phi, -1, [ 1, 1 ]);
  %disp(r);
  %disp(n);
  %disp(rotsym);

  % Test a superellipsoid
  shape = Shape.simple('superellipsoid', [ 0.2, 0.5, 1, 1, 1 ]);
  r = shape.radii(theta, phi);
  n = shape.normals(theta, phi);
  [~, ~, rotsym] = shape.axialSymmetry();
  rTarget = [0.3780    0.2000    0.3780    0.5714    0.3288, ...
      0.5714    0.5714    0.3288    0.5714 ].';
  nTarget = [
    0.5587    0.8294         0
    1.0000    0.0000         0
    0.5587   -0.8294         0
    0.4680    0.5460    0.6949
    0.7131    0.0000    0.7010
    0.4680   -0.5460    0.6949
    0.4680    0.5460   -0.6949
    0.7131    0.0000   -0.7010
    0.4680   -0.5460   -0.6949 ];
  testCase.verifyThat(r, IsEqualTo(rTarget, ...
      'Within', AbsoluteTolerance(tol)), ...
      'Superellipsoid radii incorrect');
  testCase.verifyThat(n, IsEqualTo(nTarget, ...
      'Within', AbsoluteTolerance(tol)), ...
      'Superellipsoid normals incorrect');
  testCase.verifyThat(rotsym, IsEqualTo(1), ...
      'Superellipsoid rotational symmetry incorrect');

  % Test a cone-tipped cylinder
  shape = Shape.simple('cone-tipped-cylinder', [ 1, 1, 1 ]);
  r = shape.radii(theta, phi);
  n = shape.normals(theta, phi);
  [~, ~, rotsym] = shape.axialSymmetry();
  rTarget = [1.0981    1.0000    1.0981    1.0981    1.0000, ...
      1.0981    1.0981    1.0000    1.0981 ].';
  nTarget = [
    0.9659    0.2588         0
    1.0000    0.0000         0
    0.9659   -0.2588         0
    0.9659    0.2588         0
    1.0000    0.0000         0
    0.9659   -0.2588         0
    0.9659    0.2588         0
    1.0000    0.0000         0
    0.9659   -0.2588         0 ];
  testCase.verifyThat(r, IsEqualTo(rTarget, ...
      'Within', AbsoluteTolerance(tol)), ...
      'Cone-tipped cylinder radii incorrect');
  testCase.verifyThat(n, IsEqualTo(nTarget, ...
      'Within', AbsoluteTolerance(tol)), ...
      'Cone-tipped cylinder normals incorrect');
  testCase.verifyThat(rotsym, IsEqualTo(0), ...
      'Cone-tipped cylinder rotational symmetry incorrect');

  % Test a cube
  shape = Shape.simple('cube', 1);
  r = shape.radii(theta, phi);
  n = shape.normals(theta, phi);
  [~, ~, rotsym] = shape.axialSymmetry();
  rTarget = [0.5774; 0.5000; 0.5774; 0.5774; 0.5774; ...
      0.5774; 0.5774; 0.5774; 0.5774];
  testCase.verifyThat(r, IsEqualTo(rTarget, ...
      'Within', AbsoluteTolerance(tol)), ...
      'Cube radii incorrect');
  nTarget = [
    0.8660   -0.5000         0
    1.0000    0.0000         0
    0.8660    0.5000         0
    0.8660   -0.5000         0
    0.8660    0.0000   -0.5000
    0.8660    0.5000         0
    0.8660   -0.5000         0
    0.8660    0.0000    0.5000
    0.8660    0.5000         0 ];
  testCase.verifyThat(n, IsEqualTo(nTarget, ...
      'Within', AbsoluteTolerance(tol)), ...
      'Cube normals incorrect');
  testCase.verifyThat(rotsym, IsEqualTo(4), ...
      'Cube rotational symmetry incorrect');

end
