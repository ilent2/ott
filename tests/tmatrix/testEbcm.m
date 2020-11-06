function tests = testEbcm
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testFromShape(testCase)

  radius = 1.0;
  shape = ott.shape.Sphere(radius);
  index = 1.2;
  tmatrix = ott.tmatrix.Ebcm.FromShape(shape, 'relative_index', index);

  testCase.verifyEqual(tmatrix.relative_index, index, 'index');
  testCase.verifyEqual(tmatrix.points(1, :), ...
      repmat(shape.radius, [1, size(tmatrix.points, 2)]), ...
      'AbsTol', 1.0e-2, 'radius');
  testCase.verifyEqual(tmatrix.invMethod, 'forwardslash', 'invMethod');
  testCase.verifyEqual(tmatrix.xySymmetry, true, 'xySymmetry');

  % Compare to Mie result
  Tmie = ott.tmatrix.Mie(radius, 'relative_index', index);
  testCase.verifyEqual(tmatrix.data, Tmie.data, 'AbsTol', 1.0e-2, 'data');

end

function testInternalMie(testCase)

  radius = 1.0;
  shape = ott.shape.Sphere(radius);
  index = 1.2;
  Nmax = 12;
  [~, tmatrix] = ott.tmatrix.Ebcm.FromShape(shape, ...
      'relative_index', index, 'Nmax', Nmax);

  testCase.verifyEqual(tmatrix.relative_index, index, 'index');
  testCase.verifyEqual(tmatrix.invMethod, 'inv', 'invMethod');
  testCase.verifyEqual(tmatrix.xySymmetry, true, 'xySymmetry');

  % Compare to Mie result
  [~, Tmie] = ott.tmatrix.Mie(radius, 'relative_index', index, ...
      'Nmax', Nmax);
  testCase.verifyEqual(tmatrix.data, Tmie.data, 'AbsTol', 2.0e-2, 'data');

end

