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
  tmatrix = ott.tmatrix.Ebcm.FromShape(shape, 'index_relative', index);

  testCase.verifyEqual(tmatrix.index_relative, index, 'index');
  testCase.verifyEqual(tmatrix.points(1, :), ...
      repmat(shape.radius, [1, size(tmatrix.points, 2)]), ...
      'AbsTol', 1.0e-2, 'radius');
  testCase.verifyEqual(tmatrix.invMethod, 'forwardslash', 'invMethod');
  testCase.verifyEqual(tmatrix.xySymmetry, true, 'xySymmetry');

  % Compare to Mie result
  Tmie = ott.tmatrix.Mie(radius, 'index_relative', index);
  testCase.verifyEqual(tmatrix.data, Tmie.data, 'AbsTol', 1.0e-2, 'data');

end

function testInternalMie(testCase)

  radius = 1.0;
  shape = ott.shape.Sphere(radius);
  index = 1.2;
  Nmax = 11;
  [~, tmatrix] = ott.tmatrix.Ebcm.FromShape(shape, ...
      'index_relative', index, 'Nmax', Nmax);

  testCase.verifyEqual(tmatrix.index_relative, index, 'index');
  testCase.verifyEqual(tmatrix.invMethod, 'inv', 'invMethod');
  testCase.verifyEqual(tmatrix.xySymmetry, true, 'xySymmetry');

  % Compare to Mie result
  [~, Tmie] = ott.tmatrix.Mie(radius, 'index_relative', index, ...
      'Nmax', Nmax);
  testCase.verifyEqual(tmatrix.data, Tmie.data, 'AbsTol', 2.0e-2, 'data');

end

