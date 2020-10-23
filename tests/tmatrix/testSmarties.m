function tests = testBsc
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testFromShape(testCase)

  shape = ott.shape.Sphere(1.0);
  index = 1.2;
  tmatrix = ott.tmatrix.Smarties.FromShape(shape, 'relative_index', index);
  Tmie = ott.tmatrix.Mie.FromShape(shape, 'relative_index', index);

  testCase.assertEqual(tmatrix.Nmax, Tmie.Nmax, 'Nmax');
  testCase.verifyEqual(tmatrix.data, Tmie.data, ...
      'AbsTol', 1.0e-6, 'data');

end

function testEllipsoid(testCase)

  radii = [1, 1, 0.8]*0.1;
  shape = ott.shape.Ellipsoid(radii);
  index = 1.2;
  Nmax = 5;
  tmatrix = ott.tmatrix.Smarties.FromShape(shape, ...
    'relative_index', index, 'Nmax', Nmax);
  Tpm = ott.tmatrix.Pointmatch.FromShape(shape, ...
    'relative_index', index, 'Nmax', Nmax);

  testCase.assertEqual(tmatrix.Nmax, Tpm.Nmax, 'Nmax');
  testCase.verifyEqual(tmatrix.data, Tpm.data, ...
      'AbsTol', 1.0e-3, 'data');

end

