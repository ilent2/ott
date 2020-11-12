function tests = testMie
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testFromShape(testCase)

  n = 1.2;
  r = 1.0;

  shape = ott.shape.Sphere(r);
  tmatrix = ott.tmatrix.Mie.FromShape(shape, 'index_relative', n);

  testCase.verifyEqual(tmatrix.index_relative, n, 'n');
  testCase.verifyEqual(tmatrix.radius, r, 'r');

end

function testShapeVolume(testCase)

  n = 1.2;
  r = 1.0;

  shape = ott.shape.Sphere(r);
  tmatrix = ott.tmatrix.Mie.ShapeVolume(shape, 'index_relative', n);

  testCase.verifyEqual(tmatrix.index_relative, n, 'n');
  testCase.verifyEqual(tmatrix.radius, r, 'r');

end

function testNmax1(testCase)

  tmatrix1 = ott.tmatrix.Mie(1, 'index_relative', 2, 'Nmax', 1);
  tmatrix2 = ott.tmatrix.Mie(1, 'index_relative', 2, 'Nmax', 2);
  
  D2 = diag(tmatrix2);
  
  testCase.verifyEqual(diag(tmatrix1), D2([1:3, 9:11]));

end

function testValuesFromOtt15(testCase)

  % Compare to values from OTTv1.5
  v15data = [-0.9074 + 0.2898i, -0.8633 + 0.3435i, ...
    -0.8665 + 0.3401i, -0.8390 + 0.3675i];
  v15data = diag([[1, 1, 1].*v15data(1), [1, 1, 1, 1, 1].*v15data(2), ...
                  [1, 1, 1].*v15data(3), [1, 1, 1, 1, 1].*v15data(4)]);
  
  tmatrix = ott.tmatrix.Mie(1.0, 'index_relative', 1.2, 'Nmax', 2);
  
  testCase.verifyEqual(full(tmatrix.data), v15data, ...
      'AbsTol', 1e-4, 'v15data');
                
end
