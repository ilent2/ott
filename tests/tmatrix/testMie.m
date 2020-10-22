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
  tmatrix = ott.tmatrix.Mie.FromShape(shape, 'relative_index', n);

  testCase.verifyEqual(tmatrix.relative_index, n, 'n');
  testCase.verifyEqual(tmatrix.radius, r, 'r');

end

function testShapeVolume(testCase)

  n = 1.2;
  r = 1.0;

  shape = ott.shape.Sphere(r);
  tmatrix = ott.tmatrix.Mie.ShapeVolume(shape, 'relative_index', n);

  testCase.verifyEqual(tmatrix.relative_index, n, 'n');
  testCase.verifyEqual(tmatrix.radius, r, 'r');

end


