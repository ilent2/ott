function tests = testTriangularMesh
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  verts = randn(3, 3);
  faces = [1; 2; 3];

  shape = ott.shape.TriangularMesh(verts, faces);

  testCase.verifyEqual(shape.verts, verts, 'verts');
  testCase.verifyEqual(shape.faces, faces, 'faces');

end

function testInsideXyz(testCase)

  cube = ott.shape.Cube();
  cube = ott.shape.TriangularMesh(cube);
  testCase.assertClass(cube, 'ott.shape.TriangularMesh', 'type');
  radius = cube.maxRadius;

  % Choose three points inside the shape and one outside
  btarget = [true, true, true, false].';
  x = [0.5*radius.*rand(1, 3), 4.0];
  y = [0.5*radius.*rand(1, 3), 4.0];
  z = [0.5*radius.*rand(1, 3), 4.0];

  xyz = [x(:), y(:), z(:)].';

  b = cube.insideXyz(xyz);
  testCase.verifyEqual(b, btarget, 'insideXyz');

end

