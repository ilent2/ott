function tests = testObjLoader
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  filename = 'cube.obj';

  % Create an OBJ
  cube = ott.shape.Cube();
  cube = ott.shape.TriangularMesh(cube);
  cube.writeWavefrontObj(filename);

  shape = ott.shape.ObjLoader(filename);

  testCase.verifyEqual(shape.filename, filename, 'filename');
  testCase.verifyEqual(shape.verts, cube.verts, 'verts');
  testCase.verifyEqual(shape.faces, cube.faces, 'faces');
  testCase.verifyEqual(shape.maxRadius, cube.maxRadius, 'maxRadius');
  testCase.verifyEqual(shape.volume, cube.volume, 'volume');

  % Clean up file
  delete(filename);

end

