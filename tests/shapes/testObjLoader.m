function tests = testObjLoader
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  filename = 'cube.obj';

  % Create an OBJ
  cube = ott.shapes.Cube();
  cube = ott.shapes.TriangularMesh(cube);
  cube.writeWavefrontObj(filename);

  shape = ott.shapes.ObjLoader(filename);

  testCase.verifyEqual(shape.filename, filename, 'filename');
  testCase.verifyEqual(shape.verts, cube.verts, 'verts');
  testCase.verifyEqual(shape.faces, cube.faces, 'faces');
  testCase.verifyEqual(shape.maxRadius, cube.maxRadius, 'maxRadius');
  testCase.verifyEqual(shape.volume, cube.volume, 'volume');

  % Clean up file
  delete(filename);

end

