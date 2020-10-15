function tests = testPatchMesh
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  cube = ott.shape.Cube();

  shape = ott.shape.PatchMesh(cube.verts, cube.faces);

  testCase.verifyEqual(shape.verts, cube.verts, 'verts');
  testCase.verifyEqual(shape.faces, cube.faces, 'faces');

end

