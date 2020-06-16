function tests = testPatchMesh
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  cube = ott.shapes.Cube();

  shape = ott.shapes.PatchMesh(cube.verts, cube.faces);

  testCase.verifyEqual(shape.verts, cube.verts, 'verts');
  testCase.verifyEqual(shape.faces, cube.faces, 'faces');

end

