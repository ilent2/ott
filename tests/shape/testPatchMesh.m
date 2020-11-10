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
  testCase.verifyEqual(shape.volume, cube.volume, 'volume');
  testCase.verifyEqual(shape.maxRadius, cube.maxRadius, 'maxR');
  testCase.verifyEqual(shape.boundingBox, cube.boundingBox, 'bb');

  % Scale
  shape = shape ./ 2;
  cube = cube * 0.5;
  testCase.verifyEqual(shape.verts, cube.verts, 'scale');

end

function testInsideXyz(testCase)

  cube = ott.shape.Cube();
  shape = ott.shape.PatchMesh(cube);
  
  testCase.verifyEqual(shape.insideXyz([0;0;0]), true);
end
