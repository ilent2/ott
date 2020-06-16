function tests = testStlLoader
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  filename = 'cube.stl';

  shape = ott.shapes.StlLoader(filename);

  testCase.verifyEqual(shape.filename, filename, 'filename');
  testCase.verifyEqual(size(shape.verts), [3, 8], 'sz verts');
  testCase.verifyEqual(size(shape.faces), [3, 2*6], 'sz faces');

end

