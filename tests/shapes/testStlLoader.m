function tests = stlloader
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testLoad(testCase)

  % Load the shape
  shape = ott.shapes.StlLoader('cube.stl');

  import matlab.unittest.constraints.IsEqualTo;

  testCase.verifyThat(numel(shape.verts), IsEqualTo(3*8), ...
    'Incorrect number of vertices');

  testCase.verifyThat(numel(shape.faces), IsEqualTo(3*2*6), ...
    'Incorrect number of faces');

end
