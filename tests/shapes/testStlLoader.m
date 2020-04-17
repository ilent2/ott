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

function testInsideXyz(testCase)

  shape = ott.shapes.StlLoader('cube.stl');
  radius = shape.maxRadius;
  
  % Choose three points inside the shape and one outside
  b = [true, true, true, false].';
  x = [0.5*radius.*rand(1, 3), 4.0];
  y = [0.5*radius.*rand(1, 3), 4.0];
  z = [0.5*radius.*rand(1, 3), 4.0];
  
  xyz = [x(:), y(:), z(:)].';
  
  testCase.verifyEqual(shape.insideXyz(x, y, z), b, ...
    'insideXyz with 3 arguments failed');
  testCase.verifyEqual(shape.insideXyz(xyz), b, ...
    'insideXyz with 1 argument failed');
  
  testCase.verifyEqual(shape.insideXyz(x, y, z, 'origin', 'world'), b, ...
    'insideXyz with 3 arguments failed and optional arg');
  testCase.verifyEqual(shape.insideXyz(xyz, 'origin', 'world'), b, ...
    'insideXyz with 1 argument failed and optional arg');
end
