function tests = shape
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testSimple(testCase)

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.IsOfClass;
  import matlab.unittest.constraints.RelativeTolerance;
  tol = 1.0e-3;

  r = 1.0;
  shape = ott.shapes.Shape.simple('sphere', r);
  testCase.verifyThat(shape, IsOfClass(?ott.shapes.Sphere), ...
      'Incorrect shape type');
  testCase.verifyThat(shape.perimiter, IsEqualTo(2*pi*r, ...
      'Within', RelativeTolerance(tol)), ...
      'Incorrect sphere radius');

end

function testInsideXyz(testCase)

  radius = 1.0;
  shape = ott.shapes.Shape.simple('sphere', radius);
  
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
