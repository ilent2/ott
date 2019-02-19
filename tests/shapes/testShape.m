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
