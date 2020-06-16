function tests = testSet
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  % Can't test directly, use Union

  a = ott.shapes.Cube();
  b = ott.shapes.Cylinder();

  shape = ott.shapes.Intersection([a, b]);

  testCase.verifyEqual(shape.shapes, [a, b], 'shapes');

end

