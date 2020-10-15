function tests = testIntersection
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  a = ott.shape.Cube();
  b = ott.shape.Cylinder();

  shape = ott.shape.Intersection([a, b]);

  testCase.verifyEqual(shape.shapes, [a, b], 'shapes');

end

