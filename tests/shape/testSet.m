function tests = testSet
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  % Can't test directly, use Intersection

  a = ott.shape.Cube();
  b = ott.shape.Cylinder();

  shape = ott.shape.Intersection([a, b]);

  testCase.verifyEqual(shape.shapes, [a, b], 'shapes');

end

