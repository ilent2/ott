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
  
  shape.maxRadius;  % Test coverage (not sure if its right)

  shape2 = a & b;
  testCase.verifyEqual(shape2, shape, '& op');
  
  shape3 = shape2 & a;
  testCase.verifyEqual(shape3.shapes(3), a, '& op interp');

end

function testMembers(testCase)

  cube = ott.shape.Cube();
  shape = ott.shape.Intersection(cube);
  
  testCase.verifyEqual(shape.boundingBox, cube.boundingBox, 'bb');
  testCase.verifyEqual(shape.insideRtp([0;0;0]), true, 'insideRtp');
  testCase.verifyEqual(shape.insideXyz([0;0;0]), true, 'insideXyz');

end