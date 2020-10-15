function tests = testSphere
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  radius = 1.5;
  shape = ott.shape.Sphere(radius);

  testCase.verifyEqual(shape.radius, radius, 'radius');
  testCase.verifyEqual(shape.volume, 4/3*pi*radius.^3, 'AbsTol', 1.0e-14, 'vol');
  testCase.verifyEqual(shape.perimeter, 2*pi*radius, 'perimeter');

end

