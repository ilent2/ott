function tests = testEllipsoid
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  radii = [1; 2; 3];
  shape = ott.shapes.Ellipsoid(radii);

  testCase.verifyEqual(shape.radii, radii, 'radii');

  % test volume where sphere
  shape = ott.shapes.Ellipsoid([1, 1, 1]);
  testCase.verifyEqual(shape.volume, 4/3*pi, 'AbsTol', 1.0e-14, 'volume');

end

