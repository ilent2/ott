function tests = testSuperellipsoid
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  radii = [1, 2, 3];
  ew = 1;
  ns = 1.4;
  shape = ott.shape.Superellipsoid(radii, ew, ns);

  testCase.verifyEqual(shape.position, [0;0;0], 'position');
  testCase.verifyEqual(shape.radii, radii(:), 'radii');
  testCase.verifyEqual(shape.ew, ew, 'ew');
  testCase.verifyEqual(shape.ns, ns, 'ns');
  testCase.verifyEqual(shape.xySymmetry, true, 'xy');
  testCase.verifyEqual(shape.zRotSymmetry, 2, 'zrot');
  testCase.verifyEqual(shape.isEllipsoid, false, 'isellipse');
  testCase.verifyEqual(shape.isSphere, false, 'issphere');

end

