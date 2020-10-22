function tests = testMieLayered
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testFromShape(testCase)

  shape = [ott.shape.Sphere(0.5), ott.shape.Sphere(1.0)];
  tmatrix = ott.tmatrix.MieLayered.FromShape(shape, [1.2, 1.1]);

  testCase.verifyEqual(tmatrix.radii, [0.5, 1.0]);
  testCase.verifyEqual(tmatrix.relative_indices, [1.2, 1.1]);

end

