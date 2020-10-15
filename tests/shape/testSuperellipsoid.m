function tests = testSuperellipsoid
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  shape = ott.shape.Superellipsoid();

  testCase.verifyEqual(shape.position, [0;0;0], 'position');

end

