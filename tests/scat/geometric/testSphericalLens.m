function tests = testSphericalLens
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstruct(testCase)

  position = [0;0;0];
  focal_length = 2.0;

  L = ott.scat.geometric.SphericalLens(focal_length, 'position', position);
  testCase.verifyEqual(L.focal_length, focal_length, 'fl');

end

