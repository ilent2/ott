function tests = testPlane
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstruct(testCase)

  position = [0;0;0];
  focal_length = 2.0;

  L = ott.scat.geometric.ThinLens(focal_length, 'position', position);
  testCase.verifyEqual(L.focal_length, focal_length, 'fl');

end

function testScatter(testCase)

  position = [0;0;0];
  focal_length = 2.0;
  ray = ott.beam.Ray('direction', [0;0;1], 'origin', [0;0;-1]);

  L = ott.scat.geometric.ThinLens(focal_length, 'position', position);
  sray = L.scatter(ray);

  testCase.verifyClass(sray, 'ott.beam.ScatteredRay');

end
