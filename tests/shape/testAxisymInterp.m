function tests = testAxisymInterp
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  points = [1, 1; 2, 3];

  shape = ott.shape.AxisymInterp(points);

  testCase.verifyEqual(shape.points, points);

end

function testInside(testCase)

  shape = ott.shape.AxisymInterp.ConeTippedCylinder();

  xyz = [0;0;0];
  testCase.verifyEqual(shape.insideXyz(xyz), true, 'xyz');

  rtp = [0;0;0];
  testCase.verifyEqual(shape.insideRtp(xyz), true, 'rtp');

end

function testNormals(testCase)

  height = 1; radius = 2; code = 3;
  shape = ott.shape.AxisymInterp.ConeTippedCylinder(height, radius, cone);

  rtp = [radius;pi/2;0];
  testCase.verifyEqual(shape.normalsRtp(rtp), [1;0;0], 'rtp');

  xyz = [1;0;0];
  testCase.verifyEqual(shape.normalsXyz(xyz), [1;0;0], 'xyz');

end

function testConeTippedCylinder(testCase)

  shape = ott.shape.AxisymInterp.ConeTippedCylinder();
  testCase.verifyEqual(shape.xySymmetry, true, 'xysym');
  testCase.verifyEqual(shape.starShaped, true, 'starShaped');

end

function testBicone(testCase)

  shape = ott.shape.AxisymInterp.Bicone();
  testCase.verifyEqual(shape.xySymmetry, true, 'xysym');
  testCase.verifyEqual(shape.starShaped, true, 'starShaped');

end

