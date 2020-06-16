function tests = testAxisymInterp
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  points = [1, 1; 2, 3];

  shape = ott.shapes.AxisymInterp(points);

  testCase.verifyEqual(shape.points, points);

end

function testConeTippedCylinder(testCase)

  shape = ott.shapes.AxisymInterp.ConeTippedCylinder();
  testCase.verifyEqual(shape.xySymmetry, true, 'xysym');
  testCase.verifyEqual(shape.starShaped, true, 'starShaped');

end

function testBicone(testCase)

  shape = ott.shapes.AxisymInterp.Bicone();
  testCase.verifyEqual(shape.xySymmetry, true, 'xysym');
  testCase.verifyEqual(shape.starShaped, true, 'starShaped');

end

