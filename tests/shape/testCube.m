function tests = testCube
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  width = 1.0;

  shape = ott.shape.Cube(width);

  testCase.verifyEqual(shape.width, width, 'width');
  testCase.verifyEqual(shape.starShaped, true, 'star');
  testCase.verifyEqual(shape.volume, width.^3, 'volume');
  testCase.verifyEqual(shape.xySymmetry, true, 'xysym');
  testCase.verifyEqual(shape.zRotSymmetry, 4, 'zsym');
  testCase.verifyEqual(shape.maxRadius, sqrt(3*width.^2./4), ...
      'AbsTol', 1.0e-10, 'maxRadius');

end

function testIntersectInternal(testCase)

  width = 1;
  shape = ott.shape.Cube(width);
  
  t1 = shape.intersect([0;0;0], [2;0;0]);
  testCase.verifyEqual(t1, [0.5;0;0], 't1');
  
  t2 = shape.intersectInternal([0, 0; 0, 0; 0, 0], [2, 0; 0, 2; 0, 0]);
  testCase.verifyEqual(t2, [0.5, 0; 0, 0.5; 0, 0], 't2');

end

