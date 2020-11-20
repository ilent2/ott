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
  
  % Set volume
  shape.volume = 2*2*2;
  testCase.verifyEqual(shape.width, 2, 'set volume');
  
  % Set maxRadius
  shape.maxRadius = sqrt(3*0.5^2);
  testCase.verifyEqual(shape.width, 1, 'RelTol', 1e-14, 'set maxR');

end

function testIntersectInternal(testCase)

  width = 1;
  shape = ott.shape.Cube(width);
  
  t1 = shape.intersect([0;0;0], [2;0;0]);
  testCase.verifyEqual(t1, [0.5;0;0], 't1');
  
  t2 = shape.intersectInternal([0, 0; 0, 0; 0, 0], [2, 0; 0, 2; 0, 0]);
  testCase.verifyEqual(t2, [0.5, 0; 0, 0.5; 0, 0], 't2');
  
  t3 = shape.intersectAll([0;0;0], [0;0;2]);
  testCase.verifyEqual(t3, [[0;0;0.5], nan(3, 1)], 't3');

end

function testNormals(testCase)

  width = 1;
  shape = ott.shape.Cube(width);
  
  point = [width/2;0;0];
  target = point ./ vecnorm(point);
  
  testCase.verifyEqual(shape.normalsXyz(point), target);

end


