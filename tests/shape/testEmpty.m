function tests = testEmpty
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  shape = ott.shape.Empty();

  testCase.verifyEqual(shape.volume, 0, 'volume');
  testCase.verifyEqual(shape.xySymmetry, true, 'xysym');
  testCase.verifyEqual(shape.zRotSymmetry, 0, 'zsym');
  testCase.verifyEqual(shape.starShaped, true, 'star');
  testCase.verifyEqual(shape.maxRadius, 0, 'maxr');
  testCase.verifyEqual(shape.boundingBox, zeros(3, 2), 'bb');
  
  % Scale
  Sshape = shape * 1000;
  testCase.verifyEqual(Sshape, shape, 'scale');

end

function testMethods(testCase)

  shape = ott.shape.Empty();

  testCase.verifyEqual(shape.insideRtp([0;0;0]), false, 'insidertp');
  testCase.verifyEqual(shape.insideXyz([0;0;0]), false, 'insidexyz');
  testCase.verifyEqual(shape.normalsXyz([0;0;0]), nan(3, 1), 'insidexyz');
  testCase.verifyEqual(shape.normalsRtp([0;0;0]), nan(3, 1), 'insidertp');
  testCase.verifyEqual(shape.intersect([0;0;0], [1;0;0]), ...
    nan(3,1), 'isect');
  
  % Coverage
  shape.surf('visualise', false);
  shape.surfPoints();
  testCase.verifyEqual(shape.intersectAll([0;0;0], [1;0;0]), nan(3, 1), 'isectAll');

end

