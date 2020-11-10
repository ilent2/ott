function tests = testCylinder
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  radius = 1.0;
  height = 2.0;

  shape = ott.shape.Cylinder(radius, height);

  testCase.verifyEqual(shape.radius, radius, 'rad');
  testCase.verifyEqual(shape.height, height, 'height');
  testCase.verifyEqual(shape.maxRadius, ...
      sqrt(radius^2 + (height/2)^2), 'maxRadius');
  testCase.verifyEqual(shape.volume, pi*radius.^2*height, 'volume');
  testCase.verifyEqual(shape.zRotSymmetry, 0, 'zrot');
  testCase.verifyEqual(shape.xySymmetry, true, 'xysym');
  testCase.verifyEqual(shape.boundingBox, [-1, 1; -1, 1; -1, 1]);
  testCase.verifyEqual(shape.perimeter, 4*radius+2*height, 'perimeter');

  % Scale
  shape = shape / 2;
  testCase.verifyEqual(shape.radius, radius/2, 'scale radius');
  testCase.verifyEqual(shape.height, height/2, 'scale height');

end

function testCasts(testCase)

  shape = ott.shape.Cylinder();
  shape = ott.shape.AxisymInterp(shape);
  testCase.verifyInstanceOf(shape, 'ott.shape.AxisymInterp');

end

function testAxisymStarShapeFunctions(testCase)

  shape = ott.shape.Cylinder();
  shape.normalsRtInternal([0;0]);  % Coverage;
  shape.insideRtInternal([0;0]);  % Coverage

end

function testSurf(testCase)

  f = figure();
  testCase.addTeardown(@() close(f));
  shape = ott.shape.Cylinder();
  a = shape.surf('visualise', true, 'showNormals', true);

end

function testIntersect(testCase)

  radius = 1;
  height = 2;
  shape = ott.shape.Cylinder(radius, height);
  
  testCase.verifyEqual(shape.intersect([0;0;0], [0;0;10]), [0;0;1], 'top');
  testCase.verifyEqual(shape.intersect([0;0;-10], [0;0;0]), [0;0;-1], 'bot');
  testCase.verifyEqual(shape.intersect([0;0;0], [2;0;0]), [1;0;0], 'rad');

end


