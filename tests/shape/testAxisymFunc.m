function tests = testAxisymFunc
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  func = @(x) ones(size(x));
  type = 'angular';
  range = [0, 1];

  shape = ott.shape.AxisymFunc(func, type, 'range', range);

  testCase.verifyEqual(shape.func, func, 'func');
  testCase.verifyEqual(shape.type, type, 'type');
  testCase.verifyEqual(shape.range, range, 'range');
  testCase.verifyEqual(shape.starShaped, true, 'star');
  testCase.verifyEqual(shape.maxRadius, 1, 'maxR');
  
  shape = shape * 2;
  testCase.verifyEqual(shape.maxRadius, 2, 'maxR');

end

function testPill(testCase)

  shape = ott.shape.AxisymFunc.Pill();

  testCase.verifyEqual(shape.boundingBox, ...
    [-0.5, 0.5; -0.5, 0.5; -1, 1], 'AbsTol', 1.0e-6, 'bb');
  testCase.verifyEqual(shape.insideRtp([0;0;0]), true, 'inside');
  
  shape.normalsRtInternal([1;0]);  % Coverage

end

function testBiconcaveDisc(testCase)

  shape = ott.shape.AxisymFunc.BiconcaveDisc();

  testCase.verifyEqual(shape.type, 'axialSym', 'type');
  testCase.verifyEqual(shape.xySymmetry, true, 'xy-sym');
  testCase.verifyEqual(shape.zRotSymmetry, 0, 'z-sym');
  testCase.verifyEqual(shape.insideRtp([0;0;0]), true, 'inside');

end

function testSurf(testCase)

  shape = ott.shape.AxisymFunc.BiconcaveDisc();

  h = figure();
  testCase.addTeardown(@close, h);

  p = shape.surf('showNormals', true);
  testCase.verifyClass(p, 'matlab.graphics.primitive.Patch');
  axis equal;
  
  shape.surfPoints();

end


