function tests = testInverse
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  a = ott.shape.Cube();

  shape = ott.shape.Inverse(a);

  testCase.verifyEqual(shape.internal, a, 'shape');
  testCase.verifyEqual(shape.volume, Inf, 'vol');
  testCase.verifyEqual(shape.maxRadius, Inf, 'R');
  testCase.verifyEqual(shape.starShaped, a.starShaped, 'star');
  testCase.verifyEqual(shape.zRotSymmetry, a.zRotSymmetry', 'z');
  testCase.verifyEqual(shape.xySymmetry, a.xySymmetry', 'xy');
  testCase.verifyEqual(shape.boundingBox, Inf(3, 2).*[-1,1], 'bb');
  
  % Scale
  shape = shape * 2;
  testCase.verifyEqual(shape.internal.width, a.width*2, 'scale');

end

function testInverseInverse(testCase)

  a = ott.shape.Cube();
  b = ~a;
  c = ~b;

  testCase.verifyEqual(c, a, 'cube');
end

function testInside(testCase)

  a = ott.shape.Cube();
  shape = ott.shape.Inverse(a);
  
  testCase.verifyEqual(shape.insideXyz([0;0;0]), false, 'xyz');
  testCase.verifyEqual(shape.insideRtp([0;0;0]), false, 'rtp');

end

function testNormals(testCase)

  a = ott.shape.Cube(1);
  shape = ott.shape.Inverse(a);
  
  p = [0.5;0;0];
  rtpp = ott.utils.xyz2rtp(p);
  t = [-1;0;0];
  
  testCase.verifyEqual(shape.normalsXyz(p), t, 'xyz');
  testCase.verifyEqual(shape.normalsRtp(rtpp), t, 'rtp');

end

function testIntersect(testCase)

  a = ott.shape.Cube(1);
  shape = ott.shape.Inverse(a);
  
  testCase.verifyEqual(shape.intersect([0;0;0], [0;0;10]), ...
      [0;0;0.5]);
  testCase.verifyEqual(shape.intersectAll([0;0;-10], [0;0;10]), ...
      [0,0;0,0;-0.5,0.5]);
end

function testSurf(testCase)

  a = ott.shape.Cube(1);
  shape = ott.shape.Inverse(a);
  
  shape.surf('visualise', false);
  shape.surfPoints();

end

function testSetMethods(testCase)

  a = ott.shape.Cube();
  a = ott.shape.PatchMesh(a);
  shape = ott.shape.Inverse(a);
  
  % Test set methods (coverage)
  shape.zRotSymmetry = 0;
  shape.xySymmetry = true;
  shape.starShaped = true;

end
