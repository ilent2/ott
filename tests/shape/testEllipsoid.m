function tests = testEllipsoid
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  radii = [1; 2; 3];
  shape = ott.shape.Ellipsoid(radii);

  testCase.verifyEqual(shape.radii, radii, 'radii');

  % test volume where sphere
  shape = ott.shape.Ellipsoid([1, 1, 1]);
  testCase.verifyEqual(shape.volume, 4/3*pi, 'AbsTol', 1.0e-14, 'volume');
  testCase.verifyEqual(shape.isSphere, true, 'isSphere');
  testCase.verifyEqual(shape.maxRadius, 1, 'maxRadius');
  testCase.verifyEqual(shape.xySymmetry, true, 'xy');
  testCase.verifyEqual(shape.zRotSymmetry, 0, 'xy');
  testCase.verifyEqual(shape.boundingBox, [-1,1;-1,1;-1,1], 'bb');
  testCase.verifyEqual(shape.starShaped, true, 'starShaped');
  
  % Scale
  shape = shape * 2;
  testCase.verifyEqual(shape.maxRadius, 2, 'scale');

end

function testCasts(testCase)

  sph = ott.shape.Sphere();
  elp = ott.shape.Ellipsoid(sph);
  testCase.assertEqual(elp.isSphere, true, 'isSphere');

  sph = ott.shape.Sphere(elp);
  testCase.verifyInstanceOf(sph, 'ott.shape.Sphere');

end

function testMethods(testCase)

  sph = ott.shape.Sphere();
  elp = ott.shape.Ellipsoid(sph);
  
  testCase.verifyEqual(elp.starRadii(0, 0), sph.starRadii(0, 0), 'starRadii');
  testCase.verifyEqual(elp.normalsRtp([1;0;0]), sph.normalsRtp([1;0;0]), 'normalsRtp');
  testCase.verifyEqual(elp.intersectAll([0;0;0], [2;0;0]), ...
      sph.intersectAll([0;0;0], [2;0;0]), 'intersectAll');
    
  testCase.verifyEqual(elp.starRadii(0, 0), sph.radius, 'starRadii');
  testCase.verifyEqual(elp.insideXyz([0;0;0]), true, 'inside');
  testCase.verifyEqual(elp.normalsXyz([1;0;0]), [1;0;0], ...
    'AbsTol', 1e-15, 'inside');
  

end
