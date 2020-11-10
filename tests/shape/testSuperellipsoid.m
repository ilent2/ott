function tests = testSuperellipsoid
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  radii = [1, 2, 3];
  ew = 1;
  ns = 1.4;
  shape = ott.shape.Superellipsoid(radii, ew, ns);

  testCase.verifyEqual(shape.position, [0;0;0], 'position');
  testCase.verifyEqual(shape.radii, radii(:), 'radii');
  testCase.verifyEqual(shape.ew, ew, 'ew');
  testCase.verifyEqual(shape.ns, ns, 'ns');
  testCase.verifyEqual(shape.xySymmetry, true, 'xy');
  testCase.verifyEqual(shape.zRotSymmetry, 2, 'zrot');
  testCase.verifyEqual(shape.isEllipsoid, false, 'isellipse');
  testCase.verifyEqual(shape.isSphere, false, 'issphere');
  
  % Test sphere
  sph = ott.shape.Sphere();
  shape = ott.shape.Superellipsoid(sph);
  testCase.verifyEqual(shape.boundingBox, sph.boundingBox, 'bb');
  testCase.verifyEqual(shape.maxRadius, sph.maxRadius, 'maxR');
  
  % Scale
  shape = shape / 2;
  testCase.verifyEqual(shape.radii, sph.radius/2*[1;1;1], 'scale');

end

function testSurfMethods(testCase)

  sph = ott.shape.Sphere();
  shape = ott.shape.Superellipsoid(sph);
  
  testCase.verifyEqual(shape.starRadii(0, 0), ...
      sph.starRadii(0, 0), 'starRadii');
  testCase.verifyEqual(shape.normalsRtp([1;0;0]), ...
      sph.normalsRtp([1;0;0]), 'normalsRtp');
end

function testCasts(testCase)

  sph = ott.shape.Sphere();
  shape = ott.shape.Superellipsoid(sph);
  
  sph1 = ott.shape.Sphere(shape);
  testCase.verifyInstanceOf(sph1, 'ott.shape.Sphere');
  
  ellips = ott.shape.Ellipsoid(shape);
  testCase.verifyInstanceOf(ellips, 'ott.shape.Ellipsoid');
end
