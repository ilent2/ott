function tests = testFieldVector
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testCartesianToSpherical(testCase)

  loc = randn(3, 5);
  val = randn(3, 5);
  S = ott.utils.FieldVector(loc, val, 'cartesian');

  rtp = ott.utils.xyzv2rtpv(val, loc);
  
  testCase.verifyEqual(S.locations, loc, 'loc');
  testCase.verifyEqual(S.vxyz, val, 'xyz');
  testCase.verifyEqual(S.vrtp, rtp, 'xyz');
  
end

function testSphericalToCartesian(testCase)

  loc = randn(3, 5);
  val = randn(3, 5);
  S = ott.utils.FieldVector(loc, val, 'spherical');

  xyz = ott.utils.rtpv2xyzv(val, loc);
  
  testCase.verifyEqual(S.locations, loc, 'loc');
  testCase.verifyEqual(S.vrtp, val, 'xyz');
  testCase.verifyEqual(S.vxyz, xyz, 'xyz');
  
end
